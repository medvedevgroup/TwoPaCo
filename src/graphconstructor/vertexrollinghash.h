#ifndef _VERTEX_ROLLING_HASH_H_
#define _VERTEX_ROLLING_HASH_H_

#include "common.h"
#include "concurrentbitvector.h"
#include "ngramhashing/cyclichash.h"

namespace TwoPaCo
{
	typedef CyclicHash<uint64_t> HashFunction;
	typedef std::unique_ptr<HashFunction> HashFunctionPtr;

	class VertexRollingHashSeed
	{
	public:
		VertexRollingHashSeed(const VertexRollingHashSeed & s)
		{
			hashFunction_.resize(s.HashFunctionsNumber());
			for (HashFunctionPtr & ptr : hashFunction_)
			{
				ptr = HashFunctionPtr(new HashFunction(s.VertexLength(), s.BitsNumber()));
			}
		}

		VertexRollingHashSeed(size_t numberOfFunctions, size_t vertexLength, size_t bits)
		{
			hashFunction_.resize(numberOfFunctions);
			for (HashFunctionPtr & ptr : hashFunction_)
			{
				ptr = HashFunctionPtr(new HashFunction(vertexLength, bits));
			}
		}

		size_t VertexLength() const
		{
			return hashFunction_[0]->n;
		}

		size_t BitsNumber() const
		{
			return hashFunction_[0]->wordsize;
		}

		size_t HashFunctionsNumber() const
		{
			return hashFunction_.size();
		}

	private:		
		std::vector<HashFunctionPtr> hashFunction_;
		friend class VertexRollingHash;
	};

	class VertexRollingHash
	{
	public:
		VertexRollingHash(const VertexRollingHashSeed & seed, std::string::const_iterator begin, size_t hashFunctions)
		{
			size_t size = seed.hashFunction_[0]->n;
			posVertexHash_.resize(hashFunctions);
			negVertexHash_.resize(hashFunctions);
			for (size_t i = 0; i < posVertexHash_.size(); i++)
			{
				posVertexHash_[i] = HashFunctionPtr(new HashFunction(*seed.hashFunction_[i]));
				negVertexHash_[i] = HashFunctionPtr(new HashFunction(*seed.hashFunction_[i]));
				for (auto it = begin; it != begin + size; ++it)
				{
					posVertexHash_[i]->eat(*it);
				}

				assert(posVertexHash_[i]->hashvalue == posVertexHash_[i]->hash(std::string(begin, begin + size)));
				for (std::string::const_reverse_iterator it(begin + size); it != std::string::const_reverse_iterator(begin); ++it)
				{
					char ch = DnaChar::ReverseChar(*it);
					negVertexHash_[i]->eat(DnaChar::ReverseChar(*it));
				}

				assert(negVertexHash_[i]->hashvalue == negVertexHash_[i]->hash(DnaChar::ReverseCompliment(std::string(begin, begin + size))));
			}
		}

		void Update(char positivePreviousChar, char positiveNextChar)
		{
			char negativeNextChar = DnaChar::ReverseChar(positiveNextChar);
			char negativePreviousChar = DnaChar::ReverseChar(positivePreviousChar);
			for (size_t i = 0; i < posVertexHash_.size(); i++)
			{
				posVertexHash_[i]->update(positivePreviousChar, positiveNextChar);
				negVertexHash_[i]->reverse_update(negativeNextChar, negativePreviousChar);
			}
		}

		bool Assert(std::string::const_iterator begin) const
		{
			size_t size = posVertexHash_[0]->n;
			for (size_t i = 0; i < posVertexHash_.size(); i++)
			{
				assert(posVertexHash_[i]->hashvalue == posVertexHash_[i]->hash(std::string(begin, begin + size)));
				assert(negVertexHash_[i]->hashvalue == negVertexHash_[i]->hash(DnaChar::ReverseCompliment(std::string(begin, begin + size))));
			}

			return true;
		}

		uint64_t RawPositiveHash(size_t hf) const
		{
			return posVertexHash_[hf]->hashvalue;
		}

		uint64_t RawNegativeHash(size_t hf) const
		{
			return negVertexHash_[hf]->hashvalue;
		}

		uint64_t GetVertexHash() const
		{			
			uint64_t posHash = posVertexHash_[0]->hashvalue;
			uint64_t negHash = negVertexHash_[0]->hashvalue;
			return min(posHash, negHash);
		}

		void GetIngoingEdgeHash(char previousPositiveCharacter, std::vector<uint64_t> & value) const
		{
			char previousNegativeCharacter = DnaChar::ReverseChar(previousPositiveCharacter);
			StrandComparisonResult result = DetermineStrandPrepend(previousPositiveCharacter, previousNegativeCharacter);
			if (result == positiveLess || result == tie)
			{
				GetPrependValues(previousPositiveCharacter, posVertexHash_, value);
			}
			else
			{
				GetExtendValues(previousNegativeCharacter, negVertexHash_, value);
			}
		}

		void GetOutgoingEdgeHash(char nextPositiveCharacter, std::vector<uint64_t> & value) const
		{
			char nextNegativeCharacter = DnaChar::ReverseChar(nextPositiveCharacter);
			
			StrandComparisonResult result = DetermineStrandExtend(nextPositiveCharacter, nextNegativeCharacter);
			if (result == positiveLess || result == tie)
			{
				GetExtendValues(nextPositiveCharacter, posVertexHash_, value);
			}
			else
			{
				GetPrependValues(nextNegativeCharacter, negVertexHash_, value);
			}
		}

	private:
		DISALLOW_COPY_AND_ASSIGN(VertexRollingHash);		
		std::vector<HashFunctionPtr> posVertexHash_;
		std::vector<HashFunctionPtr> negVertexHash_;

		enum StrandComparisonResult
		{
			positiveLess,
			negativeLess,
			tie
		};

		void GetPrependValues(char nextCh, const std::vector<HashFunctionPtr> & hashFunction, std::vector<uint64_t> & value) const
		{
			for (size_t i = 0; i < hashFunction.size(); i++)
			{
				value.push_back(hashFunction[i]->hash_prepend(nextCh));
			}
		}

		void GetExtendValues(char nextCh, const std::vector<HashFunctionPtr> & hashFunction, std::vector<uint64_t> & value) const
		{
			for (size_t i = 0; i < hashFunction.size(); i++)
			{
				value.push_back(hashFunction[i]->hash_extend(nextCh));
			}			
		}

		StrandComparisonResult DetermineStrandExtend(char nextCh, char revNextCh) const
		{
			for (size_t i = 0; i < posVertexHash_.size(); i++)
			{
				uint64_t posHash = posVertexHash_[i]->hash_extend(nextCh);
				uint64_t negHash = negVertexHash_[i]->hash_prepend(revNextCh);
				if (posHash != negHash)
				{
					return posHash < negHash ? positiveLess : negativeLess;
				}
			}

			return tie;
		}

		StrandComparisonResult DetermineStrandPrepend(char prevCh, char revPrevCh) const
		{
			for (size_t i = 0; i < posVertexHash_.size(); i++)
			{
				uint64_t posHash = posVertexHash_[i]->hash_prepend(prevCh);
				uint64_t negHash = negVertexHash_[i]->hash_extend(revPrevCh);
				if (posHash != negHash)
				{
					return posHash < negHash ? positiveLess : negativeLess;
				}
			}

			return tie;
		}
	};

	inline bool IsOutgoingEdgeInBloomFilter(const ConcurrentBitVector & filter, std::vector<uint64_t> & temp, const VertexRollingHash & hf, char farg)
	{
		temp.clear();
		hf.GetOutgoingEdgeHash(farg, temp);
		for (size_t i = 0; i < temp.size(); i++)
		{
			if (!filter.GetBit(temp[i]))
			{
				return false;
			}
		}

		return true;
	}

	inline bool IsIngoingEdgeInBloomFilter(const ConcurrentBitVector & filter, std::vector<uint64_t> & temp, const VertexRollingHash & hf, char farg)
	{
		temp.clear();
		hf.GetIngoingEdgeHash(farg, temp);
		for (size_t i = 0; i < temp.size(); i++)
		{
			if (!filter.GetBit(temp[i]))
			{
				return false;
			}
		}

		return true;
	}

}


#endif