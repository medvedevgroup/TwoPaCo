#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "common.h"
#include "compressedstring.h"


namespace TwoPaCo
{

	template<size_t CAPACITY>
	class BifurcationStorage
	{
	public:
		typedef CompressedString<CAPACITY> DnaString;
		BifurcationStorage(){}

		uint64_t GetDistinctVerticesCount() const
		{
			return bifurcationKey_.size();
		}

		uint64_t GetTotalVerticesCount() const
		{
			return (bifurcationKey_.size() - selfRevCompCount_) * 2 + selfRevCompCount_;
		}

		void Init(std::istream & bifurcationTempRead, uint64_t verticesCount, uint64_t vertexLength, size_t threads)
		{
			selfRevCompCount_ = 0;
			uint64_t bitsPower = 0;
			vertexLength_ = vertexLength;
			while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
			{
				++bitsPower;
			}

			size_t hashFunctionNumber = 3;
			bitsPower = max(bitsPower, size_t(24));
			bifurcationFilter_.assign(uint64_t(1) << bitsPower, false);
			hashFunction_.resize(hashFunctionNumber);
			for (HashFunctionPtr & ptr : hashFunction_)
			{
				ptr.reset(new HashFunction(vertexLength, bitsPower));
			}

			DnaString buf;
			std::string stringBuf(vertexLength, ' ');
			for (size_t i = 0; i < verticesCount; i++)
			{
				buf.ReadFromFile(bifurcationTempRead);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't read from a temporary file");
				}

				buf.ToString(stringBuf, vertexLength);
				if (DnaChar::IsSelfReverseCompliment(stringBuf.begin(), vertexLength))
				{
					selfRevCompCount_++;
				}

				bifurcationKey_.push_back(buf);
				for (HashFunctionPtr & ptr : hashFunction_)
				{
					uint64_t hf = ptr->hash(stringBuf);
					bifurcationFilter_[hf] = true;
				}
			}

			tbb::task_scheduler_init init(threads);
			tbb::parallel_sort(bifurcationKey_.begin(), bifurcationKey_.end(), DnaString::Less);
		}

		std::pair<uint64_t, uint64_t> GetId(std::string::const_iterator pos) const
		{
			return GetId(pos, true, true);
		}

		std::pair<uint64_t, uint64_t> GetId(std::string::const_iterator pos, const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash) const
		{
			bool posFound = true;
			bool negFound = true;
			uint64_t ret = INVALID_VERTEX;
			for (size_t i = 0; i < posVertexHash.size() && (posFound || negFound); i++)
			{
				if (!bifurcationFilter_[posVertexHash[i]->hashvalue])
				{
					posFound = false;
				}

				if (!bifurcationFilter_[negVertexHash[i]->hashvalue])
				{
					negFound = false;
				}
			}

			return GetId(pos, posFound, negFound);
		}

		const std::vector<HashFunctionPtr>& GetHashFunctions() const
		{
			return hashFunction_;
		}

	private:
		std::pair<uint64_t, uint64_t> GetId(std::string::const_iterator pos, bool posFound, bool negFound) const
		{
			DnaString bitBuf;
			uint64_t ret = INVALID_VERTEX;
			if (posFound)
			{
				posFound = false;
				bitBuf.Clear();
				bitBuf.CopyFromString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
				if (it != bifurcationKey_.end() && *it == bitBuf)
				{
					posFound = true;
					ret = it - bifurcationKey_.begin();
				}

			}

			if (negFound && !posFound)
			{
				negFound = false;
				bitBuf.Clear();
				bitBuf.CopyFromReverseString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
				if (it != bifurcationKey_.end() && *it == bitBuf)
				{
					negFound = true;
					ret = it - bifurcationKey_.begin();
				}
			}
#ifdef _DEBUG
			bool found = false;
			for (size_t strand = 0; strand < 2; ++strand)
			{
				bitBuf.Clear();
				if (strand == 0)
				{
					bitBuf.CopyFromString(pos, vertexLength_);
				}
				else
				{
					bitBuf.CopyFromReverseString(pos, vertexLength_);
				}

				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
				if (it != bifurcationKey_.end() && *it == bitBuf)
				{
					found = true;
				}
			}

			assert(found == (posFound || negFound));
#endif
			if (posFound || negFound)
			{
				if (DnaChar::IsSelfReverseCompliment(pos, vertexLength_))
				{
					return std::make_pair(ret, ret);
				}

				if (posFound)
				{
					return std::make_pair(ret, ret + bifurcationKey_.size());
				}

				return std::make_pair(ret + bifurcationKey_.size(), ret);
			}

			return std::make_pair(INVALID_VERTEX, INVALID_VERTEX);
		}

		size_t vertexLength_;
		uint64_t selfRevCompCount_;
		std::vector<bool> bifurcationFilter_;
		std::vector<DnaString> bifurcationKey_;
		std::vector<HashFunctionPtr> hashFunction_;
	};
}

#endif