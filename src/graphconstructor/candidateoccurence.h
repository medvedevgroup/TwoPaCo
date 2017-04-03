#ifndef _CANDIDATE_OCCURENCE_
#define _CANDIDATE_OCCURENCE_

#include "vertexrollinghash.h"
#include "compressedstring.h"

namespace TwoPaCo
{
	template<size_t CAPACITY>
	class CandidateOccurence
	{
	public:
		static const size_t ADDITIONAL_CHAR = 4;
		static const size_t MAX_SIZE = CAPACITY * 32;
		static const size_t MASK_POS = MAX_SIZE - ADDITIONAL_CHAR;
		static const size_t VERTEX_SIZE = MAX_SIZE - ADDITIONAL_CHAR;		

		CandidateOccurence(){}

		//EPIC GOVNOCODE STARTS
		bool Init(const VertexRollingHash & hash,
			std::string::const_iterator pos,
			size_t vertexLength,
			char prev,
			char next)
		{	
			bool ret = true;
			char setPrev = prev;
			char setNext = next;
			uint64_t posHash0 = hash.RawPositiveHash(0);
			uint64_t negHash0 = hash.RawNegativeHash(0);
			if (posHash0 < negHash0 || (posHash0 == negHash0 && DnaChar::LessSelfReverseComplement(pos, vertexLength)))
			{
				body_.CopyFromString(pos, vertexLength);
			}
			else
			{
				ret = false;
				body_.CopyFromReverseString(pos, vertexLength);
				setPrev = DnaChar::ReverseChar(next);
				setNext = DnaChar::ReverseChar(prev);
			}

			if (setPrev != 'N')
			{
				SetIngoingEdge(setPrev);
			}

			if (setNext != 'N')
			{
				SetOutgoingEdge(setNext);
			}

			return ret;
		}
		//EPIC GOVNOCODE ENDS		

		void Merge(const CandidateOccurence & other)
		{
			body_.OrInto(other.body_);
		}

		bool EqualBase(const CandidateOccurence & occurence) const
		{
			return CompressedString<CAPACITY>::EqualPrefix(VERTEX_SIZE, occurence.body_, body_);
		}

		uint64_t Hash() const
		{
			return body_.HashPrefix(VERTEX_SIZE);
		}

		CompressedString<CAPACITY> GetBase() const
		{
			CompressedString<CAPACITY> ret;
			ret.CopyPrefixFrom(body_, VERTEX_SIZE);
			return ret;
		}

		bool operator < (const CandidateOccurence & other) const
		{
			return CompressedString<CAPACITY>::LessPrefix(body_, other.body_, VERTEX_SIZE);
		}

		bool GetIngoingEdge(char ch) const
		{
			uint64_t charIdx = DnaChar::MakeUpChar(ch);
			uint64_t pos = MASK_POS * 2 + charIdx;
			return body_.GetBit(pos);
		}

		bool GetOutgoingEdge(char ch) const
		{
			uint64_t charIdx = DnaChar::MakeUpChar(ch);
			uint64_t pos = (MASK_POS + 2) * 2 + charIdx;
			return body_.GetBit(pos);
		}
		
	private:
		void SetIngoingEdge(char ch)
		{
			uint64_t charIdx = DnaChar::MakeUpChar(ch);
			uint64_t pos = MASK_POS * 2 + charIdx;
			body_.SetBit(pos);
			assert(GetIngoingEdge(ch));
		}

		void SetOutgoingEdge(char ch)
		{
			uint64_t charIdx = DnaChar::MakeUpChar(ch);
			uint64_t pos = (MASK_POS + 2) * 2 + charIdx;
			body_.SetBit(pos);
			assert(GetOutgoingEdge(ch));			
		}		

		CompressedString<CAPACITY> body_;
	};


	inline size_t CalculateNeededCapacity(size_t vertexLength)
	{
		size_t bufSize = vertexLength + CandidateOccurence<1>::ADDITIONAL_CHAR;
		return bufSize / UNIT_CAPACITY + (bufSize % UNIT_CAPACITY == 0 ? 0 : 1);
	}
}

#endif
