#ifndef _CANDIDATE_OCCURENCE_
#define _CANDIDATE_OCCURENCE_

#include "compressedstring.h"

namespace Sibelia
{
	template<size_t CAPACITY>
	class CandidateOccurence
	{
	public:
		static const size_t ADDITIONAL_CHAR = 3;
		static const size_t MAX_SIZE = CAPACITY * 32;
		static const size_t NEXT_POS = MAX_SIZE - ADDITIONAL_CHAR;
		static const size_t PREV_POS = NEXT_POS + 1;
		static const size_t NMASK_POS = NEXT_POS + 2;
		static const size_t VERTEX_SIZE = MAX_SIZE - ADDITIONAL_CHAR;

		CandidateOccurence(){}
		void Set(uint64_t sequenceId,
			uint64_t position,
			uint64_t posHash0,
			uint64_t negHash0,			
			std::string::const_iterator pos,
			size_t vertexLength,
			char posExtend,
			char posPrev)
		{	
			position_ = position;
			sequenceId_ = sequenceId;
			if (posHash0 < negHash0 || (posHash0 == negHash0 && DnaChar::LessSelfReverseComplement(pos, vertexLength)))
			{
				body_.CopyFromString(pos, vertexLength);
				body_.SetChar(NEXT_POS, posExtend);
				body_.SetChar(PREV_POS, posPrev);
				body_.SetChar(NMASK_POS, EncodeNmask(posExtend, posPrev));
			}
			else
			{
				body_.CopyFromReverseString(pos, vertexLength);
				body_.SetChar(NEXT_POS, DnaChar::ReverseChar(posPrev));
				body_.SetChar(PREV_POS, DnaChar::ReverseChar(posExtend));
				body_.SetChar(NMASK_POS, EncodeNmask(posPrev, posExtend));
			}
		}

		char Prev() const
		{
			char nmask = body_.GetChar(NMASK_POS);
			return nmask & 2 ? 'N' : body_.GetChar(PREV_POS);
		}

		char Next() const
		{
			char nmask = body_.GetChar(NMASK_POS);
			return nmask & 1 ? 'N' : body_.GetChar(NEXT_POS);
		}

		uint32_t GetSequenceId() const
		{
			return sequenceId_;
		}

		uint32_t GetPosition() const
		{
			return position_;
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

		bool EqualBase(const CandidateOccurence & other) const
		{
			return CompressedString<CAPACITY>::EqualPrefix(VERTEX_SIZE, body_, other.body_);
		}

		bool IsSelfReverseCompliment(size_t vertexSize) const
		{
			for (size_t i = 0; i < vertexSize; i++)
			{
				if (body_.GetChar(i) != DnaChar::ReverseChar(body_.GetChar(vertexSize - i - 1)))
				{
					return false;
				}
			}

			return true;
		}
		
	private:
		char EncodeNmask(char next, char prev)
		{
			return DnaChar::LITERAL[(next == 'N' ? 0 : 1) | (prev == 'N' ? 0 : 2)];
		}

		uint32_t sequenceId_;
		uint32_t position_;
		CompressedString<CAPACITY> body_;
	};

	inline size_t CalculateNeededCapacity(size_t vertexLength)
	{
		size_t bufSize = vertexLength + CandidateOccurence<1>::ADDITIONAL_CHAR;
		return bufSize / UNIT_CAPACITY + (bufSize % UNIT_CAPACITY == 0 ? 0 : 1);
	}
}

#endif