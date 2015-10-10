#ifndef _COMPRESSED_STRING_H_
#define _COMPRESSED_STRING_H_

#include <string>
#include <cstdint>
#include <algorithm>

#include "SpookyV2.h"

namespace Sibelia
{
	extern const size_t UNIT_CAPACITY;

	template<size_t CAPACITY>
	class CompressedString
	{
	public:
		CompressedString()
		{
			std::fill(str_, str_ + CAPACITY, 0);
		}

		size_t Bytes() const
		{
			return sizeof(str_[0] * CAPACITY);
		}

		CompressedString(const CompressedString & toCopy)
		{
			std::copy(toCopy.str_, toCopy.str_ + CAPACITY, str_);
		}

		static uint64_t Mask(size_t prefix)
		{
			return (uint64_t(1) << (prefix * 2)) - uint64_t(1);
		}

		static bool EqualPrefix(size_t prefix, const CompressedString & it1, const CompressedString & it2)
		{
			size_t remain = prefix;
			for (size_t i = 0; remain > 0; i++)
			{
				size_t current = std::min(remain, UNIT_CAPACITY);
				uint64_t apiece = it1.str_[i];
				uint64_t bpiece = it2.str_[i];
				if (current != UNIT_CAPACITY)
				{
					uint64_t mask = Mask(prefix);
					apiece &= mask;
					bpiece &= mask;
				}

				if (apiece != bpiece)
				{
					return false;
				}

				remain -= current;
			}

			return true;
		}

		static bool LessPrefix(const CompressedString & v1, const CompressedString & v2, size_t prefix)
		{
			size_t remain = prefix;
			for (size_t i = 0; remain > 0; i++)
			{
				size_t current = std::min(remain, UNIT_CAPACITY);
				uint64_t apiece = v1.str_[i];
				uint64_t bpiece = v2.str_[i];
				if (current != UNIT_CAPACITY)
				{
					uint64_t mask = Mask(prefix);
					apiece &= mask;
					bpiece &= mask;
				}

				if (apiece != bpiece)
				{
					return apiece < bpiece;
				}

				remain -= current;
			}

			return false;
		}

		void CopyPrefixFrom(const CompressedString & copy, size_t prefix)
		{
			size_t remain = prefix;
			for (size_t i = 0; remain > 0; i++)
			{
				uint64_t piece = copy.str_[i];
				size_t current = std::min(remain, UNIT_CAPACITY);
				if (current != UNIT_CAPACITY)
				{
					piece &= Mask(current);
				}

				str_[i] = piece;
				remain -= current;
			}
		}

		bool operator == (const CompressedString & other) const
		{
			return std::equal(str_, str_ + CAPACITY, other.str_);
		}

		bool operator != (const CompressedString & other) const
		{
			return !(*this == other);
		}

		static char Id(char ch)
		{
			return ch;
		}

		uint64_t HashPrefix(size_t prefix) const
		{
			CompressedString buf;
			buf.CopyPrefixFrom(*this, prefix);
			return SpookyHash::Hash64(buf.str_, sizeof(str_[0]) * CAPACITY, 0);
		}

		CompressedString ReverseComplement(size_t stringSize) const
		{
			CompressedString ret;
			for (size_t i = 0; i < stringSize; i++)
			{
				ret.SetChar(i, DnaChar::ReverseChar(GetChar(stringSize - i - 1)));
			}

			return ret;
		}

		void SetChar(uint64_t idx, char ch)
		{
			uint64_t element = TranslateIdx(idx);
			uint64_t charIdx = str_[element] >> (2 * idx);
			uint64_t mask = uint64_t(0x3) << (idx * 2);
			str_[element] &= ~(mask);
			str_[element] |= DnaChar::MakeUpChar(ch) << (2 * idx++);
		}

		char GetChar(uint64_t idx) const
		{
			uint64_t element = TranslateIdx(idx);
			uint64_t charIdx = str_[element] >> (2 * idx);
			return DnaChar::LITERAL[charIdx & 0x3];
		}

		char RawChar(uint64_t idx) const
		{
			uint64_t element = TranslateIdx(idx);
			return str_[element] >> (2 * idx);			
		}

		std::string ToString(size_t size) const
		{
			std::string ret(size, ' ');
			for (size_t i = 0; i < size; i++)
			{
				ret[i] = GetChar(i);
			}

			return ret;
		}

		void CopyFromString(std::string::const_iterator it, size_t size)
		{
			StrCpy(it, 0, 0, size, Id);
		}

		void CopyFromReverseString(std::string::const_iterator it, size_t size)
		{
			StrCpy(std::string::const_reverse_iterator(it + size), 0, 0, size, DnaChar::ReverseChar);
		}

	private:
		uint64_t str_[CAPACITY];

		template<class T, class F>
		void StrCpy(T src, size_t element, size_t idx, size_t size, F f)
		{
			for (size_t i = 0; i < size; i++)
			{
				str_[element] |= DnaChar::MakeUpChar(f(*src++)) << (2 * idx++);
				if (idx >= UNIT_CAPACITY)
				{
					idx = 0;
					++element;
				}
			}
		}

		uint64_t TranslateIdx(uint64_t & idx) const
		{
			uint64_t ret = idx >> 5;
			idx = idx & ((uint64_t(1) << uint64_t(5)) - 1);
			return ret;
		}		
	};

}


#endif