#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#include <algorithm>
#include <tbb/concurrent_vector.h>

#include "streamfastaparser.h"
#include "concurrentbitvector.h"

namespace Sibelia
{
	class VertexEnumerator
	{
	public:
		const static size_t INVALID_VERTEX;
		size_t GetVerticesCount() const;
		size_t GetId(const DnaString & vertex) const;
		
		template<class Iterator>
			void Dump(Iterator out)
			{/*
				for (CompressedString & str : bifurcation_)
				{
					*out++ = str.ToString();
				}*/
			}
			

		VertexEnumerator(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t aggregationThreads,
			const std::string & tmpFileName);

		static const size_t capacity = 1;

		class CompressedString
		{
		public:
			CompressedString()
			{
				std::fill(str, str + capacity, 0);
			}

			CompressedString(uint64_t init)
			{
				str[0] = init;
			}

			CompressedString(const CompressedString & toCopy)
			{
				std::copy(toCopy.str, toCopy.str + capacity, str);
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
					size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
					uint64_t apiece = it1.str[i];
					uint64_t bpiece = it2.str[i];
					if (current != DnaString::UNIT_CAPACITY)
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

			static bool LessPrefix(const VertexEnumerator::CompressedString & v1, const VertexEnumerator::CompressedString & v2, size_t prefix)
			{
				size_t remain = prefix;
				for (size_t i = 0; remain > 0; i++)
				{
					size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
					uint64_t apiece = v1.str[i];
					uint64_t bpiece = v2.str[i];
					if (current != DnaString::UNIT_CAPACITY)
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

			void StrCpyPrefix(const CompressedString & cpy, size_t prefix)
			{
				size_t remain = prefix;
				for (size_t i = 0; remain > 0; i++)
				{
					uint64_t piece = cpy.str[i];
					size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
					if (current != DnaString::UNIT_CAPACITY)
					{						
						piece &= Mask(current);
					}

					str[i] = piece;
					remain -= current;
				}
			}

			bool operator == (const CompressedString & other) const
			{
				return std::equal(str, str + capacity, other.str);
			}

			bool operator != (const CompressedString & other) const
			{
				return !(*this == other);
			}

			char GetChar(uint64_t idx) const
			{
				uint64_t element = TranslateIdx(idx);
				uint64_t charIdx = str[element] >> (2 * idx);
				return DnaString::LITERAL[charIdx & 0x3];
			}

			template<class T, class F>
			void StrCpy(T src, size_t & element, size_t & idx, size_t size, F f)
			{
				for (size_t i = 0; i < size; i++)
				{
					str[element] |= DnaString::MakeUp(f(*src++)) << (2 * idx++);
					if (idx >= DnaString::UNIT_CAPACITY)
					{
						idx = 0;
						++element;
					}
				}
			}

		private:
			uint64_t str[capacity];

			uint64_t TranslateIdx(uint64_t & idx) const
			{
				uint64_t ret = idx >> 5;
				idx = idx & ((uint64_t(1) << uint64_t(5)) - 1);
				return ret;
			}
		};

	private:
		size_t vertexSize_;
		std::vector<CompressedString> bifurcation_;
	};
}

#endif	