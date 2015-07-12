#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

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
			{
				for (uint64_t body : bifurcation_)
				{
					DnaString str(vertexSize_, &body);
					*out++ = str.ToString();
				}
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

		struct CompressedString
		{
			
			uint64_t str[capacity];

			CompressedString() {}
			CompressedString(uint64_t init)
			{
				str[0] = init;
			}

			bool operator == (const CompressedString & other) const
			{
				std::equal(str, str + capacity, other.str);
			}

			bool operator != (const CompressedString & other) const
			{
				return !(*this == other);
			}

			uint64_t TranslateIdx(uint64_t & idx) const
			{
				uint64_t ret = idx >> 5;
				idx = idx & ((uint64_t(1) << uint64_t(5)) - 1);
				return ret;
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
		};

		

	private:
		size_t vertexSize_;
		

		

		std::vector<CompressedString> bifurcation_;
	};
}

#endif	