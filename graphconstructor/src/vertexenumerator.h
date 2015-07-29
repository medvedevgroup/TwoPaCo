#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#include <algorithm>
#include <tbb/concurrent_vector.h>

#include "compressedstring.h"
#include "streamfastaparser.h"
#include "concurrentbitvector.h"

namespace Sibelia
{
	std::string RevComp(const std::string & str);

	class VertexEnumerator
	{
	public:
		const static size_t INVALID_VERTEX;
		static const size_t capacity = 2;
		typedef CompressedString<capacity> DnaString;

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


	private:
		size_t vertexSize_;
		std::vector<DnaString> bifurcation_;
	};
}

#endif	