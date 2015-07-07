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
		/*
		template<class Iterator>
			void Dump(Iterator out)
			{
				for (uint64_t body : bifurcation_)
				{
					DnaString str(vertexSize_, body);
					*out++ = str.ToString();
				}
			}*/
			

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
		std::vector<uint64_t> bifurcation_;
	};
}

#endif	