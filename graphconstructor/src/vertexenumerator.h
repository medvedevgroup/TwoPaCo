#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#include <vector>

#include "streamfastaparser.h"

namespace Sibelia
{
	class VertexEnumerator
	{
	public:
		const static size_t INVALID_VERTEX;
		size_t GetVerticesCount() const;
		size_t GetId(const DnaString & vertex) const;		
		VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize);
	private:
		
		std::vector<uint64_t> bifurcation_;
	};
}

#endif	