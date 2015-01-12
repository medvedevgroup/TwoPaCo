#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#include <vector>

#include "streamfastaparser.h"

namespace Sibelia
{
	class VertexEnumerator
	{
	public:
		VertexEnumerator(const std::vector<std::string> & fileName, size_t k, size_t filterSize);
	private:
	};
}

#endif	