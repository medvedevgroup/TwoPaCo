#ifndef _TEST_H_

#include "dnastring.h"
#include "vertexenumerator.h"

namespace Sibelia
{
	bool Runtests();
	void DnaStringTest(size_t n, std::ostream & log);
	void VertexEnumeratorTest(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, std::ostream & log);
}

#endif