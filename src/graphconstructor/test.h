#ifndef _TEST_H_
#define _TEST_H_

#include <vector>
#include <string>
#include <iostream>

namespace TwoPaCo
{
	bool Runtests();
	void DnaStringTest(size_t n, std::ostream & log);
	void VertexEnumeratorTest(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, std::ostream & log);
}

#endif