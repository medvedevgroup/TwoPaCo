#include <ctime>
#include <vector>
#include <cassert>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iostream>

#include "test.h"
#include "vertexenumerator.h"

int main(int argc, char * argv[])
{
	assert(Sibelia::Runtests());

	try
	{		
		std::vector<std::string> fileName(argv + 1, argv + argc);
		Sibelia::VertexEnumerator vid(fileName, 25, (1 << 28) + 1);
		std::cout << vid.GetVerticesCount() << std::endl;
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}
	

	return 0;
}