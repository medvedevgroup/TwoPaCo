#include <ctime>
#include <vector>
#include <cassert>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "test.h"
#include "vertexenumerator.h"

int main(int argc, char * argv[])

{
	assert(Sibelia::Runtests());
	
	try
	{		
		std::vector<std::string> fileName(argv + 2, argv + argc);
		std::stringstream ss(*(argv + 1));
		size_t filterSize;
		ss >> filterSize;
		Sibelia::VertexEnumerator vid(fileName, 25, filterSize);
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}
	

	return 0;
}