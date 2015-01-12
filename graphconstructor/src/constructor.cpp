#include <ctime>
#include <vector>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vertexenumerator.h"

int main(int argc, char * argv[])
{

	try
	{
		std::vector<std::string> fileName;
		fileName.push_back(argv[1]);
		Sibelia::VertexEnumerator vid(fileName, 21, 1 << 28);
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}
	

	return 0;
}