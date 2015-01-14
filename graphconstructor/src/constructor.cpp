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
		
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}
	

	return 0;
}