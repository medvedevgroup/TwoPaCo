#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>

#include "vertexenumerator.h"

int main(int argc, char * argv[])
{
	try
	{
		Sibelia::StreamFastaParser parser("test.fasta");
		while (parser.ReadRecord())
		{
			std::cout << '>' << parser.GetCurrentHeader() << std::endl;
			for (char ch; parser.GetChar(ch); std::cout << ch);
			std::cout << std::endl;
		}
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}
	

	return 0;
}