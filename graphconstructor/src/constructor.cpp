#include <ctime>
#include <vector>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vertexenumerator.h"

int main(int argc, char * argv[])
{
	/*
	size_t k = 28;
	std::string test = "CCTGACTGACTGACCTCCTGACTGACTGAGATTTACCTGACTGACTGACCTATG";
	Sibelia::DnaString str;
	for (size_t i = 0; i < k; i++)
	{
		str.AppendBack(test[i]);
	}

	for (size_t i = 0; i <= test.size() - k; i++)
	{
		std::string trueKmer = std::string(test.begin() + i, test.begin() + i + k);
		std::cout << trueKmer << "\t" << str.ToString() << std::endl;
		assert(trueKmer == str.ToString());
		str.PopFront();
		std::cout << "Pop\t" << str.ToString() << std::endl;		
		str.AppendBack(test[i + k]);
		std::cout << "Append\t" << str.ToString() << std::endl;
	}*/


	try
	{
		std::vector<std::string> fileName;
		fileName.push_back("test.fasta");
		Sibelia::VertexEnumerator vid(fileName, 4, 1 << 28);
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}
	

	return 0;
}