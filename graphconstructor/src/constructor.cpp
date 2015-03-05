#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <cassert>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sstream>

#include "test.h"
#include "vertexenumerator.h"

size_t Atoi(const char * str)
{
	size_t ret;
	std::stringstream ss(str);
	ss >> ret;
	return ret;
}



int main(int argc, char * argv[])
{
	/*
	CyclicHash<uint64>  hf(5, 64);
	hf.
	string input = "ABCDE";
	hf.eat(input[0]);//A
	hf.eat(input[1]);//B
	hf.eat(input[2]);//C
	hf.eat(input[3]);//D
	cout << "Hash value of ABCD is " << hf.hashvalue << endl;
	// we check the answer going the long way...
	const std::vector<unsigned char> charvectslice(input.begin(), input.begin() + 4);
	uint64_t trueanswerslice = hf.hash(charvectslice);
	if (trueanswerslice != hf.hashvalue) throw runtime_error("bug");
	// we continue
	hf.eat(input[4]);//E
	cout << "Hash value of ABCDE is " << hf.hashvalue << endl;
	// we check the answer going the long way
	const std::vector<unsigned char> charvect(input.begin(), input.end());
	uint64_t trueanswer = hf.hash(charvect);
	if (trueanswer != hf.hashvalue) throw runtime_error("bug");
	return 0;*/

	assert(Sibelia::Runtests());

	try
	{
		std::vector<std::string> fileName(argv + 4, argv + argc);

		size_t filterSize = Atoi(*(argv + 1));
		size_t hashFunctions = Atoi(*(argv + 3));
		Sibelia::VertexEnumerator vid(fileName, 25, filterSize, hashFunctions);
		std::cout << "Distinct = " << vid.GetVerticesCount();
		if (std::string(*(argv + 2)) == "1")
		{
			std::vector<std::string> all;
			vid.Dump(std::back_inserter(all));
			std::vector<std::string> rev;
			std::transform(all.begin(), all.end(), std::back_inserter(rev), static_cast<std::string(*)(const std::string&)>(Sibelia::DnaString::RevComp));
			all.insert(all.begin(), rev.begin(), rev.end());
			std::sort(all.begin(), all.end());
			all.erase(std::unique(all.begin(), all.end()), all.end());
			std::cout << "All = " << all.size() << std::endl;
		}
	}
	catch (const std::runtime_error & msg)
	{
		std::cerr << msg.what() << std::endl;
	}

	return 0;
}