#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sstream>

#include <tclap/CmdLine.h>

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
	//assert(TwoPaCo::Runtests());
	try
	{
		TCLAP::CmdLine cmd("Program for construction of the condensed de Bruijn graph from complete genomes", ' ', "0.0.0");
		
		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			false,
			25,
			"integer",
			cmd);

		TCLAP::ValueArg<uint64_t> filterSize("f",
			"filtersize",
			"Size of the filter",
			true,
			0,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> hashFunctions("q",
			"hashfnumber",
			"Number of hash functions",
			false,
			5,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> rounds("r",
			"rounds",
			"Number of computation rounds",
			false,
			1,
			"integer",
			cmd);

		TCLAP::ValueArg<unsigned int> threads("t",
			"threads",
			"Number of worker threads",
			false,
			1,
			"integer",
			cmd);

		TCLAP::ValueArg<std::string> tmpDirName("",
			"tmpdir",
			"Temporary directory name",
			false,
			".",
			"directory name",
			cmd);

		TCLAP::UnlabeledMultiArg<std::string> fileName("filenames",
			"FASTA file(s) with nucleotide sequences.",
			true,
			"fasta files with genomes",
			cmd);

		TCLAP::ValueArg<std::string> outFileName("o",
			"outfile",
			"Output file name prefix",
			false,
			"de_bruijn",
			"file name",
			cmd);

		cmd.parse(argc, argv);
		
		std::unique_ptr<TwoPaCo::VertexEnumerator> vid = TwoPaCo::CreateEnumerator(fileName.getValue(),
			kvalue.getValue(), filterSize.getValue(),
			hashFunctions.getValue(),
			rounds.getValue(),
			threads.getValue(),
			tmpDirName.getValue(),
			outFileName.getValue());
			
		std::cout << "Distinct junctions = " << vid->GetVerticesCount() << std::endl;
		std::cout << std::endl;
	}
	catch (TCLAP::ArgException & e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (std::runtime_error & e)
	{
		std::cerr << "error: " << e.what() << std::endl;
		return 1;
	}
	
	return 0;
}