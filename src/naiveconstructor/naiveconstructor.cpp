#include <string> 
#include <vector>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <exception>
#include <tclap/CmdLine.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>

#include <streamfastaparser.h>

uint32_t currentUnknownChar = 4;

typedef std::vector<uint32_t> DnaString;

uint32_t MakeUpChar(char ch)
{
	if (!Sibelia::DnaChar::IsDefinite(ch))
	{
		if (currentUnknownChar == UINT32_MAX)
		{
			throw std::runtime_error("Too many unknown characters");
		}

		return currentUnknownChar++;
	}

	return Sibelia::DnaChar::MakeUpChar(ch);
}

uint32_t ReverseMadeUpChar(uint32_t ch)
{
	char before = Sibelia::DnaChar::UnMakeUpChar(ch);
	return MakeUpChar(Sibelia::DnaChar::ReverseChar(before));	
}

void MakeDeBruijnGraph(const std::vector<std::string> & fileName, std::ostream & outFile)
{
	std::vector<DnaString> posChar;
	std::vector<DnaString> negChar;
	for (auto name : fileName)
	{
		for (Sibelia::StreamFastaParser parser(name); parser.ReadRecord(); )
		{
			posChar.push_back(DnaString());
			negChar.push_back(DnaString());			
			for (char ch; !parser.GetChar(ch);)
			{
				posChar.back().push_back(MakeUpChar(ch));
			}

			for (auto it = posChar.back().rbegin(); it != posChar.back().rend(); ++it)
			{
				negChar.back().push_back(ReverseMadeUpChar(*it));
			}
		}
	}


}

int main(int argc, char * argv[])
{	
	return 0;
}