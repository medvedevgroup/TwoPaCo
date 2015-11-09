#include <set>
#include <map>
#include <deque>
#include <string> 
#include <vector>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <exception>
#include <tclap/CmdLine.h>

#include <streamfastaparser.h>

typedef uint32_t DnaChar;
DnaChar currentUnknownChar = 4;
typedef std::vector<DnaChar> DnaString;

struct DnaStringComparer
{
public:
	bool operator()(const DnaString & a, const DnaString  & b) const
	{
		return std::lexicographical_compare(a.begin() + 1, a.end() - 1, b.begin() + 1, b.end() - 1);
	}
};

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

DnaString SubStr(DnaString str, size_t pos, size_t k)
{
	return DnaString(str.begin() + pos, str.begin() + pos + k);
}

void InsertJunction(std::map<DnaString, uint64_t> & junctionMap, const DnaString & newJunction)
{
	if (junctionMap.count(newJunction) == 0)
	{
		size_t id = junctionMap.size();
		junctionMap[newJunction] = id;
	}
}

bool OccurenceEqual(const DnaString & a, const DnaString & b)
{
	return std::equal(a.begin() + 1, a.end() - 1, b.begin() + 1);
}

void MakeDeBruijnGraph(const std::vector<std::string> & fileName, size_t k, std::ostream & outFile)
{
	std::vector<DnaString> strand[2];
	for (auto name : fileName)
	{
		for (Sibelia::StreamFastaParser parser(name); parser.ReadRecord(); )
		{
			strand[0].push_back(DnaString());
			strand[1].push_back(DnaString());
			for (char ch; parser.GetChar(ch);)
			{
				strand[0].back().push_back(MakeUpChar(ch));
			}

			for (auto it = strand[0].back().rbegin(); it != strand[0].back().rend(); ++it)
			{
				strand[1].back().push_back(ReverseMadeUpChar(*it));
			}
		}
	}
	
	std::map<DnaString, uint64_t> junctionMap;	
	std::multiset<DnaString, DnaStringComparer> junctionOccurence;
	for (auto currentStrand : strand)
	{
		for (auto str : currentStrand)
		{
			if (str.size() >= k)
			{
				InsertJunction(junctionMap, SubStr(str, 0, k));
				InsertJunction(junctionMap, SubStr(str, str.size() - k, k));
			}

			for (size_t pos = 0; pos + k + 2 <= str.size(); ++pos)
			{
				junctionOccurence.insert(SubStr(str, pos, k + 2));
			}
		}
	}

	DnaStringComparer cmp;
	for (auto it = junctionOccurence.begin(); it != junctionOccurence.end(); )
	{
		auto jt = it;
		std::set<DnaChar> inGoing;
		std::set<DnaChar> outGoing;
		for (; jt != junctionOccurence.end() && OccurenceEqual(*it, *jt); ++jt)
		{
			inGoing.insert((*jt)[0]);
			outGoing.insert((*jt)[k + 1]);
		}

		if (inGoing.size() > 1 || outGoing.size() > 1)
		{
			InsertJunction(junctionMap, SubStr(*it, 1, k));
		}

		it = jt;
	}

	for (auto str : strand[0])
	{
		if (str.size() >= k)
		{
			for (size_t i = 0; i <= str.size() - k; i++)
			{
				DnaString kmer = SubStr(str, i, k);
				auto it = junctionMap.find(kmer);
				if (it != junctionMap.end())
				{
					outFile << i << ' ' << it->second << std::endl;
				}
			}
		}		
	}
}

int main(int argc, char * argv[])
{	
	try
	{
		TCLAP::CmdLine cmd("A really naive program for condensed de Bruijn graph construction", ' ', "0");

		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			false,
			25,
			"integer",
			cmd);

		TCLAP::UnlabeledMultiArg<std::string> fileName("filenames",
			"FASTA file(s) with nucleotide sequences.",
			true,
			"fasta files with genomes",
			cmd);

		TCLAP::ValueArg<std::string> outFileName("o",
			"outfile",
			"Output file name",
			false,
			"de_bruijn.bin",
			"file name",
			cmd);

		cmd.parse(argc, argv);

		std::ofstream out(outFileName.getValue().c_str());
		if (!out)
		{
			throw std::runtime_error("Can't open output file");
		}

		MakeDeBruijnGraph(fileName.getValue(), kvalue.getValue(), out);

	}
	catch (TCLAP::ArgException &e)
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