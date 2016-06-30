#include <deque>
#include <vector>
#include <bitset>
#include <iostream>
#include <algorithm>

#include <tclap/CmdLine.h>
#include <tbb/parallel_sort.h>

#include <dnachar.h>
#include <streamfastaparser.h>
#include <junctionapi/junctionapi.h>

class ChrReader
{
public:
	ChrReader(const std::vector<std::string> & fileName) : currentFile_(0), fileName_(fileName)
	{
		if (fileName.size() > 0)
		{
			parser_.reset(new TwoPaCo::StreamFastaParser(fileName[0]));
		}
	}

	bool NextChr(std::string & buf)
	{
		buf.clear();
		while (currentFile_ < fileName_.size())
		{
			if (parser_->ReadRecord())
			{
				char ch;
				while (parser_->GetChar(ch))
				{
					buf.push_back(ch);
				}

				return true;
			}
			else
			{
				if (++currentFile_ < fileName_.size())
				{
					parser_.reset(new TwoPaCo::StreamFastaParser(fileName_[currentFile_]));
				}
			}
		}

		return false;
	}

private:
	size_t currentFile_;
	std::vector<std::string> fileName_;
	std::auto_ptr<TwoPaCo::StreamFastaParser> parser_;
};

bool CompareJunctionsById(const TwoPaCo::JunctionPosition & a, const TwoPaCo::JunctionPosition & b)
{
	return a.GetId() < b.GetId();
}

bool CompareJunctionsByPos(const TwoPaCo::JunctionPosition & a, const TwoPaCo::JunctionPosition & b)
{
	return std::make_pair(a.GetChr(), a.GetPos()) < std::make_pair(b.GetChr(), b.GetPos());
}

struct EqClass
{
	int64_t label;
	std::vector<TwoPaCo::JunctionPosition> position;
};

int64_t Abs(int64_t x)
{
	return x > 0 ? x : -x;
}

int64_t CanonicalBeginId(TwoPaCo::JunctionPosition begin, TwoPaCo::JunctionPosition end)
{
	return begin.GetId() > 0 ? begin.GetId() : -end.GetId();
}

int64_t CanonicalEndId(TwoPaCo::JunctionPosition begin, TwoPaCo::JunctionPosition end)
{
	return begin.GetId() > 0 ? end.GetId() : -begin.GetId();
}

int64_t SegmentId(TwoPaCo::JunctionPosition begin, TwoPaCo::JunctionPosition end, char ch)
{
	int64_t beginId = CanonicalBeginId(begin, end);
	int64_t endId = CanonicalEndId(begin, end);
	if (Abs(beginId) >= INT32_MAX || Abs(endId) >= INT32_MAX)
	{
		throw std::runtime_error("A vertex id is too large, cannot generate GFA");
	}

	int64_t ret = (beginId << 2) || TwoPaCo::DnaChar::MakeUpChar(ch);
	return beginId == begin.GetId() ? ret : -ret;
}

bool CompareJunctionClasses(const EqClass & a, const EqClass & b)
{
	return CompareJunctionsByPos(a.position[0], b.position[0]);
}

void GenerateGroupOutupt(const std::string & inputFileName)
{
	TwoPaCo::JunctionPosition pos;
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::vector<EqClass> eqClass;
	std::vector<TwoPaCo::JunctionPosition> junction;
	while (reader.NextJunctionPosition(pos))
	{
		junction.push_back(pos);
	}

	std::sort(junction.begin(), junction.end(), CompareJunctionsById);
	for (size_t i = 0; i < junction.size();)
	{
		size_t j = i;
		for (; j < junction.size() && junction[i].GetId() == junction[j].GetId(); j++);
		std::sort(junction.begin() + i, junction.begin() + j, CompareJunctionsByPos);
		eqClass.push_back(EqClass());
		eqClass.back().label = junction[i].GetId();
		for (size_t k = i; k < j; k++)
		{
			eqClass.back().position.push_back(junction[k]);
		}

		i = j;
	}

	tbb::parallel_sort(eqClass.begin(), eqClass.end(), CompareJunctionClasses);
	for (auto junctionClass : eqClass)
	{
		for (auto j : junctionClass.position)
		{
			std::cout << j.GetChr() << ' ' << j.GetPos() << "; ";
		}

		std::cout << std::endl;
	}

}

void GenerateOrdinaryOutput(const std::string & inputFileName)
{
	TwoPaCo::JunctionPosition pos;
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	while (reader.NextJunctionPosition(pos))
	{
		std::cout << pos.GetChr() << ' ' << pos.GetPos() << ' ' << pos.GetId() << std::endl;
	}

}


void GenerateGfaOutupt(const std::string & inputFileName, size_t k)
{		
	bool start = true;
	std::deque<char> edge;
	int64_t prevSegmentId = 0;
	TwoPaCo::JunctionPosition end;
	TwoPaCo::JunctionPosition begin;
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::bitset<int64_t(1) << int64_t(32)> seen;
	seen.reset();
	int64_t previousId = 0;
	std::cout << "H\tVN:Z:1.0" << std::endl;
	while (reader.NextJunctionPosition(end))
	{
		if (!start)
		{
			if (begin.GetChr() == end.GetChr())
			{
				int64_t segmentId = SegmentId(begin, end, 'A');
				if (!seen[segmentId])
				{
					std::cout << "S\t" << segmentId << "\t" << "*" << std::endl;
					seen[segmentId] = true;
				}

			}
			else
			{

			}
			
			begin = end;
		}
		else
		{
			start = false;
			begin = end;
		}
	}
}

int main(int argc, char * argv[])
{
	try
	{
		TCLAP::CmdLine cmd("This utility converts the binary format into human readable one", ' ', "0");
		TCLAP::SwitchArg group("g", "group", "Group together positions of the same junctions", cmd, false);
		TCLAP::SwitchArg gfa("f", "gfa", "Output the grap in GFA format", cmd, false);

		TCLAP::UnlabeledValueArg<std::string> inputFileName("infile",
			"input file name",
			true,
			"",
			"file name",
			cmd);

		TCLAP::ValueArg<std::string> seqFileName("s",
			"seqfile",
			"sequences file name",
			false,
			"",
			"file name",
			cmd);

		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			true,
			25,
			"integer",
			cmd);

		cmd.parse(argc, argv);		

		if (group.isSet())
		{
			GenerateGroupOutupt(inputFileName.getValue());
		}
		else if (gfa.isSet())
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			GenerateGfaOutupt(inputFileName.getValue(), kvalue.getValue());
		}
		else
		{
			GenerateOrdinaryOutput(inputFileName.getValue());
		}
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
