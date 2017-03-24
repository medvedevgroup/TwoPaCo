#include <map>
#include <deque>
#include <string>
#include <vector>
#include <bitset>
#include <memory>
#include <cassert>
#include <sstream>
#include <iterator>
#include <iostream>
#include <algorithm>

#include <tclap/CmdLine.h>
#include <tbb/parallel_sort.h>

#include <dnachar.h>
#include <streamfastaparser.h>
#include <junctionapi/junctionapi.h>


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

const int64_t ID_POWER = 35;
int64_t reservedPath = int64_t(1) << (ID_POWER - 1);
const int64_t MAX_JUNCTION_ID = int64_t(1) << (ID_POWER - 4);
const int64_t MAX_SEGMENT_NUMBER = int64_t(1) << ID_POWER;

class Segment
{
public:
	Segment() {}
	Segment(TwoPaCo::JunctionPosition begin, TwoPaCo::JunctionPosition end, char posEdgeCh, char negEdgeCh)
	{
		bool uniquePath = false;
		int64_t absBeginId = Abs(begin.GetId());
		int64_t absEndId = Abs(end.GetId());
		if (absBeginId >= MAX_JUNCTION_ID || absEndId >= MAX_JUNCTION_ID)
		{
			throw std::runtime_error("A vertex id is too large, cannot generate GFA");
		}

		if (absBeginId < absEndId || (absBeginId == absEndId && absBeginId > 0))
		{
			uniquePath = posEdgeCh == 'N';
			segmentId_ = TwoPaCo::DnaChar::MakeUpChar(posEdgeCh);
			begin_ = begin;
			end_ = end;
		}
		else
		{
			uniquePath = negEdgeCh == 'N';
			segmentId_ = TwoPaCo::DnaChar::MakeUpChar(negEdgeCh);
			begin_ = TwoPaCo::JunctionPosition(begin.GetChr(), begin.GetPos(), -end.GetId());
			end_ = TwoPaCo::JunctionPosition(end.GetChr(), end.GetPos(), -begin.GetId());
		}
		
		if (!uniquePath)
		{
			if (begin_.GetId() < 0)
			{
				segmentId_ |= 1 << 2;
				segmentId_ |= Abs(begin_.GetId()) << 3;
			}
			else
			{
				segmentId_ |= begin_.GetId() << 3;
			}

			if (begin.GetId() != begin_.GetId())
			{
				segmentId_ = -segmentId_;
			}
		}
		else
		{
			segmentId_ = reservedPath++;
		}
	}

	int64_t GetSegmentId() const
	{
		return segmentId_;
	}	

	int64_t GetAbsSegmentId() const
	{
		return Abs(segmentId_);
	}

private:
	int64_t segmentId_;
	TwoPaCo::JunctionPosition begin_;
	TwoPaCo::JunctionPosition end_;	
};

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

char Sign(int64_t arg)
{
	return arg >= 0 ? '+' : '-';
}



void ReadInputSequences(const std::vector<std::string> & genomes, std::vector<std::string> & chrSegmentId, std::vector<uint64_t> & chrSegmentLength, std::map<std::string, std::string> & fileName, bool noPrefix)
{
	size_t chrCount = 0;
	chrSegmentId.clear();	
	chrSegmentLength.clear();
	for (const std::string & chrFileName : genomes)
	{
		TwoPaCo::StreamFastaParser parser(chrFileName);
		while (parser.ReadRecord())
		{
			std::stringstream ssId;
			if (noPrefix)
			{
				ssId << parser.GetCurrentHeader();
			}
			else
			{
				ssId << "s" << chrCount << "_" << parser.GetCurrentHeader();
			}

			chrSegmentId.push_back(ssId.str());
			fileName[ssId.str()] = chrFileName;
			/*
			if (!silent)
			{
				std::cout << "S\t" << ssId.str() << "\t*\tUR:Z:" << chrFileName << std::endl;
			}			

			*/

			uint64_t size = 0;
			for (char ch; parser.GetChar(ch); ++size);
			chrSegmentLength.push_back(size);
		}
	}
}

class Gfa1Generator
{
public:
	void Header(std::ostream & out) const
	{
		out << "H\tVN:Z:1.0";
	}

	void ListInputSequences(const std::vector<std::string> & seq, std::map<std::string, std::string> & fileName, std::ostream & out) const
	{
		for (const auto & it : seq)
		{
			out << "S\t"
				<< it
				<< "\t*\tUR:Z:"
				<< fileName[it]
				<< std::endl;
		}
	}

	void Segment(int64_t segmentId, uint64_t segmentSize, const std::string & body, std::ostream & out) const
	{
		out << "S\t"
			<< Abs(segmentId) << "\t" 
			<< body << std::endl;
	}

	void Occurrence(int64_t segmentId, uint64_t segmentSize, const std::string & chrSegmentId, uint64_t chrSegmentSize, uint64_t begin, uint64_t end, uint64_t k, std::ostream & out) const
	{
		out << "C\t" 
			<< Abs(segmentId) << '\t' 
			<< Sign(segmentId) << '\t'
			<< chrSegmentId << "\t+\t" 
			<< end << std::endl;
	}

	void Edge(int64_t prevSegmentId, uint64_t prevSegmentSize, int64_t segmentId, uint64_t segmentSize, uint64_t k, std::ostream & out) const
	{
		out << "L\t" 
			<< Abs(prevSegmentId) << '\t' 
			<< Sign(prevSegmentId) << '\t' 
			<< Abs(segmentId) << '\t' 
			<< Sign(segmentId) << '\t' 
			<< k << 'M' << std::endl;
	}

	void FlushPath(std::vector<int64_t> & currentPath, const std::string & seqId, size_t k, std::ostream & out) const
	{
		if (currentPath.size() > 0)
		{
			out << "P\t" << seqId << '\t';
			std::copy(currentPath.begin(), currentPath.end() - 1, std::ostream_iterator<int64_t>(out, ","));
			out << currentPath.back();
			out << '\t';

			if (currentPath.size() > 1)
			{
				for (size_t i = 0; i < currentPath.size() - 2; i++)
				{
					out << k << "M,";
				}

				out << k << "M";
			}

			out << std::endl;
			currentPath.clear();
		}
	}
};

std::string Gfa2Position(size_t pos, size_t length)
{
	std::stringstream ss;
	if (pos == length)
	{
		ss << pos << "$";
	}
	else
	{
		ss << pos;
	}

	return ss.str();
}

std::string Gfa2Segment(int64_t segment)
{
	std::stringstream ss;
	ss << Abs(segment) << Sign(segment);
	return ss.str();
}

class Gfa2Generator
{
public:
	void Header(std::ostream & out) const
	{
		out << "H\tVN:Z:2.0" << std::endl;
	}

	void ListInputSequences(const std::vector<std::string> & seq, std::map<std::string, std::string> & fileName, std::ostream & out) const
	{

	}

	void Segment(int64_t segmentId, uint64_t segmentSize, const std::string & body, std::ostream & out) const
	{
		out << "S\t" 
			<< Abs(segmentId) << "\t" 
			<< segmentSize << "\t" 
			<< body << std::endl;
	}

	void Occurrence(int64_t segmentId, uint64_t segmentSize, const std::string & chrSegmentId, uint64_t chrSegmentSize, uint64_t begin, uint64_t end, uint64_t k, std::ostream & out) const
	{
		std::cout << "F\t"
			<< Abs(segmentId) << '\t'
			<< chrSegmentId << Sign(segmentId)
			<< "0\t"
			<< segmentSize << "$" << "\t"
			<< Gfa2Position(begin, chrSegmentSize)
			<< Gfa2Position(end, chrSegmentSize) << std::endl;
	}

	void Edge(int64_t prevSegmentId, uint64_t prevSegmentSize, int64_t segmentId, uint64_t segmentSize, uint64_t k, std::ostream & out) const
	{
		out << "E\t" 
			<< Gfa2Segment(prevSegmentId)
			<< '\t' << Gfa2Segment(segmentId) << '\t' 
			<< Gfa2Position(prevSegmentSize - k, prevSegmentSize) << '\t'
			<< Gfa2Position(prevSegmentSize, prevSegmentSize) << '\t' 
			<< Gfa2Position(segmentSize - k, segmentSize) << '\t'
			<< Gfa2Position(segmentSize, segmentSize) << '\t'
			<< k << 'M' << std::endl;
	}

	void FlushPath(std::vector<int64_t> & currentPath, const std::string & seqId, size_t k, std::ostream & out) const
	{
		if (currentPath.size() > 0)
		{
			out << "O\t" << seqId << "p" << '\t';
			std::copy(currentPath.begin(), currentPath.end(), std::ostream_iterator<int64_t>(out, " "));
			out << std::endl;
			currentPath.clear();
		}
	}
};

template<class G>
void GenerateGfaOutput(const std::string & inputFileName, const std::vector<std::string> & genomes, size_t k, bool prefix, G g)
{	
	std::vector<uint64_t> chrSegmentLength;
	std::vector<std::string> chrSegmentId;
	std::map<std::string, std::string> chrFileName;

	//std::cout << "H\tVN:Z:1.0" << std::endl;
	g.Header(std::cout);

	ReadInputSequences(genomes, chrSegmentId, chrSegmentLength, chrFileName, !prefix);
	g.ListInputSequences(chrSegmentId, chrFileName, std::cout);

	std::vector<int64_t> currentPath;
	const int64_t NO_SEGMENT = 0;
	std::string chr;	
	int64_t seqId = NO_SEGMENT;
	int64_t prevSegmentId = NO_SEGMENT;
	int64_t prevSegmentSize = -1;
	TwoPaCo::JunctionPosition end;
	TwoPaCo::JunctionPosition begin;
	TwoPaCo::ChrReader chrReader(genomes);
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::vector<bool> seen(MAX_SEGMENT_NUMBER, 0);
	int64_t previousId = 0;
	
#ifdef _DEBUG
	std::map<int64_t, std::string> segmentBody;
#endif
	if (reader.NextJunctionPosition(begin))
	{
		chrReader.NextChr(chr);
		while (reader.NextJunctionPosition(end))
		{
			if (begin.GetChr() == end.GetChr())
			{
				Segment nowSegment(begin, end, chr[begin.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]));				
				int64_t segmentId = nowSegment.GetSegmentId();
				currentPath.push_back(segmentId);
				uint64_t segmentSize = end.GetPos() + k - begin.GetPos();
				if (!seen[Abs(segmentId)])
				{
					//std::cout << "S\t" << Abs(segmentId) << "\t";
					std::stringstream ss;					
					if (segmentId > 0)
					{
						std::copy(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k, std::ostream_iterator<char>(ss));

					}
					else
					{
						std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
						std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(ss));
					}

					g.Segment(segmentId, chrSegmentLength[seqId], ss.str(), std::cout);
					seen[Abs(segmentId)] = true;
				}

#ifdef _DEBUG
				int64_t absSegmentId = Abs(segmentId);
				std::string buf = segmentId > 0 ? std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k) : 
					TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));				
				if (segmentBody.count(absSegmentId) == 0)
				{
					segmentBody[absSegmentId] = buf;
				}
				else
				{
					assert(segmentBody[absSegmentId] == buf);
				}
#endif
				g.Occurrence(segmentId, segmentSize, chrSegmentId[seqId], chrSegmentLength[seqId], begin.GetPos(), end.GetPos(), k, std::cout);
				//std::cout << "C\t" << Abs(segmentId) << '\t' << Sign(segmentId) << '\t' << chrSegmentId[seqId] << "\t+\t" << begin.GetPos() << std::endl;

				if (prevSegmentId != NO_SEGMENT)
				{
					//std::cout << "L\t" << Abs(prevSegmentId) << '\t' << Sign(prevSegmentId) << '\t' << Abs(segmentId) << '\t' << Sign(segmentId) << '\t' << k << 'M' << std::endl;
					g.Edge(prevSegmentId, prevSegmentSize, segmentId, segmentSize, k, std::cout);
				}				

				prevSegmentId = segmentId;
				prevSegmentSize = segmentSize;
				begin = end;
			}
			else
			{
				g.FlushPath(currentPath, chrSegmentId[seqId], k, std::cout);
				chrReader.NextChr(chr);
				prevSegmentId = 0;
				begin = end;

				if (begin.GetChr() != ++seqId)
				{
					throw std::runtime_error("The input is corrupted");
				}
			}
		}
	}

	g.FlushPath(currentPath, chrSegmentId[seqId], k, std::cout);
}

/*


void FlushPathGfa2(std::vector<int64_t> & currentPath, int64_t seqId, size_t k)
{
	if (currentPath.size() > 0)
	{
		std::cout << "P\t" << seqId + 1 << '\t';
		std::copy(currentPath.begin(), currentPath.end() - 1, std::ostream_iterator<int64_t>(std::cout, ","));
		std::cout << currentPath.back();
		std::cout << '\t';

		if (currentPath.size() > 1)
		{
			for (size_t i = 0; i < currentPath.size() - 2; i++)
			{
				std::cout << k << "M,";
			}

			std::cout << k << "M";
		}

		std::cout << std::endl;
		currentPath.clear();
	}
}


void GenerateGfa2Output(const std::string & inputFileName, const std::vector<std::string> & genomes, size_t k)
{	
	std::vector<std::string> chrSegmentId;
	std::vector<uint64_t> chrSegmentLength;
	std::cout << "H\tVN:Z:2.0" << std::endl;
	ReadInputSequences(genomes, chrSegmentId, chrSegmentLength, true, true);
	std::vector<int64_t> currentPath;
	const int64_t NO_SEGMENT = 0;
	std::string chr;
	int64_t seqId = NO_SEGMENT;
	int64_t prevSegmentId = NO_SEGMENT;
	TwoPaCo::JunctionPosition end;
	TwoPaCo::JunctionPosition begin;
	TwoPaCo::ChrReader chrReader(genomes);
	TwoPaCo::JunctionPositionReader reader(inputFileName.c_str());
	std::vector<bool> seen(MAX_SEGMENT_NUMBER, 0);
	int64_t previousId = 0;

#ifdef _DEBUG
	std::map<int64_t, std::string> segmentBody;
#endif
	if (reader.NextJunctionPosition(begin))
	{
		chrReader.NextChr(chr);
		while (reader.NextJunctionPosition(end))
		{
			if (begin.GetChr() == end.GetChr())
			{
				Segment nowSegment(begin, end, chr[begin.GetPos() + k], TwoPaCo::DnaChar::ReverseChar(chr[end.GetPos() - 1]));
				int64_t segmentId = nowSegment.GetSegmentId();
				currentPath.push_back(segmentId);
				uint64_t segmentSize = end.GetPos() + k - begin.GetPos();
				if (!seen[Abs(segmentId)])
				{
					std::cout << "S\t" << Abs(segmentId) << "\t" << segmentSize << "\t";
					if (segmentId > 0)
					{
						std::copy(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k, std::ostream_iterator<char>(std::cout));

					}
					else
					{
						std::string buf = TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
						std::copy(buf.begin(), buf.end(), std::ostream_iterator<char>(std::cout));
					}

					std::cout << std::endl;
					seen[Abs(segmentId)] = true;
				}

#ifdef _DEBUG
				int64_t absSegmentId = Abs(segmentId);
				std::string buf = segmentId > 0 ? std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k) :
					TwoPaCo::DnaChar::ReverseCompliment(std::string(chr.begin() + begin.GetPos(), chr.begin() + end.GetPos() + k));
				if (segmentBody.count(absSegmentId) == 0)
				{
					segmentBody[absSegmentId] = buf;
				}
				else
				{
					assert(segmentBody[absSegmentId] == buf);
				}
#endif

				std::cout << "F\t" 
					<< absSegmentId << '\t'
					<< chrSegmentId[seqId] << Sign(segmentId)
					<< "0\t" 
					<< segmentSize << "$" << "\t" 
					<< Gfa2Position(begin.GetPos(), chrSegmentLength[seqId]) << std::endl;

				if (prevSegmentId != NO_SEGMENT)
				{
					std::cout << "E\t" << Abs(prevSegmentId) << '\t' << Sign(prevSegmentId) << '\t' << Abs(segmentId) << '\t' << Sign(segmentId) << '\t' << k << 'M' << std::endl;
				}

				prevSegmentId = segmentId;
				begin = end;
			}
			else
			{
				FlushPathGfa2(currentPath, seqId, k);
				chrReader.NextChr(chr);
				prevSegmentId = 0;
				begin = end;

				if (begin.GetChr() != ++seqId)
				{
					throw std::runtime_error("The input is corrupted");
				}
			}
		}
	}

	FlushPathGfa2(currentPath, seqId, k);
}
*/

int main(int argc, char * argv[])
{
	std::vector<std::string> format;
	format.push_back("seq");
	format.push_back("group");
	format.push_back("dot");
	format.push_back("gfa1");
	format.push_back("gfa2");
	std::stringstream formatString;
	std::copy(format.begin(), format.begin(), std::ostream_iterator<std::string>(formatString, "|"));
	try
	{
		TCLAP::CmdLine cmd("This utility converts the binary output of TwoPaCo to another format", ' ', "0.9.0");
		TCLAP::SwitchArg prefix("", "prefix", "Add a prefix to segments in GFA (in case if you have genomes with identical FASTA headers)", cmd, false);

		TCLAP::UnlabeledValueArg<std::string> inputFileName("infile",
			"input file name",
			true,
			"",
			"file name",
			cmd);

		TCLAP::ValuesConstraint<std::string> formatConstraint(format);
		TCLAP::ValueArg<std::string> outputFileFormat("f",
			"format",
			"Output format",
			true,
			format[0],
			&formatConstraint,
			cmd);

		TCLAP::MultiArg<std::string> seqFileName("s",
			"seqfile",
			"sequences file name",
			false,
			"",
			cmd);

		TCLAP::ValueArg<unsigned int> kvalue("k",
			"kvalue",
			"Value of k",
			true,
			25,
			"integer",
			cmd);

		cmd.parse(argc, argv);
		if (outputFileFormat.getValue() == format[0])
		{
			GenerateOrdinaryOutput(inputFileName.getValue());
		}
		else if (outputFileFormat.getValue() == format[1])
		{
			GenerateGroupOutupt(inputFileName.getValue());
		}
		else if (outputFileFormat.getValue() == format[2])
		{

		}
		else if (outputFileFormat.getValue() == format[3])
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			GenerateGfaOutput(inputFileName.getValue(), seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), Gfa1Generator());
		}
		else if (outputFileFormat.getValue() == format[4])
		{
			if (!seqFileName.isSet())
			{
				throw TCLAP::ArgParseException("Required argument missing\n", "seqfilename");
			}

			//GenerateGfaOutput(inputFileName.getValue(), seqFileName.getValue(), kvalue.getValue(), prefix.getValue(), Gfa2Generator());
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
