#include <vector>
#include <iostream>
#include <algorithm>

#include <tclap/CmdLine.h>
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
	size_t label;
	std::vector<TwoPaCo::JunctionPosition> position;
};

bool CompareJunctionClasses(const EqClass & a, const EqClass & b)
{
	return CompareJunctionsByPos(a.position[0], b.position[0]);
}

int main(int argc, char * argv[])
{
	try
	{
		TCLAP::CmdLine cmd("This utility converts the binary format into human readable one", ' ', "0");

		TCLAP::UnlabeledValueArg<std::string> inputFileName("infile",
			"input file name",
			true,
			"",
			"file name",
			cmd);
		
		cmd.parse(argc, argv);
		TwoPaCo::JunctionPosition pos;
		TwoPaCo::JunctionPositionReader reader(inputFileName.getValue().c_str());
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

		std::sort(eqClass.begin(), eqClass.end(), CompareJunctionClasses);
		for (auto junctionClass : eqClass)
		{
			for (auto j : junctionClass.position)
			{
				std::cout << j.GetChr() << ' ' << j.GetPos() << "; ";
			}

			std::cout << std::endl;
		}
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}	

	return 0;
}
