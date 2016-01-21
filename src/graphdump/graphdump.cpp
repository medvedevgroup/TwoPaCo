#include <vector>
#include <iostream>
#include <algorithm>
#include <spooky/SpookyV2.h>
#include <junctionapi/junctionapi.h>

bool CompareJunctionsById(const Sibelia::JunctionPosition & a, const Sibelia::JunctionPosition & b)
{
	return a.GetId() < b.GetId();
}

bool CompareJunctionsByPos(const Sibelia::JunctionPosition & a, const Sibelia::JunctionPosition & b)
{
	return std::make_pair(a.GetChr(), a.GetPos()) < std::make_pair(b.GetChr(), b.GetPos());
}

struct EqClass
{
	size_t label;
	std::vector<Sibelia::JunctionPosition> position;
};

bool CompareJunctionClasses(const EqClass & a, const EqClass & b)
{
	return CompareJunctionsByPos(a.position[0], b.position[0]);
}

int main(int argc, char * argv[])
{
	Sibelia::JunctionPosition pos;
	Sibelia::JunctionPositionReader reader(argv[1]);
	std::vector<EqClass> eqClass;
	std::vector<Sibelia::JunctionPosition> junction;	
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

	return 0;
}