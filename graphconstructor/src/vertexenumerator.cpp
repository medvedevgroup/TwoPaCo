#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <numeric>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vertexenumerator.h"

namespace Sibelia
{
	std::string RevComp(const std::string & str)
	{
		std::string ret;
		for (auto it = str.rbegin(); it != str.rend(); ++it)
		{
			ret.push_back(VertexEnumerator::DnaString::ReverseChar(*it));
		}

		return ret;
	}
}