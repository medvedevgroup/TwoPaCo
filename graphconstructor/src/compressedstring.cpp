#include "compressedstring.h"

namespace Sibelia
{
	const std::string LITERAL = "ACGT";
	extern const size_t UNIT_CAPACITY = 32;

	std::string RevComp(const std::string & str)
	{
		std::string ret;
		for (auto it = str.rbegin(); it != str.rend(); ++it)
		{
			ret.push_back(CompressedString<1>::ReverseChar(*it));
		}

		return ret;
	}
}