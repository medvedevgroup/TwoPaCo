#include "compressedstring.h"

namespace Sibelia
{
	const std::string LITERAL = "ACGT";
	extern const size_t UNIT_CAPACITY = 32;

	namespace
	{
		class ReverseTable
		{
		public:
			static char table[1 << 8];
			ReverseTable()
			{
				std::fill(table, table + (1 << 8), 'N');
				table['A'] = 'T';
				table['T'] = 'A';
				table['C'] = 'G';
				table['G'] = 'C';
			}

		};

		char ReverseTable::table[256];
		ReverseTable table;
	}

	char ReverseChar(char ch)
	{
		return ReverseTable::table[ch];
	}

	std::string RevComp(const std::string & str)
	{
		std::string ret;
		for (auto it = str.rbegin(); it != str.rend(); ++it)
		{
			ret.push_back(ReverseChar(*it));
		}

		return ret;
	}
}