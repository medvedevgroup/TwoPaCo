#include "dnachar.h"

namespace Sibelia
{
	bool DnaChar::isValid_[CHAR_SIZE];
	bool DnaChar::isDefinite_[CHAR_SIZE];
	char DnaChar::reverseTable_[CHAR_SIZE];
	const std::string DnaChar::LITERAL = "ACGT";
	const std::string DnaChar::VALID_CHARS = "ACGTURYKMSWBDHWNX";
	
	namespace
	{
		DnaChar helper;
	}
	

	DnaChar::DnaChar()
	{
		std::fill(reverseTable_, reverseTable_ + CHAR_SIZE, 'N');
		reverseTable_['A'] = 'T';
		reverseTable_['T'] = 'A';
		reverseTable_['C'] = 'G';
		reverseTable_['G'] = 'C';
		std::fill(isValid_, isValid_ + CHAR_SIZE, false);
		std::fill(isDefinite_, isDefinite_ + CHAR_SIZE, false);
		for (char ch : LITERAL)
		{
			isDefinite_[ch] = true;
		}

		for (char ch : VALID_CHARS)
		{
			isValid_[ch] = true;
		}
	}

	bool DnaChar::IsValid(char ch)
	{
		return isValid_[ch];
	}

	bool DnaChar::IsDefinite(char ch)
	{
		return isDefinite_[ch];
	}

	char DnaChar::ReverseChar(char ch)
	{
		return reverseTable_[ch];
	}

	std::string DnaChar::ReverseCompliment(const std::string & str)
	{
		std::string ret;
		for (auto it = str.rbegin(); it != str.rend(); ++it)
		{
			ret.push_back(DnaChar::ReverseChar(*it));
		}

		return ret;
	}

	bool DnaChar::LessSelfReverseComplement(std::string::const_iterator pit, size_t vertexSize)
	{
		std::string::const_reverse_iterator nit(pit + vertexSize);
		for (size_t i = 0; i < vertexSize; i++)
		{
			char reverse = DnaChar::ReverseChar(*nit);
			if (*pit != reverse)
			{
				return *pit < reverse;
			}

			++nit;
			++pit;
		}

		return false;
	}
}