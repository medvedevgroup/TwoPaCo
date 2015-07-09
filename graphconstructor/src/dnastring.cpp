#include <sstream>
#include <algorithm>

#include "DnaString.h"

namespace Sibelia
{
	const std::string DnaString::LITERAL = "ACGT";

	/*
	bool operator == (const DnaString & a, const DnaString & b)
	{
		if (a.GetSize() == b.GetSize())
		{
			for (size_t i = 0; i < a.capacity_; i++)
			{
				if (a.str_[i] != b.str_[i])
				{
					return false;
				}
			}

			return true;
		}

		return false;
	}

	bool operator != (const DnaString & a, const DnaString & b)
	{
		return !(a == b);
	}
	
	bool operator < (const DnaString & a, const DnaString & b)
	{
		for (size_t i = 0; i < a.capacity_; i++)
		{
			if (a.str_[i] != b.str_[i])
			{
				return a.str_[i] != b.str_[i];
			}
		}

		return false;
	}	*/
}
