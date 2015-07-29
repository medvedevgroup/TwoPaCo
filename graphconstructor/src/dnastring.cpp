#include <sstream>
#include <algorithm>

#include "DnaString.h"

namespace Sibelia
{/*
	
	
	bool operator == (const DnaString & a, const DnaString & b)
	{
		if (a.GetSize() == b.GetSize())
		{
			size_t remain = a.size_;			
			for (size_t i = 0; remain > 0; i++)
			{
				size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
				uint64_t apiece = a.str_[i];
				uint64_t bpiece = b.str_[i];				
				if (apiece != bpiece)
				{
					return false;
				}

				remain -= current;
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
		if (a.GetSize() == b.GetSize())
		{
			size_t remain = a.size_;
			for (size_t i = 0; remain > 0; i++)
			{
				size_t current = std::min(remain, DnaString::UNIT_CAPACITY);
				uint64_t apiece = a.str_[i];
				uint64_t bpiece = b.str_[i];
				if (apiece != bpiece)
				{
					return apiece < bpiece;
				}

				remain -= current;
			}
		}

		return false;
	}*/
}
