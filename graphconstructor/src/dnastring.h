#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

namespace Sibelia
{
	class DnaString
	{
	public:
		std::string ToString() const;

	private:
		uint64_t body;

	};
}

#endif