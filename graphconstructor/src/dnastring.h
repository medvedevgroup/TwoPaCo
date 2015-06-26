#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

#include "streamfastaparser.h"

namespace Sibelia
{
	class DnaString
	{
	public:		
		char PopBack();
		char PopFront();
		void AppendBack(uint64_t ch);
		void AppendFront(uint64_t ch);
		DnaString(uint64_t size = 0);
		DnaString(const std::string & body);
		DnaString(uint64_t size, uint64_t body);		
		char GetChar(uint64_t idx) const;
		void SetChar(uint64_t idx, uint64_t ch);
		DnaString Prefix(uint64_t size) const;
		static uint64_t Prefix(uint64_t body, uint64_t size);

		size_t GetSize() const;
		uint64_t GetBody() const;
		std::string ToString(bool unMakeUp = false) const;
		DnaString RevComp() const;
		static std::string RevComp(const std::string & str);
	private:
		uint64_t size_;
		uint64_t body_;
		friend bool operator == (const DnaString & a, const DnaString & b);
		friend bool operator != (const DnaString & a, const DnaString & b);
	};

	bool operator == (const DnaString & a, const DnaString & b);
	bool operator != (const DnaString & a, const DnaString & b);
}

#endif