#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

namespace Sibelia
{
	class DnaString
	{
	public:		
		char PopBack();
		char PopFront();
		void AppendBack(char ch);
		void AppendFront(char ch);
		DnaString(uint64_t size = 0);
		DnaString(uint64_t size, uint64_t body);
		char GetChar(uint64_t idx) const;
		void SetChar(uint64_t idx, char ch);
		DnaString Prefix(uint64_t size) const;
		size_t GetSize() const;
		uint64_t GetBody() const;
		std::string ToString() const;
		static const std::string LITERAL;
		static char Reverse(char ch);
		DnaString RevComp() const;
		static std::string RevComp(const std::string & str);
	private:
		uint64_t size_;
		uint64_t body_;
		static uint64_t MakeUp(char ch);
		friend bool operator == (const DnaString & a, const DnaString & b);
		friend bool operator != (const DnaString & a, const DnaString & b);
	};

	bool operator == (const DnaString & a, const DnaString & b);
	bool operator != (const DnaString & a, const DnaString & b);
}

#endif