#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

namespace Sibelia
{
	class DnaString
	{
	public:		
		void PopBack();
		void PopFront();
		void AppendBack(char ch);
		void AppendFront(char ch);
		DnaString(uint64_t size = 0);
		char GetChar(uint64_t idx) const;
		void SetChar(uint64_t idx, char ch);
		size_t GetSize() const;
		uint64_t GetBody() const;
		std::string ToString() const;
		static const std::string LITERAL;
	private:
		uint64_t size_;
		uint64_t body_;
		static uint64_t MakeUp(char ch);
	};
}

#endif