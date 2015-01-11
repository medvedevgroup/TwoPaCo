#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

namespace Sibelia
{
	class DnaString
	{
	public:
		DnaString();
		void PopBack();
		void PopFront();
		void AppendBack(char ch);
		void AppendFront(char ch);
		char GetChar(size_t idx) const;
		size_t GetSize() const;
		uint64_t GetBody() const;
		std::string ToString() const;
		static const std::string LITERAL;
	private:
		size_t size_;
		uint64_t body_;
		static size_t MakeUp(char ch);
	};
}

#endif