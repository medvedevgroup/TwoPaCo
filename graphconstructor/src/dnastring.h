#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

namespace Sibelia
{
	class DnaPiece
	{
	public:		
		char PopFront();
		void AppendFront(char ch);
		DnaPiece();		
		DnaPiece(uint64_t size, uint64_t body);
		DnaPiece(std::string::const_iterator begin, std::string::const_iterator end);
		char GetChar(uint64_t idx) const;
		void SetChar(uint64_t idx, char ch);
		static uint64_t Prefix(uint64_t body, uint64_t size);
		size_t GetSize() const;
		uint64_t GetBody() const;
		std::string ToString(size_t size) const;
		static const std::string LITERAL;
		static const size_t CAPACITY = 32;
	private:
		uint64_t body_;
		static uint64_t MakeUp(char ch);
		friend bool operator == (const DnaPiece & a, const DnaPiece & b);
		friend bool operator != (const DnaPiece & a, const DnaPiece & b);
	};

	bool operator == (const DnaPiece & a, const DnaPiece & b);
	bool operator != (const DnaPiece & a, const DnaPiece & b);

	class DnaString
	{
	public:
		char PopBack();
		char PopFront();		
		void AppendBack(char ch);
		void AppendFront(char ch);
		size_t GetSize() const;
		size_t MaxSize() const;
		uint64_t* GetBody() const;
		size_t GetCapacity() const;
		char GetChar(uint64_t idx) const;
		void SetChar(uint64_t idx, char ch);
		~DnaString();		
		DnaString(DnaString && str);
		DnaString(const DnaString & str);		
		DnaString(const std::string & str);		
		DnaString(size_t size, const uint64_t * body);
		DnaString(size_t size, bool capacity = false);
		
		DnaString RevComp() const;
		static char Reverse(char ch);
		std::string ToString() const;
		static std::string RevComp(const std::string & str);
		

		static const std::string LITERAL;
	private:		
		size_t size_;
		size_t capacity_;
		DnaPiece * str_;
		static size_t CalculateCapacity(size_t size);
		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & pos) const;
	};
}

#endif