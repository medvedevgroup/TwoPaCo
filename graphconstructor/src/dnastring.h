#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>
#include <cstdint>

namespace Sibelia
{
	class DnaString
	{
	public:
		void Clear();
		char PopBack();
		char PopFront();	
		void AppendBack(char ch);
		void AppendFront(char ch);
		size_t GetSize() const;
		size_t MaxSize() const;
		size_t BytesInBody() const;
		uint64_t* GetBody() const;		
		void Assign(const uint64_t * str);
		char GetChar(uint64_t idx) const;
		void SetChar(uint64_t idx, char ch);
		~DnaString();		
		DnaString(DnaString && str);
		DnaString(const DnaString & str);	
		DnaString(const std::string & str);
		DnaString(size_t maxSize, size_t size);
		DnaString(size_t size, const uint64_t * body);
				
		static char Reverse(char ch);
		std::string ToString() const;
		static uint64_t MakeUp(char ch);		
		DnaString RevComp() const;
		static std::string RevComp(const std::string & str);
		static const std::string LITERAL;
	private:		
		size_t size_;
		size_t capacity_;
		uint64_t * str_;		
		static const size_t UNIT_CAPACITY = 32;		
		uint64_t TranslateIdx(uint64_t & idx) const;
		static size_t CalculateCapacity(size_t size);
		char GetChar(uint64_t element, uint64_t idx) const;
		void SetChar(uint64_t element, uint64_t idx, char ch);
		friend bool operator < (const DnaString & a, const DnaString & b);
		friend bool operator == (const DnaString & a, const DnaString & b);
		friend bool operator != (const DnaString & a, const DnaString & b);
	};

	bool operator < (const DnaString & a, const DnaString & b);
	bool operator == (const DnaString & a, const DnaString & b);
	bool operator != (const DnaString & a, const DnaString & b);
}

#endif