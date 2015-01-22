#include <algorithm>

#include "dnastring.h"
#include "lib/SpookyV2.h"

namespace Sibelia
{
	const std::string DnaString::LITERAL = "ACGT";

	DnaString::DnaString(uint64_t size) : size_(size), body_(0)
	{
	}

	uint64_t DnaString::MakeUp(char ch)
	{
		ch = toupper(ch);
		size_t idx = std::find(LITERAL.begin(), LITERAL.end(), ch) - LITERAL.begin();
		if (idx != LITERAL.size())
		{
			return idx;
		}

		return rand() % LITERAL.size();
	}

	void DnaString::PopBack()
	{
		uint64_t mask = uint64_t(0x3) << (--size_ * 2);
		body_ &= ~(mask);
	}

	void DnaString::PopFront()
	{		
		--size_;
		body_ &= ~0x3;
		body_ >>= 2;
	}

	void DnaString::AppendFront(char ch)
	{
		++size_ ;
		body_ <<= 2;
		body_ |= MakeUp(ch);
	}

	void DnaString::AppendBack(char ch)
	{
		body_ |= MakeUp(ch) << (2 * size_++);
	}

	size_t DnaString::GetSize() const 
	{
		return size_;
	}

	uint64_t DnaString::GetBody() const
	{
		return body_;
	}

	char DnaString::GetChar(uint64 idx) const
	{
		uint64_t charIdx = body_ >> (2 * idx);
		return LITERAL[charIdx & 0x3];
	}

	void DnaString::SetChar(uint64_t idx, char ch)
	{
		uint64_t mask = uint64_t(0x3) << (idx * 2);
		body_ &= ~(mask);
		body_ |= MakeUp(ch) << (2 * idx++);
	}

	std::string DnaString::ToString() const
	{
		std::string ret(size_, ' ');
		for (size_t i = 0; i < size_; i++)
		{
			ret[i] = GetChar(i);
		}

		return ret;
	}
}
