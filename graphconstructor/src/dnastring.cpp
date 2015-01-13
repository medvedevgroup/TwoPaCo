#include "dnastring.h"
#include "lib/SpookyV2.h"

namespace Sibelia
{
	const std::string DnaString::LITERAL = "ACGT";

	DnaString::DnaString() : size_(0), body_(0)
	{
	}

	size_t DnaString::MakeUp(char ch)
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
		uint64_t mask = 0x3 << (--size_ * 2);
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

	char DnaString::GetChar(size_t idx) const
	{
		size_t charIdx = body_ >> (2 * idx);
		return LITERAL[charIdx & 0x3];
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
