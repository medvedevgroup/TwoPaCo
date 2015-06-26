#include <algorithm>

#include "dnastring.h"
#include "lib/SpookyV2.h"

namespace Sibelia
{
	DnaString::DnaString(const std::string & body) : size_(0), body_(0)
	{
		for (size_t i = 0; i < body.size(); i++)
		{
			AppendBack(body[i]);
		}
	}

	bool operator == (const DnaString & a, const DnaString & b)
	{
		return a.size_ == b.size_ && a.body_ == b.body_;
	}

	bool operator != (const DnaString & a, const DnaString & b)
	{
		return !(a == b);
	}

	DnaString::DnaString(uint64_t size) : size_(size), body_(0)
	{
	}

	DnaString::DnaString(uint64_t size, uint64_t body) : size_(size), body_(body)
	{
		if (size != sizeof(body) * 4)
		{
			uint64_t mask = (uint64_t(1) << (size * 2));
			body_ &= (mask - 1);
		}
	}

	char DnaString::PopBack()
	{
		char ret = GetChar(size_ - 1);
		uint64_t mask = uint64_t(0x3) << (--size_ * 2);
		body_ &= ~(mask);
		return ret;
	}

	char DnaString::PopFront()
	{	
		char ret = GetChar(0);
		--size_;
		body_ &= ~0x3;
		body_ >>= 2;
		return ret;
	}

	void DnaString::AppendFront(uint64_t ch)
	{
		++size_ ;
		body_ <<= 2;
		body_ |= ch;
	}

	void DnaString::AppendBack(uint64_t ch)
	{
		body_ |=  ch << (2 * size_++);
	}

	size_t DnaString::GetSize() const 
	{
		return size_;
	}

	uint64_t DnaString::GetBody() const
	{
		return body_;
	}

	uint64_t DnaString::Prefix(uint64_t body, uint64_t size)
	{
		uint64_t mask = (uint64_t(1) << (size * 2));
		return body & (mask - 1);
	}

	DnaString DnaString::Prefix(uint64_t size) const
	{
		if (size == 0)
		{
			return DnaString();
		}

		if (size >= size_)
		{
			return *this;
		}

		uint64_t mask = (uint64_t(1) << (size * 2));
		return DnaString(size, body_ & (mask - 1));
	}

	char DnaString::GetChar(uint64 idx) const
	{
		uint64_t charIdx = body_ >> (2 * idx);
		return charIdx & 0x3;
	}

	void DnaString::SetChar(uint64_t idx, uint64_t ch)
	{
		uint64_t mask = uint64_t(0x3) << (idx * 2);
		body_ &= ~(mask);
		body_ |= ch << (2 * idx++);
	}

	std::string DnaString::ToString(bool unMakeUp) const
	{
		std::string ret(size_, ' ');
		for (size_t i = 0; i < size_; i++)
		{
			ret[i] = unMakeUp ? StreamFastaParser::UnMakeUp(GetChar(i)) : GetChar(i);
		}

		return ret;
	}

	DnaString DnaString::RevComp() const
	{
		DnaString ret;
		for (size_t i = 0; i < size_; i++)
		{
			ret.AppendFront(StreamFastaParser::Reverse(GetChar(i)));
		}

		return ret;
	}

	std::string DnaString::RevComp(const std::string & str)
	{
		std::string ret;
		std::transform(str.rbegin(), str.rend(), std::back_inserter(ret), StreamFastaParser::Reverse);
		return ret;
	}
}
