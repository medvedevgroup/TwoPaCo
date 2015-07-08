#include <sstream>
#include <algorithm>

#include "DnaString.h"

namespace Sibelia
{
	const std::string DnaString::LITERAL = "ACGT";
//
//	DnaPiece::DnaPiece() : body_(0)
//	{
//	}
//
//	DnaPiece::DnaPiece(uint64_t size, uint64_t body) : body_(body)
//	{
//		if (size != sizeof(body) * 4)
//		{
//			uint64_t mask = (uint64_t(1) << (size * 2));
//			body_ &= (mask - 1);
//		}
//	}
//
	uint64_t DnaString::MakeUp(char ch)
	{
		switch (ch)
		{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		}

		return 0;
	}
//
//	char DnaPiece::PopFront()
//	{
//		char ret = GetChar(0);
//		body_ &= ~0x3;
//		body_ >>= 2;
//		return ret;
//	}
//
//	void DnaPiece::AppendFront(char ch)
//	{
//		body_ <<= 2;
//		body_ |= DnaString::MakeUp(ch);
//	}
//
//	uint64_t DnaPiece::GetBody() const
//	{
//		return body_;
//	}
//
//

//
	
//
//

//
//	bool operator == (const DnaPiece & a, const DnaPiece & b)
//	{
//		return a.GetBody() == b.GetBody();
//	}
//
//	bool operator != (const DnaPiece & a, const DnaPiece & b)
//	{
//		return !(a == b);
//	}
//
//	//-------------------------------------------------------------------------

	void DnaString::Clear()
	{
		size_ = 0;
	}

	DnaString::~DnaString()
	{
		delete[] str_;
	}

	DnaString::DnaString(DnaString && str) : size_(str.size_), capacity_(capacity_), str_(0)
	{
		std::swap(str_, str.str_);
	}

	DnaString::DnaString(size_t maxSize, size_t size) : size_(size), capacity_(CalculateCapacity(maxSize)),
		str_(new uint64_t[capacity_])
	{

	}

	uint64_t DnaString::TranslateIdx(uint64_t & idx) const
	{
		uint64_t ret = idx >> 5;
		idx = idx & ((uint64_t(1) << uint64_t(5)) - 1);
		return ret;
	}

	char DnaString::GetChar(uint64_t idx) const
	{		
		uint64_t element = TranslateIdx(idx);
		return GetChar(element, idx);
	}

	void DnaString::SetChar(uint64_t idx, char ch)
	{
		uint64_t element = TranslateIdx(idx);
		SetChar(element, idx, ch);
	}

	char DnaString::GetChar(uint64_t element, uint64_t idx) const
	{
		uint64_t charIdx = str_[element] >> (2 * idx);
		return DnaString::LITERAL[charIdx & 0x3];
	}

	void DnaString::SetChar(uint64_t element, uint64_t idx, char ch)
	{
		uint64_t mask = uint64_t(0x3) << (idx * 2);
		str_[element] &= ~(mask);
		str_[element] |= DnaString::MakeUp(ch) << (2 * idx++);
	}

	size_t DnaString::GetSize() const
	{
		return size_;
	}

	size_t DnaString::MaxSize() const
	{
		return capacity_ * UNIT_CAPACITY;
	}
	
	uint64_t* DnaString::GetBody() const
	{
		return reinterpret_cast<uint64_t*>(str_);
	}

	DnaString::DnaString(const std::string & str) :size_(str.size()), capacity_(CalculateCapacity(str.size()))
	{
		for (size_t i = 0; i < str.size(); i++)
		{
			SetChar(i, str[i]);
		}
	}
	
	DnaString::DnaString(size_t size, const uint64_t * body) :size_(size), capacity_(CalculateCapacity(size))
	{
		std::copy(body, body + capacity_, reinterpret_cast<uint64_t*>(str_));
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

	std::string DnaString::RevComp(const std::string & str)
	{
		std::string ret;
		std::transform(str.rbegin(), str.rend(), std::back_inserter(ret), DnaString::Reverse);
		return ret;
	}

	size_t DnaString::CalculateCapacity(size_t size)
	{
		size_t mainPart = size / UNIT_CAPACITY;
		return mainPart + (mainPart * UNIT_CAPACITY == size ? 0 : 1);
	}

	char DnaString::Reverse(char ch)
	{
		switch (ch)
		{
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		}

		return 'N';
	}

	size_t DnaString::BytesInBody() const
	{
		return capacity_ * sizeof(str_[0]);
	}

	char DnaString::PopBack()
	{
		return GetChar(size_--);
	}

	char DnaString::PopFront()
	{		
		size_t remain = size_;
		char ret = GetChar(0);
		for (size_t pos = 0; remain > 0; ++pos)
		{
			size_t current = std::min(remain, UNIT_CAPACITY);
			char prev = GetChar(pos, 0);
			str_[pos] &= ~0x3;
			str_[pos] >>= 2;
			remain -= current;
			if (pos > 0)
			{				
				SetChar(pos - 1, UNIT_CAPACITY - 1, prev);
			}			
		}

		--size_;
		return ret;
	}

	void DnaString::AppendBack(char ch)
	{
		SetChar(size_++, ch);
	}

	void DnaString::AppendFront(char ch)
	{
		size_t remain = ++size_;
		for (size_t pos = 0; remain > 0; ++pos)
		{			
			size_t current = std::min(remain, UNIT_CAPACITY);
			char nextCh = GetChar(pos,  current - 1);
			str_[pos] <<= 2;
			str_[pos] |= DnaString::MakeUp(ch);
			ch = nextCh;
			remain -= current;
		}
	}

	DnaString DnaString::RevComp() const
	{
		DnaString ret(size_, size_);
		for (size_t i = 0; i < size_; i++)
		{
			ret.SetChar(size_ - i - 1, Reverse(GetChar(i)));
		}

		return ret;
	}

	bool operator == (const DnaString & a, const DnaString & b)
	{
		if (a.GetSize() == b.GetSize())
		{
			for (size_t i = 0; i < a.capacity_; i++)
			{
				if (a.str_[i] != b.str_[i])
				{
					return false;
				}
			}

			return true;
		}

		return false;
	}

	bool operator != (const DnaString & a, const DnaString & b)
	{
		return !(a == b);
	}
	
	bool operator < (const DnaString & a, const DnaString & b)
	{
		for (size_t i = 0; i < a.capacity_; i++)
		{
			if (a.str_[i] != b.str_[i])
			{
				return a.str_[i] != b.str_[i];
			}
		}

		return false;
	}	
}
