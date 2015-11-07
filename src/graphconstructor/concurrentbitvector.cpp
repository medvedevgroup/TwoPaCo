#include <cassert>
#include <iostream>

#include "concurrentbitvector.h"

namespace Sibelia
{
	ConcurrentBitVector::ConcurrentBitVector(size_t size)
		: size_(size), realSize_(size / 32 + 1), filter_(new UInt[realSize_])
	{
		Init();
	}

	void ConcurrentBitVector::Init()
	{
		for (size_t i = 0; i < realSize_; i++)
		{
			filter_[i] = 0;
		}
	}

	size_t ConcurrentBitVector::Size() const
	{
		return size_;
	}

	void ConcurrentBitVector::SetConcurrently(size_t idx)
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		filter_[element].fetch_or(uint32_t(1) << uint32_t(bit));
	}

	bool ConcurrentBitVector::Get(size_t idx) const
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		return (filter_[element] & (uint32_t(1) << uint32_t(bit))) != 0;
	}

	void ConcurrentBitVector::GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
	{
		bit = idx & ((uint32_t(1) << uint64_t(5)) - 1);
		element = idx >> 5;
		assert(element < size_ / 32 + 1);
	}

	ConcurrentBitVector::~ConcurrentBitVector()
	{
		delete[] filter_;
	}
}