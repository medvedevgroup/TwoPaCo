#include <cassert>
#include <iostream>

#include "concurrentbitvector.h"

namespace Sibelia
{
	ConcurrentBitVector::ConcurrentBitVector(size_t size)
		: size_(size), realSize_(size / 64 + 1), filter_(new UInt[realSize_])
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
		filter_[element].fetch_or(uint64_t(1) << uint64_t(bit));
	}

	bool ConcurrentBitVector::Get(size_t idx) const
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		return (filter_[element] & (uint64_t(1) << uint64_t(bit))) != 0;
	}

	void ConcurrentBitVector::GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
	{
		bit = idx & ((uint64_t(1) << uint64_t(6)) - 1);
		element = idx >> 6;
		assert(element < size_ / 64 + 1);
	}

	ConcurrentBitVector::~ConcurrentBitVector()
	{
		delete[] filter_;
	}
}