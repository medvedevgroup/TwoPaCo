#include <cassert>
#include <iostream>

#include "concurrentbitvector.h"

namespace Sibelia
{
	ConcurrentBitVector::ConcurrentBitVector(size_t size)
		: size_(size), realSize_(size_t(1) << (size - 5)), filter_(new std::atomic<UInt>[realSize_])
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

	void ConcurrentBitVector::GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
	{
		bit = idx & ((uint32_t(1) << uint64_t(5)) - 1);
		element = idx >> 5;
	}

	ConcurrentBitVector::~ConcurrentBitVector()
	{
		delete[] filter_;
	}
}