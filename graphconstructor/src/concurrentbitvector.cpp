#include <cassert>
#include <iostream>

#include "concurrentbitvector.h"

namespace Sibelia
{
	ConcurrentBitVector::ConcurrentBitVector(size_t size) : ConcurrentBase(size)
	{

	}

	void ConcurrentBitVector::Set(size_t idx)
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		filter_[element].fetch_or(uint32_t(1) << uint32_t(bit));
	}

	NonConcurrentBitVector::NonConcurrentBitVector(size_t size) : NonConcurrentBase(size)
	{

	}

	bool NonConcurrentBitVector::Get(size_t idx) const
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		return (filter_[element] & (uint32_t(1) << uint32_t(bit))) != 0;
	}
}