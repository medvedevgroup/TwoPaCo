#include <cstdio>
#include <fstream>
#include <cassert>
#include <iostream>
#include <stdexcept>

#include "concurrentbitvector.h"

namespace TwoPaCo
{
	ConcurrentBitVector::ConcurrentBitVector(size_t size)
		: size_(size), realSize_(size / BASIC_TYPE_BITS + 1), filter_(new UInt[realSize_])
	{
		Reset();		
	}

	ConcurrentBitVector::BASIC_TYPE ConcurrentBitVector::GetElement(size_t idx) const
	{
		return filter_[idx];
	}

	void ConcurrentBitVector::OrElementCouncerrently(size_t idx, BASIC_TYPE mask)
	{
		filter_[idx].fetch_or(mask);
	}

	void ConcurrentBitVector::Reset()
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

	void ConcurrentBitVector::SetBitConcurrently(size_t idx)
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		filter_[element].fetch_or(BASIC_TYPE(1) << BASIC_TYPE(bit));
	}

	bool ConcurrentBitVector::GetBit(size_t idx) const
	{
		uint64_t bit;
		uint64_t element;
		GetCoord(idx, element, bit);
		return (filter_[element] & (BASIC_TYPE(1) << BASIC_TYPE(bit))) != 0;
	}

	void ConcurrentBitVector::GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
	{
		bit = idx & ((uint64_t(1) << uint64_t(BASIC_TYPE_POWER)) - uint64_t(1));
		element = idx >> BASIC_TYPE_POWER;
		assert(element < size_ / BASIC_TYPE_BITS + 1);
	}

	ConcurrentBitVector::~ConcurrentBitVector()
	{
		delete[] filter_;
	}

	void ConcurrentBitVector::WriteToFile(const std::string & fileName) const
	{
		std::ofstream candidateMaskFile(fileName.c_str(), std::ios::binary);
		if (!candidateMaskFile)
		{
			throw std::runtime_error("Can't open a temporary file");
		}

		if (!candidateMaskFile.write(reinterpret_cast<const char*>(filter_), realSize_ * sizeof(BASIC_TYPE)))
		{
			throw std::runtime_error("Can't write to a temporary file");
		}
	}

	void ConcurrentBitVector::ReadFromFile(const std::string & fileName, bool cleanUp)
	{
		{
			std::ifstream candidateMaskFile(fileName.c_str(), std::ios::binary);
			if (!candidateMaskFile)
			{
				throw std::runtime_error("Can't open a temporary file");
			}

			if (!candidateMaskFile.read(reinterpret_cast<char*>(filter_), realSize_ * sizeof(BASIC_TYPE)))
			{
				throw std::runtime_error("Can't read from a temporary file");
			}
		}

		if (cleanUp)
		{
			std::remove(fileName.c_str());
		}
	}
}
