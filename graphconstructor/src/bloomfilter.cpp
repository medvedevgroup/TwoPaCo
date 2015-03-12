#include "bloomfilter.h"
#include <cassert>
#include <iostream>
namespace Sibelia
{
	BloomFilter::BloomFilter(size_t size) : size_(size), realSize_(size / 64 + 1), filter_(new UInt[realSize_])
	{
		Init();
	}

	void BloomFilter::Init()
	{
		for (size_t i = 0; i < realSize_; i++)
		{
			filter_[i] = 0;
		}
	}

	size_t BloomFilter::Size() const
	{
		return size_;
	}

	void BloomFilter::Put(const std::vector<size_t> & hf)
	{
		for (size_t idx : hf)
		{ 
			uint64_t bit;
			uint64_t element;
			uint64_t oldValue;
			uint64_t newValue;
			GetCoord(idx, element, bit);
			do
			{
				oldValue = filter_[element].load();
				newValue = oldValue | (uint64_t(1) << uint64_t(bit));
			} while (!filter_[element].compare_exchange_strong(oldValue, newValue));
		}
	}

	bool BloomFilter::Get(const std::vector<size_t> & hf) const
	{
		for (size_t idx : hf)
		{
			uint64_t bit;
			uint64_t element;
			GetCoord(idx, element, bit);
			if ((filter_[element] & (uint64_t(1) << uint64_t(bit))) == 0)
			{
				return false;
			}
		}

		return true;
	}

	void BloomFilter::GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
	{
		bit = idx & ((uint64_t(1) << uint64_t(6)) - 1);
		element = idx >> 6;
		assert(element < size_ / 64 + 1);
	}

	BloomFilter::~BloomFilter()
	{
		delete[] filter_;
	}
}