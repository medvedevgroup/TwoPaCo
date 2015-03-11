#include "bloomfilter.h"
#include <iostream>
namespace Sibelia
{
	BloomFilter::BloomFilter(size_t size) : size_(size), realSize_(size / 8 + 1), filter_(new UChar[realSize_])
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
			size_t bit;
			size_t element;
			GetCoord(idx, element, bit);
			while (true)
			{
				unsigned char oldValue = filter_[element];
				unsigned char newValue = oldValue | 1 << bit;
				if (filter_[element].compare_exchange_strong(oldValue, newValue))
				{
					break;
				}
			}
		}
	}

	bool BloomFilter::Get(const std::vector<size_t> & hf) const
	{
		for (size_t idx : hf)
		{
			size_t bit;
			size_t element;
			GetCoord(idx, element, bit);
			if ((filter_[element] & (1 << bit)) == 0)
			{
				return false;
			}
		}

		return true;
	}

	void BloomFilter::GetCoord(size_t idx, size_t & element, size_t & bit) const
	{
		bit = idx & 0x7;
		element = idx >> 3;
		assert(element < size_ / 8 + 1);
	}

	BloomFilter::~BloomFilter()
	{
		delete[] filter_;
	}
}