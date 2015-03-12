#ifndef _BLOOM_FILTER_
#define _BLOOM_FILTER

#include <cstdlib>
#include <vector>
#include <atomic>

namespace Sibelia
{
	class BloomFilter
	{
	public:
		~BloomFilter();
		BloomFilter(size_t size);
		void Init();
		size_t Size() const;
		void Put(const std::vector<size_t> & hf);
		bool Get(const std::vector<size_t> & hf) const;
	private:
		static const size_t SUCCESS = -1;		
		typedef std::atomic<uint64_t> UInt;
		size_t size_;
		size_t realSize_;
		UInt * filter_;
		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const;
	};
}

#endif