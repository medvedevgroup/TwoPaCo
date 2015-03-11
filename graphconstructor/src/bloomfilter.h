#ifndef _BLOOM_FILTER_
#define _BLOOM_FILTER

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
		typedef std::atomic<unsigned char> UChar;
		size_t size_;
		UChar * filter_;
		void GetCoord(size_t idx, size_t & element, size_t & bit) const;
	};
}

#endif