#ifndef _CONCURRENT_BIT_VECTOR_
#define _CONCURRENT_BIT_VECTOR_

#include <vector>
#include <atomic>
#include <memory>
#include <cstdlib>


namespace Sibelia
{
	class ConcurrentBitVector
	{
	public:
		~ConcurrentBitVector();
		ConcurrentBitVector(size_t size);		

		template<class HfPtr, class F>
			void PutInBloomFilter(std::vector<HfPtr> & hf, F f, char farg, uint64_t hash0)
			{
				size_t hashFunctions = hf.size() / BLOCKS;
				for (size_t b = 0; b < BLOCKS; b++)
				{
					uint64_t base = b == 0 ? hash0 : ((*hf[b * hashFunctions]).*f)(farg);
					for (size_t h = 1; h < hashFunctions; h++)
					{
						uint64_t idx = ((*hf[b * hashFunctions + h]).*f)(farg);
						uint64_t bit;
						uint64_t element;
						GetCoord(idx, element, bit);
						filter_[base + element].fetch_or(uint32_t(1) << uint32_t(bit));
					}
				}				
			}

		template<class HfPtr, class F>
			bool IsInBloomFilter(const std::vector<HfPtr> & hf, F f, char farg, uint64_t hash0) const
			{
				size_t hashFunctions = hf.size() / BLOCKS;
				for (size_t b = 0; b < BLOCKS; b++)
				{
					uint64_t base = b == 0 ? hash0 : ((*hf[b * hashFunctions]).*f)(farg);
					for (size_t h = 1; h < hashFunctions; h++)
					{
						uint64_t idx = ((*hf[b * hashFunctions + h]).*f)(farg);
						uint64_t bit;
						uint64_t element;
						GetCoord(idx, element, bit);
						if ((filter_[base + element] & (uint32_t(1) << uint32_t(bit))) == 0)
						{
							return false;
						}
					}
				}
				
				return true;
		}
		
		static const size_t BLOCKS = 3;
		static const size_t MINOR_BITS = 9;
		
	private:

		void Init();
		void SetConcurrently(size_t idx);		

		typedef uint32_t UInt;
		size_t size_;
		size_t realSize_;
		std::atomic<UInt> * filter_;
		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const;
	};

}

#endif