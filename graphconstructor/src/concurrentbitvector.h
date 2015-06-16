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
			bool PutInBloomFilter(std::vector<HfPtr> & hf, F f, char farg, uint64_t hash0)
			{
				for (size_t i = 1; i < hf.size(); i++)
				{
					uint64_t idx = ((*hf[i]).*f)(farg);
					uint64_t bit;
					uint64_t element;
					GetCoord(idx, element, bit);
					filter_[hash0 + element].fetch_or(uint32_t(1) << uint32_t(bit));
				}

				return true;
			}

		template<class HfPtr, class F>
			bool IsInBloomFilter(const std::vector<HfPtr> & hf, F f, char farg, uint64_t hash0) const
			{
				for (size_t i = 1; i < hf.size(); i++)
				{
					uint64_t idx = ((*hf[i]).*f)(farg);
					uint64_t bit;
					uint64_t element;
					GetCoord(idx, element, bit);
					if ((filter_[hash0 + element] & (uint32_t(1) << uint32_t(bit))) == 0)
					{
						return false;
					}
				}

				return true;
		}
		
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