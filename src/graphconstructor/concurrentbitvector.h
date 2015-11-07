#ifndef _CONCURRENT_BIT_VECTOR_
#define _CONCURRENT_BIT_VECTOR_

#include <cstdlib>
#include <vector>
#include <atomic>

namespace Sibelia
{
	class ConcurrentBitVector
	{
	public:
		~ConcurrentBitVector();
		ConcurrentBitVector(size_t size);
		void Init();
		size_t Size() const;
		size_t GetPower() const;
		void SetConcurrently(size_t idx);
		bool Get(size_t idx) const;
	private:
		static const size_t SUCCESS = -1;
		typedef std::atomic<uint32_t> UInt;
		size_t size_;
		size_t realSize_;
		UInt * filter_;
		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const;
	};

}

#endif