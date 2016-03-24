#ifndef _CONCURRENT_BIT_VECTOR_
#define _CONCURRENT_BIT_VECTOR_

#include <cstdlib>
#include <vector>
#include <atomic>

namespace TwoPaCo
{
	class ConcurrentBitVector
	{
	public:		
		typedef uint16_t BASIC_TYPE;
		~ConcurrentBitVector();
		ConcurrentBitVector(size_t size);
		void Reset();
		size_t Size() const;
		BASIC_TYPE GetElement(size_t idx) const;
		static int PopCount(BASIC_TYPE element)
		{
			return bitCount_[element];
		}

		void OrElementCouncerrently(size_t idx, BASIC_TYPE bit);
		void SetBitConcurrently(size_t idx);
		bool GetBit(size_t idx) const;		
		void WriteToFile(const std::string & fileName) const;
		void ReadFromFile(const std::string & fileName, bool cleanUp);
	private:
		static const size_t SUCCESS = -1;		
		typedef std::atomic<BASIC_TYPE> UInt;
		static const uint64_t BASIC_TYPE_POWER = 4;		
		static const uint64_t BASIC_TYPE_BITS = sizeof(BASIC_TYPE) * 8;
		static char bitCount_[1 << sizeof(BASIC_TYPE_BITS)];
		size_t size_;
		size_t realSize_;
		UInt * filter_;
		void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const;
	};
}

#endif