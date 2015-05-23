#ifndef _CONCURRENT_BIT_VECTOR_
#define _CONCURRENT_BIT_VECTOR_

#include <vector>
#include <atomic>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace Sibelia
{
	template<class T, size_t power>
	class GenericBitVector
	{
	public:
		GenericBitVector(size_t size) : size_(size), realSize_(size / (sizeof(T)* 4) + 1), filter_(new T[realSize_])
		{
			Init();
		}

		void Init()
		{
			for (size_t i = 0; i < realSize_; i++)
			{
				filter_[i] = 0;
			}
		}

		size_t Size() const
		{
			return size_;
		}

		void DumpToStream(const std::string & fileName) const
		{
			std::ofstream out(fileName.c_str(), std::ios_base::binary);
			if (!out)
			{
				throw std::runtime_error("Can't open the temporary file");
			}

			out.write(reinterpret_cast<const char*>(filter_), realSize_ * sizeof(T));
		}

		void ReadFromStream(const std::string & fileName)
		{
			std::ifstream in(fileName.c_str(), std::ios_base::binary);
			if (!in)
			{
				throw std::runtime_error("Can't open the temporary file");
			}

			in.read(reinterpret_cast<char*>(filter_), realSize_ * sizeof(T));
		}

		virtual ~GenericBitVector()
		{
			delete[] filter_;
		}

		protected:			
			size_t size_;
			size_t realSize_;
			T * filter_;

			void GetCoord(uint64_t idx, uint64_t & element, uint64_t & bit) const
			{
				bit = idx & ((uint64_t(1) << uint64_t(power)) - 1);
				element = idx >> power;				
			}
		};

	typedef GenericBitVector<std::atomic<uint32_t>, 5> ConcurrentBase;
	class ConcurrentBitVector : public ConcurrentBase
	{
	public:
		ConcurrentBitVector(size_t size);
		void Set(size_t idx);
	};

	typedef GenericBitVector<uint32_t, 5> NonConcurrentBase;
	class NonConcurrentBitVector : public NonConcurrentBase
	{
	public:
		NonConcurrentBitVector(size_t size);
		bool Get(size_t idx) const;
	};
}

#endif