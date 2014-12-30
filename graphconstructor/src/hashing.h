
#ifndef _HASHING_H_
#define _HASHING_H_

namespace Sibelia
{
	template<class Iterator>
	class SlidingWindow
	{
	public:
		SlidingWindow() {}
		SlidingWindow(Iterator kMerStart, Iterator end, size_t k) : k_(k), highPow_(1),
			kMerStart_(kMerStart), kMerEnd_(AdvanceForward(kMerStart, k - 1)), end_(end)
		{
			for (size_t i = 1; i < k; i++)
			{
				highPow_ *= HASH_BASE;
			}

			value_ = CalcKMerHash(kMerStart, k);
		}

		uint64_t GetValue() const
		{
			return value_;
		}

		uint64_t GetK() const
		{
			return k_;
		}

		Iterator GetBegin() const
		{
			return kMerStart_;
		}

		Iterator GetEnd() const
		{
			DNASequence::StrandIterator ret = kMerEnd_;
			return ++ret;
		}

		bool Move()
		{
			value_ = (value_ - *kMerStart_ * highPow_) * HASH_BASE;
			++kMerStart_;
			++kMerEnd_;
			if (Valid())
			{
				value_ += *kMerEnd_;
				assert(value_ == CalcKMerHash(kMerStart_, k_));
				return true;
			}

			return false;
		}

		bool Valid() const
		{
			return kMerEnd_ != end_;
		}

		static uint64_t CalcKMerHash(Iterator it, uint64_t k)
		{
			uint64_t base = 1;
			uint64_t hash = 0;
			std::advance(it, k - 1);
			for (size_t i = 0; i < k; i++)
			{
				hash += *it * base;
				base *= HASH_BASE;
				if (i != k - 1)
				{
					--it;
				}
			}

			return hash;
		}

	private:
		uint64_t k_;
		uint64_t highPow_;
		Iterator kMerStart_;
		Iterator kMerEnd_;
		Iterator end_;
		uint64_t value_;
	};
}

#endif