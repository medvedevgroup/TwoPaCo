#include "fastset.h"

#include <algorithm>

namespace Sibelia
{
	FastSet::FastSet(size_t vertexLength, size_t maxHashTableSize) :
		vertexLength_(vertexLength),
		maxHashTableSize_(maxHashTableSize),		
		hashTable_(maxHashTableSize, RecordHashFunction(vertexLength), RecordEquality(vertexLength))
	{
		array_.reserve(32);
	}

	void FastSet::Insert(Record record)
	{
		RecordLess less(vertexLength_);
		for (std::vector<uint64_t> & a : array_)
		{
			auto it = std::lower_bound(a.begin(), a.end(), record.GetVertex().GetBody(), less);
			if (it != a.end())
			{
				Record cand(vertexLength_, *it);
				if (cand.GetVertex() == record.GetVertex())
				{
					if (cand.GetStatus() == Record::CANDIDATE && (cand.GetPrev() != record.GetPrev() || cand.GetNext() != record.GetNext()))
					{
						cand.ChangeStatus(Record::BIFURCATION);
						*it = cand.GetBody();
					}
				
					return;
				}
			}
		}

		auto it = hashTable_.find(record.GetBody());
		if (it != hashTable_.end())
		{
			Record cand(vertexLength_, *it);
			if (cand.GetStatus() == Record::CANDIDATE && (cand.GetPrev() != record.GetPrev() || cand.GetNext() != record.GetNext()))
			{
				cand.ChangeStatus(Record::BIFURCATION);
				hashTable_.erase(it);
				hashTable_.insert(cand.GetBody());
			}
			
			return;
		}

		hashTable_.insert(record.GetBody());
		if (hashTable_.size() >= maxHashTableSize_)
		{
			if (array_.empty())
			{
				array_.push_back(std::vector<uint64_t>(hashTable_.begin(), hashTable_.end()));
				std::sort(array_.front().begin(), array_.front().end(), less);
			}
			else if (array_.size() * maxHashTableSize_ >= array_.front().size())
			{
				for (size_t i = array_.size() - 1; i > 0; --i)
				{
					array_.front().insert(array_.front().end(), array_[i].begin(), array_[i].end());
					array_.pop_back();
				}

				array_.front().insert(array_.front().end(), hashTable_.begin(), hashTable_.end());
				std::sort(array_.front().begin(), array_.front().end(), less);
				assert(std::unique(array_.front().begin(), array_.front().end(), RecordEquality(vertexLength_)) == array_.front().end());
			}
			else
			{
				array_.push_back(std::vector<uint64_t>(hashTable_.begin(), hashTable_.end()));
				std::sort(array_.back().begin(), array_.back().end(), less);
			}

			hashTable_.erase(hashTable_.begin(), hashTable_.end());
		}
	}
}