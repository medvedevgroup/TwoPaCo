#ifndef _SORT_SUPPORT_H_
#define _SORT_SUPPORT_H_

#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>

typedef std::vector<uint64_t>::iterator It;

class ObjectValue;

class ObjectReference
{
public:
	ObjectReference() : recordSize_(0) {}
	ObjectReference(It ptr, size_t recordSize) : ptr_(ptr), recordSize_(recordSize) {}

	void operator = (ObjectReference source) const
	{
		std::copy(source.ptr_, source.ptr_ + recordSize_, ptr_);
	}

	void operator = (const ObjectValue & source) const;

	It GetIterator() const
	{
		return ptr_;
	}

	size_t GetRecordSize() const
	{
		return recordSize_;
	}

private:
	It ptr_;
	size_t recordSize_;
};

class ObjectValue
{
public:
	ObjectValue() {}
	ObjectValue(ObjectReference prx) : object_(prx.GetIterator(), prx.GetIterator() + prx.GetRecordSize()) {}
	ObjectValue(It ptr, size_t recordSize) : object_(ptr, ptr + recordSize) {}
	std::vector<uint64_t>::const_iterator GetIterator() const
	{
		return object_.begin();
	}

private:
	std::vector<uint64_t> object_;
};

void ObjectReference::operator = (const ObjectValue & source) const
{
	std::copy(source.GetIterator(), source.GetIterator() + recordSize_, ptr_);
}

template<class Cmp>
class Comparator
{
public:
	Comparator() {}
	Comparator(Cmp cmp) : cmp_(cmp) {}

	bool operator () (const ObjectReference & a, const ObjectReference & b) const
	{
		return cmp_(&*a.GetIterator(), &*b.GetIterator());
	}

	bool operator () (const ObjectValue & a, const ObjectReference & b) const
	{
		return cmp_(&*a.GetIterator(), &*b.GetIterator());
	}

	bool operator () (const ObjectReference & a, const ObjectValue & b) const
	{
		return cmp_(&*a.GetIterator(), &*b.GetIterator());
	}

	bool operator () (const ObjectValue & a, const ObjectValue & b) const
	{
		return cmp_(&*a.GetIterator(), &*b.GetIterator());
	}

private:
	Cmp cmp_;
};

class RecordIterator : public std::iterator<std::random_access_iterator_tag, ObjectValue, size_t, RecordIterator, ObjectReference>
{
public:
	RecordIterator() : recordSize_(0) {}
	RecordIterator(It ptr, size_t recordSize) : ptr_(ptr), recordSize_(recordSize) {}
	ObjectReference operator * () const
	{
		return ObjectReference(ptr_, recordSize_);
	}

	ObjectReference operator [] (size_t diff) const
	{
		return *(*this + diff);
	}

	It GetIterator() const
	{
		return ptr_;
	}

	size_t GetRecordSize() const
	{
		return recordSize_;
	}

	RecordIterator& operator ++()
	{
		ptr_ += recordSize_;
		return *this;
	}

	RecordIterator& operator --()
	{
		ptr_ -= recordSize_;
		return *this;
	}

	RecordIterator operator ++(int)
	{
		RecordIterator ret = *this;
		ptr_ += recordSize_;
		return ret;
	}

	RecordIterator operator --(int)
	{
		RecordIterator ret = *this;
		ptr_ -= recordSize_;
		return ret;
	}

	friend bool operator < (RecordIterator it1, RecordIterator it2);
	friend bool operator > (RecordIterator it1, RecordIterator it2);
	friend bool operator == (RecordIterator it1, RecordIterator it2);
	friend bool operator != (RecordIterator it1, RecordIterator it2);
	friend size_t operator - (RecordIterator it1, RecordIterator it2);
	friend RecordIterator operator - (RecordIterator it1, size_t shift);
	friend RecordIterator operator + (RecordIterator it1, size_t shift);

private:
	It ptr_;
	size_t recordSize_;
};

bool operator < (RecordIterator it1, RecordIterator it2)
{
	return it1.ptr_ < it2.ptr_;
}

bool operator >(RecordIterator it1, RecordIterator it2)
{
	return it1.ptr_ > it2.ptr_;
}

bool operator == (RecordIterator it1, RecordIterator it2)
{
	return it1.ptr_ == it2.ptr_;
}

bool operator != (RecordIterator it1, RecordIterator it2)
{
	return !(it1 == it2);
}

RecordIterator operator - (RecordIterator it1, size_t shift)
{
	return RecordIterator(it1.ptr_ - shift * it1.recordSize_, it1.recordSize_);
}

RecordIterator operator + (RecordIterator it1, size_t shift)
{
	return RecordIterator(it1.ptr_ + shift * it1.recordSize_, it1.recordSize_);
}

size_t operator - (RecordIterator it1, RecordIterator it2)
{
	return (it1.ptr_ - it2.ptr_) / it1.recordSize_;
}

namespace std
{
	template<>
	void swap(ObjectReference & it1, ObjectReference & it2)
	{
		ObjectValue buf(it1.GetIterator(), it1.GetRecordSize());
		std::copy(it2.GetIterator(), it2.GetIterator() + it2.GetRecordSize(), it1.GetIterator());
		std::copy(buf.GetIterator(), buf.GetIterator() + it1.GetRecordSize(), it2.GetIterator());
	}
}

#endif