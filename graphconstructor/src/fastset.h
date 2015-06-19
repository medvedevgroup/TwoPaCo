#ifndef _FAST_SET_H_
#define _FAST_SET_H_

#include <vector>
#include <cassert>
#include <unordered_set>

#include "dnastring.h"

#include "lib/SpookyV2.h"

namespace Sibelia
{
	class FastSet
	{
	public:
		FastSet(size_t vertexLength, size_t maxHashTableSize);
		class Record
		{
		public:
			static const char CANDIDATE = 'A';
			static const char BIFURCATION = 'C';

			Record() {}
			Record(DnaString posVertex, DnaString negVertex, char posPrev, char posNext, char status = CANDIDATE)
			{
				if (posVertex.GetBody() < negVertex.GetBody())
				{
					*this = Record(posVertex, posPrev, posNext, status);
				}
				else
				{
					*this = Record(negVertex, DnaString::Reverse(posNext), DnaString::Reverse(posPrev), status);
				}
			}

			Record(DnaString canVertex, char prev, char next, char status = CANDIDATE) : vertex_(canVertex)
			{
				vertex_.AppendBack(prev);
				vertex_.AppendBack(next);
				vertex_.AppendBack(status);
			}

			Record(size_t vertexLength, uint64_t body) : vertex_(vertexLength + 3, body)
			{
				assert(vertexLength <= 29);
			}

			uint64_t GetBody() const
			{
				return vertex_.GetBody();
			}

			void ChangeStatus(char newStatus)
			{
				vertex_.SetChar(vertex_.GetSize() - 1, newStatus);
			}

			DnaString GetVertex() const
			{
				return vertex_.Prefix(vertex_.GetSize() - 3);
			}

			char GetPrev() const
			{
				return vertex_.GetChar(vertex_.GetSize() - 3);
			}

			char GetNext() const
			{
				return vertex_.GetChar(vertex_.GetSize() - 2);
			}

			char GetStatus() const
			{
				return vertex_.GetChar(vertex_.GetSize() - 1);
			}

		private:
			DnaString vertex_;
		};

		void Insert(Record record);
		template<class It>
			void DumpBifurcations(It begin, size_t & falsePositives) const
			{
				falsePositives = 0;
				for (const std::vector<uint64_t> & a : array_)
				{
					for (uint64_t body : a)
					{
						Record record(vertexLength_, body);
						if (record.GetStatus() == Record::BIFURCATION)
						{
							*begin++ = record.GetVertex().GetBody();
						}
						else
						{
							++falsePositives;
						}
					}
				}

				for (uint64_t body : hashTable_)
				{
					Record record(vertexLength_, body);
					if (record.GetStatus() == Record::BIFURCATION)
					{
						*begin++ = record.GetVertex().GetBody();
					}
					else
					{
						++falsePositives;
					}
				}
			}

	private:		

		class RecordHashFunction
		{
		public:
			RecordHashFunction(size_t vertexLength) : vertexLength_(vertexLength)
			{

			}

			uint64_t operator()(const uint64_t & body) const
			{
				Record record(vertexLength_, body);
				uint64_t vertexBody = record.GetVertex().GetBody();
				return SpookyHash::Hash64(&vertexBody, sizeof(vertexBody), 0);
			}

		private:
			size_t vertexLength_;
		};

		class RecordEquality
		{
		public:
			RecordEquality(size_t vertexLength) : vertexLength_(vertexLength)
			{

			}

			bool operator() (const uint64_t & body1, const uint64_t & body2) const
			{
				Record record1(vertexLength_, body1);
				Record record2(vertexLength_, body2);
				return record1.GetVertex() == record2.GetVertex();
			}

		private:
			size_t vertexLength_;
		};

		class RecordLess
		{
		public:
			RecordLess(size_t vertexLength_) : vertexLength_(vertexLength_)
			{

			}

			bool operator() (const uint64_t & body1, const uint64_t & body2) const
			{
				Record record1(vertexLength_, body1);
				Record record2(vertexLength_, body2);
				return record1.GetVertex().GetBody() < record2.GetVertex().GetBody();
			}

		private:
			size_t vertexLength_;
		};

		size_t vertexLength_;
		size_t maxHashTableSize_;		
		std::vector<std::vector<uint64_t> > array_;
		std::unordered_set<uint64_t, RecordHashFunction, RecordEquality> hashTable_;
	};
}

#endif