#include <deque>
#include <ctime>
#include <memory>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_set>

#include "lib/SpookyV2.h"
#include "vertexenumerator.h"

namespace Sibelia
{	
	const size_t VertexEnumerator::INVALID_VERTEX = -1;

	namespace
	{
		void PutInBloomFilter(std::vector<bool> & bitVector, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			for (const uint64_t & s : seed)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), s);
				bitVector[hvalue % bitVector.size()] = true;
			}
		}

		bool IsInBloomFilter(const std::vector<bool> & bitVector, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			for (const uint64_t & s : seed)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), s);
				if (!bitVector[hvalue % bitVector.size()])
				{
					return false;
				}
			}

			return true;
		}

		class VertexHashFunction
		{
		public:
			VertexHashFunction(size_t vertexSize) : vertexSize_(vertexSize)
			{

			}

			uint64_t operator () (const uint64_t & a) const
			{
				DnaString str(vertexSize_, a);
				return SpookyHash::Hash64(&a, sizeof(a), 0);
			}
		private:
			size_t vertexSize_;
		};

		class VertexEquality
		{
		public:
			VertexEquality(size_t vertexSize) : vertexSize_(vertexSize)
			{
				
			}

			bool operator () (const uint64_t & a, const uint64_t & b) const
			{
				DnaString stra(vertexSize_, a);
				DnaString strb(vertexSize_, b);
				return stra == strb;
			}
		private:
			size_t vertexSize_;
		};

		typedef std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> BifCandidateSet;

		size_t CharIndex(char ch)
		{
			return std::find(DnaString::LITERAL.begin(), DnaString::LITERAL.end(), ch) - DnaString::LITERAL.begin();
		}
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize)
	{
		std::cout << "Filter size = " << filterSize << std::endl;
		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}

		size_t q = 3;
		std::vector<uint64_t> seed(q);
		std::generate(seed.begin(), seed.end(), rand);		
		size_t edgeLength = vertexLength + 1;
		std::vector<bool> bitVector(filterSize, false);
		std::cout << "Bloom filter counting..." << std::endl;
		for (const std::string & nowFileName: fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
			{
				char ch;
				DnaString edge;
				for (size_t j = 0; j < edgeLength && parser.GetChar(ch); j++)
				{
					edge.AppendBack(ch);
				}

				if (edge.GetSize() == edgeLength)
				{
					while (true)
					{
						PutInBloomFilter(bitVector, seed, edge);
						if (parser.GetChar(ch))
						{
							edge.AppendBack(ch);
							edge.PopFront();
						}
						else
						{
							break;
						}
					}
				}
			}	
		}

		size_t mark = clock();
		std::cout << "Passed: " << double(clock()) / CLOCKS_PER_SEC << std::endl;
		std::cout << "Vertex enumeration..." << std::endl;
		std::unordered_set<uint64_t, VertexHashFunction> candidateBifSet(0, VertexHashFunction(vertexLength));
		std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> trueBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));		
		for (const std::string & nowFileName : fileName)
		{
			bool start = true;
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); start = true)
			{
				char ch;
				DnaString vertex;
				for (size_t j = 0; j < vertexLength && parser.GetChar(ch); j++)
				{
					vertex.AppendBack(ch);
				}

				if (vertex.GetSize() >= vertexLength)
				{
					for (bool go = true; go; start = false)
					{
						bool putTrueBif = false;
						bool putCandidate = false;
						bool checkCandidate = false;
						if (start)
						{
							putTrueBif = true;												
						}
						else
						{
							if (!trueBifSet.count(vertex.GetBody()))
							{
								if (!candidateBifSet.count(vertex.GetBody()))
								{
									size_t inCount = 0;
									size_t outCount = 0;
									for (char ch : DnaString::LITERAL)
									{
										DnaString inEdge = vertex;
										DnaString outEdge = vertex;
										inEdge.AppendFront(ch);
										outEdge.AppendBack(ch);
										inCount += IsInBloomFilter(bitVector, seed, inEdge) ? 1 : 0;
										outCount += IsInBloomFilter(bitVector, seed, outEdge) ? 1 : 0;
										if (inCount > 1 || outCount > 1)
										{
											putCandidate = true;
											break;
										}
									}									
								}
								else
								{
									checkCandidate = true;
								}								
							}							
						}

						if (parser.GetChar(ch))
						{
							vertex.AppendBack(ch);							
						}
						else
						{
							go = false;
							putTrueBif = true;
						}

						if (putTrueBif)
						{
							trueBifSet.insert(vertex.GetBody());
							candidateBifSet.erase(vertex.GetBody());
						}
						else if (putCandidate)
						{

						}
						else if (checkCandidate)
						{

						}
						
						vertex.PopBack();
					}					
				}
			}
		}

		
		std::cout << "Passed: " << double(clock() - mark) / CLOCKS_PER_SEC  << std::endl;
	//	std::cout << "Vertex count = " << bifSet_.size() << std::endl;
		/*
		std::cout << "Sorting and duplicates removal..." << std::endl;
		mark = clock();
		//bifurcation_.assign(bifSet_.begin(), bifSet_.end());
		std::sort(bifurcation_.begin(), bifurcation_.end());		
		std::cout << "Passed: " << double(clock() - mark) / CLOCKS_PER_SEC << std::endl;*/
		
	}

	size_t VertexEnumerator::GetVerticesCount() const
	{
		//return bifSet_.size();
		return bifurcation_.size();
	}

	size_t VertexEnumerator::GetId(const DnaString & vertex) const
	{
		std::vector<uint64_t>::const_iterator it = std::lower_bound(bifurcation_.begin(), bifurcation_.end(), vertex.GetBody());
		if (it == bifurcation_.end() || *it != vertex.GetBody())
		{
			return INVALID_VERTEX;
		}

		return it - bifurcation_.begin();

	}
}