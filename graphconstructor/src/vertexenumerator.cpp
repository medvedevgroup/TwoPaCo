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
				uint64_t body = str.GetBody();
				uint64_t hash = SpookyHash::Hash64(&body, sizeof(body), 0);
				return hash;
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

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize):
		vertexSize_(vertexLength)
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
				DnaString posEdge;
				for (size_t j = 0; j < edgeLength && parser.GetChar(ch); j++)
				{
					posEdge.AppendBack(ch);
				}

				if (posEdge.GetSize() == edgeLength)
				{
					DnaString negEdge = posEdge.RevComp();
					while (true)
					{						
						PutInBloomFilter(bitVector, seed, posEdge);
						PutInBloomFilter(bitVector, seed, negEdge);					
						if (parser.GetChar(ch))
						{
							negEdge.PopBack();
							posEdge.PopFront();
							negEdge.AppendFront(DnaString::Reverse(ch));
							posEdge.AppendBack(ch);							
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
		std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> trueBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
		std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> candidateBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
		for (const std::string & nowFileName : fileName)
		{
			for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); )
			{
				char ch;
				DnaString posVertex;
				for (size_t j = 0; j < vertexLength && parser.GetChar(ch); j++)
				{
					posVertex.AppendBack(ch);
				}

				if (posVertex.GetSize() >= vertexLength)
				{
					char posPrev;
					char negExtend;
					DnaString negVertex = posVertex.RevComp();
					trueBifSet.insert(posVertex.GetBody());
					trueBifSet.insert(negVertex.GetBody());
					for (bool go = true; go; )
					{						
						bool putCandidate = false;
						bool checkCandidate = false;
						bool putTrueBif = negVertex == posVertex;
						if (!trueBifSet.count(posVertex.GetBody()) && !putTrueBif)
						{
							if (!candidateBifSet.count(posVertex.GetBody()))
							{
								size_t inCount = 0;
								size_t outCount = 0;
								for (char ch : DnaString::LITERAL)
								{
									DnaString inEdge = posVertex;
									DnaString outEdge = posVertex;
									inEdge.AppendFront(ch);
									outEdge.AppendBack(ch);
									inCount += IsInBloomFilter(bitVector, seed, inEdge) ? 1 : 0;
									outCount += IsInBloomFilter(bitVector, seed, outEdge) ? 1 : 0;
									std::string so = outEdge.ToString();
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
						
						if (parser.GetChar(ch))
						{
							posVertex.AppendBack(ch);
							negVertex.AppendFront(DnaString::Reverse(ch));
						}
						else
						{
							go = false;
							putTrueBif = true;
						}

						if (checkCandidate)
						{
							std::unordered_set<uint64_t, VertexHashFunction>::iterator it = candidateBifSet.find(posVertex.GetBody());
							DnaString candidate(vertexLength + 2, *it);
							char candExtend = candidate.GetChar(vertexLength);
							char candPrev = candidate.GetChar(vertexLength + 1);
							putTrueBif = (candPrev != posPrev) || (candExtend != posVertex.GetChar(vertexLength));
						}

						if (putTrueBif)
						{													
							DnaString insNegVertex = negVertex;
							insNegVertex.PopFront();

							trueBifSet.insert(posVertex.GetBody());
							trueBifSet.insert(insNegVertex.GetBody());

							candidateBifSet.erase(posVertex.GetBody());
							candidateBifSet.erase(insNegVertex.GetBody());
						}
						
						if (putCandidate)
						{
							DnaString posCandidate(posVertex);
							posCandidate.AppendBack(posPrev);
							candidateBifSet.insert(posCandidate.GetBody());

							DnaString negCandidate = negVertex;
							char negPrev = negCandidate.PopFront();
							negCandidate.AppendBack(negExtend);
							negCandidate.AppendBack(negPrev);
							candidateBifSet.insert(negCandidate.GetBody());

							std::string ps = posCandidate.ToString();
							std::string ns = negCandidate.ToString();
						}
						
						posPrev = posVertex.PopFront();						
						negExtend = negVertex.PopBack();
					}					
				}
			}
		}

		
		std::cout << "Passed: " << double(clock() - mark) / CLOCKS_PER_SEC  << std::endl;
		std::cout << "Vertex count = " << trueBifSet.size() << std::endl;
		std::cout << "FP count = " << candidateBifSet.size() << std::endl;
		
		for (uint64_t vertex : trueBifSet)
		{
			DnaString v(vertexLength, vertex);
			bifurcation_.push_back(v.GetBody());
		}

		std::sort(bifurcation_.begin(), bifurcation_.end());		
		std::cout << "Passed: " << double(clock() - mark) / CLOCKS_PER_SEC << std::endl;
		
	}

	size_t VertexEnumerator::GetVerticesCount() const
	{
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