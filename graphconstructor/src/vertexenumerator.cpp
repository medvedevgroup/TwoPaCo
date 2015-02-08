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
				bool x = stra == strb;
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
					while (true)
					{						
						PutInBloomFilter(bitVector, seed, posEdge);
						if (parser.GetChar(ch))
						{
							posEdge.PopFront();
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
					for (bool start = true; ; start = false)
					{												
						if (parser.GetChar(ch))
						{							
							if (!start)
							{								
								char prev[] = { posPrev, DnaString::Reverse(ch) }; 
								char extend[] = { ch, negExtend };
								DnaString vertex[] = { posVertex, negVertex };
								for (size_t strand = 0; strand < 2; strand++)
								{
									if (trueBifSet.count(vertex[strand].GetBody()) == 0)
									{
										if (!candidateBifSet.count(vertex[strand].GetBody()))
										{
											size_t inCount = 0;
											size_t outCount = 0;
											for (char nextCh : DnaString::LITERAL)
											{
												DnaString inEdge = vertex[strand];
												DnaString outEdge = vertex[strand];
												inEdge.AppendFront(nextCh);
												outEdge.AppendBack(nextCh);
												inCount += IsInBloomFilter(bitVector, seed, inEdge) ? 1 : 0;
												outCount += IsInBloomFilter(bitVector, seed, outEdge) ? 1 : 0;
												if (inCount > 1 || outCount > 1)
												{
													DnaString candidate(vertex[strand]);
													candidate.AppendBack(extend[strand]);
													candidate.AppendBack(prev[strand]);													
													candidateBifSet.insert(candidate.GetBody());
													break;
												}
											}
										}
										else
										{											
											std::unordered_set<uint64_t, VertexHashFunction>::iterator it = candidateBifSet.find(vertex[strand].GetBody());
											DnaString candidate(vertexLength + 2, *it);
											char candExtend = candidate.GetChar(vertexLength);
											char candPrev = candidate.GetChar(vertexLength + 1);
											if ((candPrev != prev[strand]) || (candExtend != extend[strand]))
											{
												trueBifSet.insert(posVertex.GetBody());
												candidateBifSet.erase(posVertex.GetBody());
											}
										}
									}
								}
							}

							posVertex.AppendBack(ch);
							negVertex.AppendFront(DnaString::Reverse(ch));
							posPrev = posVertex.PopFront();
							negExtend = negVertex.PopBack();
						}
						else
						{							
							trueBifSet.insert(posVertex.GetBody());
							break;
						}
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
		DnaString check[2] = { vertex, vertex.RevComp() };
		for (DnaString str : check)
		{
			std::vector<uint64_t>::const_iterator it = std::lower_bound(bifurcation_.begin(), bifurcation_.end(), str.GetBody());
			if (it != bifurcation_.end() || *it == str.GetBody())
			{
				return it - bifurcation_.begin();
			}
		}

		return INVALID_VERTEX;
	}
}