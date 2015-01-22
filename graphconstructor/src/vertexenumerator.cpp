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
			uint64_t operator () (const uint64_t & a) const
			{
				 return SpookyHash::Hash64(&a, sizeof(a), 0);
			}
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

		class BifurcationCandidate
		{
		public:
			
		private:
			DnaString body;
		};
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize)
	{
		if (vertexLength > 29)
		{
			throw std::runtime_error("The vertex size is too large");
		}

		size_t q = 3;
		std::vector<uint64_t> seed(q);
		std::generate(seed.begin(), seed.end(), rand);

		size_t edgeLength = vertexLength + 1;
		std::vector<bool> bitVector(filterSize, false);
		std::cerr << "Bloom filter counting..." << std::endl;
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
		std::cerr << "Passed: " << double(clock()) / CLOCKS_PER_SEC << std::endl;
		std::cerr << "Vertex enumeration...";
		
		BifCandidateSet bifSet(0, VertexHashFunction(), VertexEquality(vertexLength));
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
					while (true)
					{
						if (start)
						{
							bifurcation_.push_back(vertex.GetBody());
						}
						else
						{
							bool isBifurcation = false;							
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
								    isBifurcation = true;
								    break;
								}
							}

							if (isBifurcation)
							{
								bifurcation_.push_back(vertex.GetBody());
							}
						}

						if (parser.GetChar(ch))
						{
							vertex.AppendBack(ch);
							vertex.PopFront();
							start = false;
						}
						else
						{
							bifurcation_.push_back(vertex.GetBody());
							break;
						}						
					}					
				}
			}
		}

		std::cout << "Passed: " << double(clock() - mark) / CLOCKS_PER_SEC  << std::endl;
		std::cout << "Sorting and duplicates removal..." << std::endl;
		mark = clock();
		std::sort(bifurcation_.begin(), bifurcation_.end());
		bifurcation_.erase(std::unique(bifurcation_.begin(), bifurcation_.end()), bifurcation_.end());
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