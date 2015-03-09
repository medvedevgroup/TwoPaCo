#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_set>

#include <boost/ref.hpp>
#include <boost/thread.hpp>
#include <boost/lockfree/spsc_queue.hpp>

#include "lib/SpookyV2.h"
#include "ngramhashing/cyclichash.h"

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

		struct Task
		{
			uint64_t start;
			std::string str;
			static const size_t TASK_SIZE = 64;
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(uint64_t start, std::string && str) : start(start), str(std::move(str)) {}
		};

		struct Result
		{
			uint64_t start;
			std::vector<bool> isCandidate;
			static const size_t GAME_OVER = SIZE_MAX;
			Result() {}
			Result(uint64_t start, std::vector<bool> && isCandidate) : start(start), isCandidate(std::move(isCandidate)) {}
			bool operator < (const Result & res) const
			{
				return start < res.start;
			}
		};

		const std::string & TEMP_FILE = "cand.bin";
		
		const size_t QUEUE_CAPACITY = 1;
		typedef boost::lockfree::spsc_queue<Task> TaskQueue;
		typedef std::unique_ptr<TaskQueue> TaskQueuePtr;
		typedef boost::lockfree::spsc_queue<Result> ResultQueue;
		typedef std::unique_ptr<ResultQueue> ResultQueuePtr;

		void CandidateCheckingWorker(uint64_t low, uint64_t high, const std::vector<uint64_t> & seed, const std::vector<bool> & bitVector, size_t vertexLength, TaskQueue & taskQueue, ResultQueue & resultQueue)
		{
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{									
					if (task.start == Task::GAME_OVER)
					{
						while (resultQueue.write_available() == 0);
						resultQueue.push(Result(Result::GAME_OVER, std::vector<bool>()));
						break;
					}

					if (task.str.size() < vertexLength)
					{
						continue;
					}

					char posExtend;
					DnaString posVertex;
					std::vector<bool> result(task.str.size() - vertexLength + 1);
					for (size_t j = 0; j < vertexLength; j++)
					{
						posVertex.AppendBack(task.str[j]);
					}

					char posPrev;
					char negExtend;
					size_t pos = 0;
					DnaString negVertex = posVertex.RevComp();
					do
					{
						size_t hit = 0;
						uint64_t hvalue = UINT64_MAX;
						DnaString kmer[] = { posVertex, negVertex };
						for (size_t i = 0; i < 2; i++)
						{
							uint64_t body = kmer[i].GetBody();
							hvalue = std::min(hvalue, SpookyHash::Hash64(&body, sizeof(body), seed[0]));
						}

						if (hvalue >= low && hvalue <= high)
						{
							size_t inCount = 0;
							size_t outCount = 0;
							for (int i = 0; i < DnaString::LITERAL.size() && inCount < 2 && outCount < 2; i++)
							{
								char nextCh = DnaString::LITERAL[i];
								DnaString posInEdge = posVertex;
								DnaString posOutEdge = posVertex;
								posInEdge.AppendFront(nextCh);
								posOutEdge.AppendBack(nextCh);
								DnaString negInEdge = negVertex;
								DnaString negOutEdge = negVertex;
								negInEdge.AppendBack(DnaString::Reverse(nextCh));
								negOutEdge.AppendFront(DnaString::Reverse(nextCh));
								assert(posInEdge.RevComp() == negInEdge);
								assert(posOutEdge.RevComp() == negOutEdge);
								if (IsInBloomFilter(bitVector, seed, posInEdge) || IsInBloomFilter(bitVector, seed, negInEdge))
								{
									inCount++;
								}

								if (IsInBloomFilter(bitVector, seed, posOutEdge) || IsInBloomFilter(bitVector, seed, negOutEdge))
								{
									outCount++;
								}
							}

							if (inCount > 1 || outCount > 1)
							{
								result[pos] = true;
							}
						}

						posExtend = task.str[pos + vertexLength];
						posVertex.AppendBack(posExtend);
						negVertex.AppendFront(DnaString::Reverse(posExtend));
						posPrev = posVertex.PopFront();
						negExtend = negVertex.PopBack();
						pos++;
					}
					while (pos <= task.str.size() - vertexLength);
					while (resultQueue.write_available() == 0);
					resultQueue.push(Result(task.start, std::move(result)));
				}				
			}
		}

		typedef unsigned long BIT_TYPE;
		const size_t BITS_COUNT = sizeof(BIT_TYPE) * 8;

		void WriterThread(size_t workerThreads, std::vector<ResultQueuePtr> & resultQueue)
		{
			uint64_t last = 0;
			std::set<Result> buffer;
			std::ofstream candid(TEMP_FILE.c_str(), std::ios_base::binary);			
			size_t bitCount = 0;
			std::bitset<BITS_COUNT> bits;
			while (workerThreads > 0)
			{
				for (ResultQueuePtr & q : resultQueue)
				{
					Result res;					
					while (q->pop(res))
					{
						if (res.start == Result::GAME_OVER)
						{
							workerThreads--;
						}
						else
						{
							buffer.insert(res);
						}
					}
				}

				while (buffer.size() > 0 && buffer.begin()->start <= last + Task::TASK_SIZE)
				{
					last = buffer.begin()->start;
					size_t cnt = 0;
					for (bool value : buffer.begin()->isCandidate)
					{
					//	std::cout << (buffer.begin()->start + cnt++) << ": " << value << std::endl;
						bits.set(bitCount++, value);
						if (bitCount == bits.size())
						{
							bitCount = 0;
							BIT_TYPE buf = bits.to_ulong();
							candid.write(reinterpret_cast<const char*>(&buf), sizeof(BIT_TYPE) / sizeof(char));
						}
					}

					buffer.erase(buffer.begin());
				}
			}

			if (last > 0)
			{
				BIT_TYPE buf = bits.to_ulong();
				candid.write(reinterpret_cast<const char*>(&buf), sizeof(BIT_TYPE) / sizeof(char));
			}
		}		
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, size_t q) :
		vertexSize_(vertexLength)
	{
		std::cout << "Filter size = " << filterSize << std::endl;
		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}
		
		std::vector<uint64_t> seed(q);
		std::generate(seed.begin(), seed.end(), rand);		
		size_t edgeLength = vertexLength + 1;
		std::vector<bool> bitVector(filterSize, false);
		std::cout << "Bloom filter counting..." << std::endl;

		uint64_t low = 0;
		const size_t MAX_ROUNDS = 1;
		const size_t WORKER_THREADS = 2;
		for (size_t round = 0; round < MAX_ROUNDS; round++)
		{
			uint64_t high = round == MAX_ROUNDS - 1 ? UINT64_MAX : (UINT64_MAX / MAX_ROUNDS) * (round + 1);
			for (const std::string & nowFileName : fileName)
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
							size_t k = 0;
							size_t hit = 0;
							DnaString kmer[2][2] = { { posEdge, negEdge }, { posEdge, negEdge } };							
							for (size_t i = 0; i < 2; i++, k++)
							{								
								kmer[i][k].PopBack();
								kmer[i][1 - k].PopFront();
								uint64_t hvalue = UINT64_MAX;
								assert(kmer[i][0] == kmer[i][1].RevComp());								
								for (size_t j = 0; j < 2; j++)
								{
									uint64_t body = kmer[i][j].GetBody();
									hvalue = std::min(hvalue, SpookyHash::Hash64(&body, sizeof(body), seed[0]));
								}
								
								hit += (hvalue >= low && hvalue <= high) ? 1 : 0;
							}
							
							if (hit)
							{
								PutInBloomFilter(bitVector, seed, posEdge);
							}

							if (parser.GetChar(ch))
							{
								posEdge.PopFront();
								posEdge.AppendBack(ch);
								negEdge.PopBack();
								negEdge.AppendFront(DnaString::Reverse(ch));
								assert(posEdge.RevComp() == negEdge);
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
			std::cout << "Vertex enumeration..." << std::endl;

			std::vector<TaskQueuePtr> taskQueue;
			std::vector<ResultQueuePtr> resultQueue;
			std::vector<boost::thread> workerThread(WORKER_THREADS);
			boost::thread writerThread(WriterThread, WORKER_THREADS, boost::ref(resultQueue));
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
				resultQueue.push_back(ResultQueuePtr(new ResultQueue(QUEUE_CAPACITY)));
				workerThread[i] = boost::thread(CandidateCheckingWorker, low, high, boost::cref(seed), boost::cref(bitVector), vertexLength, boost::ref(*taskQueue[i]), boost::ref(*resultQueue[i]));
			}

			TaskQueue q(QUEUE_CAPACITY);
			for (const std::string & nowFileName : fileName)
			{
				size_t counter = 0;
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
				{
					char ch;					
					std::string buf;
					size_t start = counter;
					bool over = false;
					do
					{
						over = !parser.GetChar(ch);
						if (!over)
						{
							counter++;
							buf.push_back(ch);
						}

						if (buf.size() >= vertexLength && (buf.size() == Task::TASK_SIZE || over))
						{
							for(bool found = false; !found; )							
							{								
								for (TaskQueuePtr & q : taskQueue)
								{
									if (q->write_available() > 0)
									{
										std::string overlap;
										if (!over)
										{
											overlap.assign(buf.end() - vertexLength + 1, buf.end());											
										}
										
										q->push(Task(start, std::move(buf)));
										start = counter - overlap.size();
										buf.swap(overlap);										
										found = true;
										break;
									}
								}
							}
						}

					} while (!over);
				}
			}

			for (TaskQueuePtr & q : taskQueue)
			{
				while (q->write_available() == 0);
				q->push(Task(Task::GAME_OVER, std::string()));
			}

			writerThread.join();
			std::ifstream candid(TEMP_FILE.c_str(), std::ios_base::binary);
			std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> trueBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
			std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> candidateBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
			for (const std::string & nowFileName : fileName)
			{
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord();)
				{
					char posExtend;
					DnaString posVertex;
					for (size_t j = 0; j < vertexLength && parser.GetChar(posExtend); j++)
					{
						posVertex.AppendBack(posExtend);
					}

					if (posVertex.GetSize() >= vertexLength)
					{
						BIT_TYPE buf;
						char posPrev;
						char negExtend;
						size_t bitCount = 0;
						DnaString negVertex = posVertex.RevComp();
						candid.read(reinterpret_cast<char*>(&buf), sizeof(buf));
						std::bitset<BITS_COUNT> candidFlag(buf);
						if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
						{
							trueBifSet.insert(posVertex.GetBody());
						}
						

						for (bool start = true;; start = false)
						{
							if (parser.GetChar(posExtend))
							{
								bool tr = posVertex.ToString() == "CTTT" || negVertex.ToString() == "CTTT";
								if (candidFlag[bitCount++])
								{
									if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
									{
										bool posFound = candidateBifSet.count(posVertex.GetBody()) > 0;
										bool negFound = candidateBifSet.count(negVertex.GetBody()) > 0;
										if (!posFound && !negFound)
										{
											DnaString candidate(posVertex);
											candidate.AppendBack(posExtend);
											candidate.AppendBack(posPrev);
											candidateBifSet.insert(candidate.GetBody());
											if (posVertex == negVertex)
											{
												negFound = true;
											}
										}

										if (posFound)
										{
											std::unordered_set<uint64_t, VertexHashFunction>::iterator it = candidateBifSet.find(posVertex.GetBody());
											DnaString candidate(vertexLength + 2, *it);
											char candExtend = candidate.GetChar(vertexLength);
											char candPrev = candidate.GetChar(vertexLength + 1);
											if ((candPrev != posPrev) || (candExtend != posExtend))
											{
												trueBifSet.insert(posVertex.GetBody());
												candidateBifSet.erase(posVertex.GetBody());
											}
										}

										if (negFound)
										{
											std::unordered_set<uint64_t, VertexHashFunction>::iterator it = candidateBifSet.find(negVertex.GetBody());
											if (it != candidateBifSet.end())
											{
												DnaString candidate(vertexLength + 2, *it);
												char candExtend = candidate.GetChar(vertexLength);
												char candPrev = candidate.GetChar(vertexLength + 1);
												if ((candPrev != DnaString::Reverse(posExtend)) || (candExtend != negExtend))
												{
													trueBifSet.insert(posVertex.GetBody());
													candidateBifSet.erase(posVertex.GetBody());
												}
											}
										}
									}									
								}

								posVertex.AppendBack(posExtend);
								negVertex.AppendFront(DnaString::Reverse(posExtend));
								posPrev = posVertex.PopFront();
								negExtend = negVertex.PopBack();

								if (bitCount >= BITS_COUNT)
								{
									candid.read(reinterpret_cast<char*>(&buf), sizeof(buf));
									candidFlag = buf;
									bitCount = 0;
								}
							}
							else
							{
								if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
								{
									trueBifSet.insert(posVertex.GetBody());
								}

								break;
							}
						}
					}
				}	
			}

			std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
			std::cout << "Vertex count = " << trueBifSet.size() << std::endl;
			std::cout << "FP count = " << candidateBifSet.size() << std::endl;			
			for (uint64_t vertex : trueBifSet)
			{
				DnaString v(vertexLength, vertex);
				bifurcation_.push_back(v.GetBody());
			}
			

			low = high + 1;
		}
		
		std::sort(bifurcation_.begin(), bifurcation_.end());
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
			if (it != bifurcation_.end() && *it == str.GetBody())
			{
				return it - bifurcation_.begin();
			}
		}

		return INVALID_VERTEX;
	}
}