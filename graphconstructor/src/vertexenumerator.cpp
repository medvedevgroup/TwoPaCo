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
		uint64_t ChopHash(uint64_t hash, uint64_t filterSize)
		{
			return hash & ((uint64_t(1) << filterSize) - 1);
		}

		void PutInBloomFilter(ConcurrentBitVector & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{			
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]);	
	
				filter.SetConcurrently(ChopHash(hvalue, filter.GetPower()));
			}
		}

		bool IsInBloomFilter(const ConcurrentBitVector & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]);
				if (!filter.Get(ChopHash(hvalue, filter.GetPower())))
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
			size_t recId;
			uint64_t start;			
			std::string str;
			static const size_t TASK_SIZE = 65536;
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(size_t recId, uint64_t start, std::string && str) : recId(recId), start(start), str(std::move(str)) {}
		};

		struct Result
		{
			size_t recId;
			uint64_t start;
			std::vector<bool> isCandidate;
			static const size_t GAME_OVER = SIZE_MAX;
			Result() {}
			Result(size_t recId, uint64_t start, std::vector<bool> && isCandidate) : recId(recId), start(start), isCandidate(std::move(isCandidate)) {}
			bool operator < (const Result & res) const
			{
				return start < res.start;
			}
		};

		const size_t QUEUE_CAPACITY = 12;
		typedef boost::lockfree::spsc_queue<Task> TaskQueue;
		typedef std::unique_ptr<TaskQueue> TaskQueuePtr;
		typedef boost::lockfree::spsc_queue<Result> ResultQueue;
		typedef std::unique_ptr<ResultQueue> ResultQueuePtr;

		std::string TempFile(size_t record)
		{
			std::stringstream ss;
			ss << record << ".bin";
			return ss.str();
		}

		void CandidateCheckingWorker(uint64_t low, uint64_t high, const std::vector<uint64_t> & seed, const ConcurrentBitVector & bitVector, size_t vertexLength, TaskQueue & taskQueue, std::vector<std::unique_ptr<ConcurrentBitVector> > & isCandidBit)
		{
			std::vector<size_t> hf(seed.size());
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{									
					if (task.start == Task::GAME_OVER)
					{						
						break;
					}

					if (task.str.size() < vertexLength)
					{
						continue;
					}

					char posExtend;
					DnaString posVertex;
					std::vector<bool> result(task.str.size() - vertexLength + 1);
					for (size_t j = 0; j < vertexLength - 1; j++)
					{
						posVertex.AppendBack(task.str[j]);
					}
					
					char negExtend;
					char posPrev = 0;
					DnaString negVertex = posVertex.RevComp();
					for (size_t pos = 0; pos + vertexLength - 1 < task.str.size(); pos++)
					{						
						posVertex.AppendBack(task.str[pos + vertexLength - 1]);
						negVertex.AppendFront(DnaString::Reverse(task.str[pos + vertexLength - 1]));
						posExtend = pos + vertexLength < task.str.size() ? task.str[pos + vertexLength] : 0;
						assert(posVertex.RevComp() == negVertex);
	
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
								char revNextCh = DnaString::Reverse(nextCh);								
								if (posPrev != 0 && nextCh == posPrev)
								{
									++inCount;
								}
								else
								{
									DnaString posInEdge = posVertex;
									DnaString negInEdge = negVertex;
									posInEdge.AppendFront(nextCh);
									negInEdge.AppendBack(DnaString::Reverse(nextCh));
									if (IsInBloomFilter(bitVector, seed, posInEdge) || IsInBloomFilter(bitVector, seed, negInEdge))
									{
										inCount++;
									}									
								}

								if (posExtend != 0 && nextCh == posExtend)
								{
									++outCount;
								}
								else
								{
									DnaString posOutEdge = posVertex;
									DnaString negOutEdge = negVertex;
									posOutEdge.AppendBack(nextCh);
									negOutEdge.AppendFront(DnaString::Reverse(nextCh));
									if (IsInBloomFilter(bitVector, seed, posOutEdge) || IsInBloomFilter(bitVector, seed, negOutEdge))
									{
										outCount++;
									}
								}
							}

							if (inCount > 1 || outCount > 1)
							{
								isCandidBit[task.recId]->SetConcurrently(task.start + pos);
							}
						}

						posPrev = posVertex.PopFront();
						negExtend = negVertex.PopBack();
					}
				}				
			}
		}

		void CountingWorker(uint64_t low, uint64_t high, const std::vector<uint64_t> & seed, ConcurrentBitVector & bitVector, size_t edgeLength, TaskQueue & taskQueue)
		{
			std::vector<size_t> hf(seed.size());
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < edgeLength)
					{
						continue;
					}

					char ch;
					DnaString posEdge;
					for (size_t j = 0; j < edgeLength - 1; j++)
					{
						posEdge.AppendBack(task.str[j]);
					}

					DnaString negEdge = posEdge.RevComp();
					for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); pos++)
					{
						ch = task.str[pos + edgeLength - 1];
						posEdge.AppendBack(ch);
						negEdge.AppendFront(DnaString::Reverse(ch));
						assert(posEdge.RevComp() == negEdge);

						size_t k = 0;
						size_t hit = 0;
						DnaString kmer[2][2] = { { posEdge, negEdge }, { posEdge, negEdge } };
						for (size_t i = 0; i < 2 && !hit; i++, k++)
						{
							kmer[i][k].PopBack();
							kmer[i][1 - k].PopFront();
							uint64_t hvalue = UINT64_MAX;
							assert(kmer[i][0] == kmer[i][1].RevComp());
							for (size_t j = 0; j < 2 && !hit; j++)
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

						posEdge.PopFront();
						negEdge.PopBack();
					}
				}
			}			
		}

		typedef unsigned long BIT_TYPE;
		const size_t BITS_COUNT = sizeof(BIT_TYPE) * 8;

		

		void WriterThread(size_t workerThreads, size_t fastaRecords, std::vector<ResultQueuePtr> & resultQueue)
		{
			//std::ofstream log("log");
			std::ostream & log = std::cerr;			
			std::vector<uint64_t> bitCount(fastaRecords, 0);
			std::vector<uint64_t> awaitStart(fastaRecords, 0);
			std::vector<std::ofstream> candid(fastaRecords);
			std::vector<std::set<Result> > buffer(fastaRecords);			
			std::vector<std::bitset<BITS_COUNT> > bits(fastaRecords);
			/*
			for (size_t i = 0; i < fastaRecords; i++)
			{
				std::string name = TempFile(i).c_str();
				candid[i].open(name.c_str(), std::ios_base::binary);
			}*/
			
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
							size_t record = res.recId;
							buffer[record].insert(res);							
							while (buffer[record].size() > 0 && buffer[record].begin()->start == awaitStart[record])
							{								
								awaitStart[record] += buffer[record].begin()->isCandidate.size();
								/*
								for (bool value : buffer[record].begin()->isCandidate)
								{
									bits[record].set(bitCount[record]++, value);
									if (bitCount[record] == bits[record].size())
									{
										bitCount[record] = 0;
										BIT_TYPE buf = bits[record].to_ulong();
										candid[record].write(reinterpret_cast<const char*>(&buf), sizeof(BIT_TYPE) / sizeof(char));
									}
								}*/

								const Result & t = *(buffer[record].begin());
								for (size_t i = 0; i < t.isCandidate.size(); i++)
								{
		//							isCandidBit[t.recId][t.start + i] = t.isCandidate[i];
								}

								buffer[record].erase(buffer[record].begin());
							}
						}
					}
				}							
			}

			/*
			for (size_t record = 0; record < fastaRecords; record++)
			{
				assert(buffer[record].empty());
				if (bitCount[record] > 0)
				{
					BIT_TYPE buf = bits[record].to_ulong();
					candid[record].write(reinterpret_cast<const char*>(&buf), sizeof(BIT_TYPE) / sizeof(char));
				}
			}*/			
		}
	}


	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, size_t hashFunctions, size_t rounds, size_t threads) :
		vertexSize_(vertexLength)
	{
		std::cout << "Threads = " << threads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << filterSize << std::endl;
		std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> trueBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));

		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}
		
		std::vector<uint64_t> seed(hashFunctions);
		std::generate(seed.begin(), seed.end(), rand);		
		size_t edgeLength = vertexLength + 1;

		uint64_t low = 0;		
		for (size_t round = 0; round < rounds; round++)
		{
			std::vector<size_t> fastaRecordsSize;
			std::vector<std::unique_ptr<ConcurrentBitVector> > isCandidBit;
			uint64_t high = round == rounds - 1 ? UINT64_MAX : (UINT64_MAX / rounds) * (round + 1);
			time_t mark = time(0);
			{
				ConcurrentBitVector bitVector(size_t(1) << filterSize, filterSize);
				std::vector<TaskQueuePtr> taskQueue;
				std::vector<boost::thread> workerThread(threads);
				std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
				std::cout << "Counting\tEnumeration\tAggregation" << std::endl;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
					workerThread[i] = boost::thread(CountingWorker, low, high, boost::cref(seed), boost::ref(bitVector), edgeLength, boost::ref(*taskQueue[i]));
				}

				for (const std::string & nowFileName : fileName)
				{
					size_t record = 0;
					for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
					{
						char ch;
						fastaRecordsSize.push_back(0);						
						std::string buf;
						uint64_t prev = 0;
						uint64_t start = 0;
						bool over = false;
						do
						{
							over = !parser.GetChar(ch);
							if (!over)
							{
								start++;
								buf.push_back(ch);
							}

							++fastaRecordsSize.back();
							if (buf.size() >= edgeLength && (buf.size() == Task::TASK_SIZE || over))
							{
								for (bool found = false; !found;)
								{
									for (TaskQueuePtr & q : taskQueue)
									{
										if (q->write_available() > 0)
										{
											std::string overlap;
											if (!over)
											{
												overlap.assign(buf.end() - edgeLength + 1, buf.end());
											}

											q->push(Task(record, prev, std::move(buf)));
											prev = start - edgeLength + 1;
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

				for (size_t sz : fastaRecordsSize)
				{
					isCandidBit.push_back(std::unique_ptr<ConcurrentBitVector>(new ConcurrentBitVector(sz)));
				}

				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					TaskQueuePtr & q = taskQueue[i];
					while (q->write_available() == 0);
					q->push(Task(0, Task::GAME_OVER, std::string()));
					workerThread[i].join();
				}

				std::cout << time(0) - mark << "\t";
				mark = time(0);			
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateCheckingWorker, low, high, boost::cref(seed), boost::cref(bitVector), vertexLength, boost::ref(*taskQueue[i]), boost::ref(isCandidBit));
				}

				for (const std::string & nowFileName : fileName)
				{
					size_t record = 0;
					for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
					{
						char ch;
						std::string buf;
						uint64_t prev = 0;
						uint64_t start = 0;
						bool over = false;
						do
						{
							over = !parser.GetChar(ch);
							if (!over)
							{
								start++;
								buf.push_back(ch);
							}

							if (buf.size() >= vertexLength && (buf.size() == Task::TASK_SIZE || over))
							{
								for (bool found = false; !found;)
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

											q->push(Task(record, prev, std::move(buf)));
											prev = start - vertexLength + 1;
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

				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					while (taskQueue[i]->write_available() == 0);
					taskQueue[i]->push(Task(0, Task::GAME_OVER, std::string()));
					workerThread[i].join();
				}
			}

			std::cout << time(0) - mark << "\t";
			mark = time(0);
			size_t nowCount = 0;
			std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> candidateBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
			for (const std::string & nowFileName : fileName)
			{
				size_t record = 0;
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
				{
					char posExtend;
					DnaString posVertex;
					for (size_t j = 0; j < vertexLength && parser.GetChar(posExtend); j++)
					{
						posVertex.AppendBack(posExtend);
					}

					if (posVertex.GetSize() >= vertexLength)
					{
						char posPrev;
						char negExtend;
						size_t bitCount = 0;
						size_t kmer = 0;
						DnaString negVertex = posVertex.RevComp();
					//	std::ifstream candid(TempFile(record).c_str(), std::ios_base::binary);
					//	candid.read(reinterpret_cast<char*>(&buf), sizeof(buf));
					//	std::bitset<BITS_COUNT> candidFlag(buf);
						if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
						{
							++nowCount;
							trueBifSet.insert(posVertex.GetBody());
						}

						for (bool go = true; go; kmer++)
						{
							if (go = parser.GetChar(posExtend))
							{
								if (isCandidBit[record]->Get(kmer))
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
												++nowCount;
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
													++nowCount;
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
							}
							else if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
							{
								++nowCount;
								trueBifSet.insert(posVertex.GetBody());
							}
							/*
							if (++bitCount >= BITS_COUNT)
							{								
								candid.read(reinterpret_cast<char*>(&buf), sizeof(buf));
								candidFlag = buf;
								bitCount -= BITS_COUNT;
							}*/
						}
					}
				}	
			}

			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << nowCount << std::endl;
			std::cout << "FP count = " << candidateBifSet.size() << std::endl;
			std::cout << std::string(80, '-') << std::endl;
			low = high + 1;
		}
		
		for (uint64_t vertex : trueBifSet)
		{
			DnaString v(vertexLength, vertex);
			bifurcation_.push_back(v.GetBody());
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