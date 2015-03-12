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
		void PutInBloomFilter(BloomFilter & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			std::vector<size_t> hf(seed.size());
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]);
				hf[i] = hvalue % filter.Size();
			}

			filter.Put(hf);
		}

		bool IsInBloomFilter(const BloomFilter & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			std::vector<size_t> hf(seed.size());
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]);
				hf[i] = hvalue % filter.Size();
			}

			return filter.Get(hf);
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

		void CandidateCheckingWorker(uint64_t low, uint64_t high, const std::vector<uint64_t> & seed, const BloomFilter & bitVector, size_t vertexLength, TaskQueue & taskQueue, ResultQueue & resultQueue)
		{
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{									
					if (task.start == Task::GAME_OVER)
					{
						while (resultQueue.write_available() == 0);
						resultQueue.push(Result(0, Result::GAME_OVER, std::vector<bool>()));
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

					char posPrev;
					char negExtend;
					DnaString negVertex = posVertex.RevComp();
					for (size_t pos = 0; pos + vertexLength - 1 < task.str.size(); pos++)
					{
						posExtend = task.str[pos + vertexLength - 1];
						posVertex.AppendBack(posExtend);
						negVertex.AppendFront(DnaString::Reverse(posExtend));
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

						posPrev = posVertex.PopFront();
						negExtend = negVertex.PopBack();
					}
					
					while (resultQueue.write_available() == 0);
					resultQueue.push(Result(task.recId, task.start, std::move(result)));
				}				
			}
		}

		void CountingWorker(uint64_t low, uint64_t high, const std::vector<uint64_t> & seed, BloomFilter & bitVector, size_t edgeLength, TaskQueue & taskQueue)
		{
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

						posEdge.PopFront();
						negEdge.PopBack();
					}
				}
			}			
		}

		typedef unsigned long BIT_TYPE;
		const size_t BITS_COUNT = sizeof(BIT_TYPE) * 8;

		std::vector<std::vector<bool> > isCandidBit;

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
									isCandidBit[t.recId][t.start + i] = t.isCandidate[i];
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
		std::cout << "Bloom filter counting..." << std::endl;

		uint64_t low = 0;
		const size_t MAX_ROUNDS = 1;
		const size_t WORKER_THREADS = 8;
		isCandidBit.clear();
		for (size_t round = 0; round < MAX_ROUNDS; round++)
		{
			size_t fastaRecords = 0;
			uint64_t high = round == MAX_ROUNDS - 1 ? UINT64_MAX : (UINT64_MAX / MAX_ROUNDS) * (round + 1);
			time_t mark = time(0);
			{
				BloomFilter bitVector(filterSize);
				std::vector<TaskQueuePtr> taskQueue;
				std::vector<ResultQueuePtr> resultQueue;
				std::vector<boost::thread> workerThread(WORKER_THREADS);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
					resultQueue.push_back(ResultQueuePtr(new ResultQueue(QUEUE_CAPACITY)));
					workerThread[i] = boost::thread(CountingWorker, low, high, boost::cref(seed), boost::ref(bitVector), edgeLength, boost::ref(*taskQueue[i]));
				}

				for (const std::string & nowFileName : fileName)
				{
					size_t record = 0;
					for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
					{
						char ch;
						fastaRecords++;
						isCandidBit.push_back(std::vector<bool>());
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

							isCandidBit.back().push_back(false);
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

				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					TaskQueuePtr & q = taskQueue[i];
					while (q->write_available() == 0);
					q->push(Task(0, Task::GAME_OVER, std::string()));
					workerThread[i].join();
				}

				std::cout << "Counting time = " << time(0) - mark << std::endl;
				std::cout << "Vertex enumeration..." << std::endl;
				mark = time(0);
				boost::thread writerThread(WriterThread, WORKER_THREADS, fastaRecords, boost::ref(resultQueue));
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateCheckingWorker, low, high, boost::cref(seed), boost::cref(bitVector), vertexLength, boost::ref(*taskQueue[i]), boost::ref(*resultQueue[i]));
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

				for (TaskQueuePtr & q : taskQueue)
				{
					while (q->write_available() == 0);
					q->push(Task(0, Task::GAME_OVER, std::string()));
				}

				writerThread.join();
			}

			std::cout << "Enumeration time = " << time(0) - mark << std::endl;
			mark = time(0);
			std::cout << "Aggregation of the results..." << std::endl;
			std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> trueBifSet(0, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
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
						BIT_TYPE buf;
						char posPrev;
						char negExtend;
						size_t bitCount = 0;
						size_t kmer = 0;
						DnaString negVertex = posVertex.RevComp();
					//	std::ifstream candid(TempFile(record).c_str(), std::ios_base::binary);
					//	candid.read(reinterpret_cast<char*>(&buf), sizeof(buf));
			//			std::bitset<BITS_COUNT> candidFlag(buf);
						if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
						{
							trueBifSet.insert(posVertex.GetBody());
						}

						for (bool go = true; go; kmer++)
						{
							if (go = parser.GetChar(posExtend))
							{
								if (isCandidBit[record][kmer])
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
							}
							else if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
							{
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

			std::cout << "Aggregation time = " << time(0) - mark << std::endl;
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