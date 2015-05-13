#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <numeric>
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
		void PutInBloomFilter(ConcurrentBitVector & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]);
				filter.SetConcurrently(hvalue % filter.Size());
			}
		}

		bool IsInBloomFilter(const ConcurrentBitVector & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]);
				if (!filter.Get(hvalue % filter.Size()))
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
			bool isFinal;
			size_t recId;
			uint64_t start;
			std::string str;
			static const size_t TASK_SIZE = 65536;
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(size_t recId, uint64_t start, bool isFinal, std::string && str) : recId(recId), start(start), isFinal(isFinal), str(std::move(str)) {}
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

		const char CHAR_GAME_OVER = 0;
		const char CHAR_SEQ_START = 1;
		const char CHAR_SEQ_END = 2;
		
		typedef std::shared_ptr<std::string> StringPtr;
		typedef boost::lockfree::spsc_queue<StringPtr> StringQueue;
		typedef std::unique_ptr<StringQueue> StringQueuePtr;

		typedef std::unordered_set<uint64_t, VertexHashFunction, VertexEquality> VertexSet;

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

	}

	uint64_t NormHash(const std::vector<uint64_t> & seed, DnaString posVertex, DnaString negVertex)
	{
		uint64_t hvalue = UINT64_MAX;
		DnaString kmer[] = { posVertex, negVertex };
		for (size_t i = 0; i < 2; i++)
		{
			uint64_t body = kmer[i].GetBody();
			hvalue = std::min(hvalue, SpookyHash::Hash64(&body, sizeof(body), seed[0]));
		}

		return hvalue;
	}

	bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
	{
		return hvalue >= low && hvalue <= high;
	}

	void AggregationWorker(size_t vertexLength,
		StringQueue & charQueue,
		uint64_t low,
		uint64_t high,
		const std::vector<uint64_t> & seed,
		std::vector<std::unique_ptr<ConcurrentBitVector> > & isCandidBit,
		VertexSet & trueBifSet,		
		size_t & falseCount)
	{
		std::string tmp;
		bool start = true;
		size_t pos = 0;
		size_t record = 0;
		DnaString posVertex;
		DnaString negVertex;
		char posPrev;
		char posExtend;
		char negPrev;
		char negExtend;
		VertexSet candidateBifSet(1024, VertexHashFunction(vertexLength), VertexEquality(vertexLength));
		StringPtr chunk;
		size_t chunkPos = 0;
		for (bool go = true; go; )
		{
			if (chunk == 0 || chunkPos == chunk->size())
			{
				chunkPos = 0;
				chunk.reset();				
				while (!charQueue.pop(chunk));				
			}
	
			switch (posExtend = (*chunk)[chunkPos++])
			{					
			case CHAR_GAME_OVER:
				go = false;
				break;
			case CHAR_SEQ_START:
				pos = 0;
				start = true;
				break;
			case CHAR_SEQ_END:
				if (posVertex.GetSize() == vertexLength && Within(NormHash(seed, posVertex, negVertex), low, high))
				{
					if (trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
					{
						trueBifSet.insert(posVertex.GetBody());
					}
				}
				
				negVertex = posVertex = DnaString();
				record++;
				break;
			default:
				if (posVertex.GetSize() == vertexLength)
				{
					if (!start)
					{
						if (isCandidBit[record]->Get(pos) && Within(NormHash(seed, posVertex, negVertex), low, high))
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
						++pos;
					}
					else
					{
						++pos;
						start = false;
						negVertex = posVertex.RevComp();
						if (Within(NormHash(seed, posVertex, negVertex), low, high) && trueBifSet.count(posVertex.GetBody()) == 0 && trueBifSet.count(negVertex.GetBody()) == 0)
						{
							trueBifSet.insert(posVertex.GetBody());
						}

						posVertex.AppendBack(posExtend);
						negVertex.AppendFront(DnaString::Reverse(posExtend));
						posPrev = posVertex.PopFront();
						negExtend = negVertex.PopBack();
					}
				}
				else
				{
					posVertex.AppendBack(posExtend);
				}

			}				
		}

		falseCount = candidateBifSet.size();
	}

	void DistributeTasks(const std::vector<std::string> & fileName, size_t overlapSize, std::vector<TaskQueuePtr> & taskQueue, std::vector<size_t> & fastaRecordsSize)
	{
		fastaRecordsSize.clear();
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
					if (buf.size() >= overlapSize && (buf.size() == Task::TASK_SIZE || over))
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
										overlap.assign(buf.end() - overlapSize + 1, buf.end());
									}

									q->push(Task(record, prev, over, std::move(buf)));
									prev = start - overlapSize + 1;
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
			while (taskQueue[i]->write_available() == 0)
			{
				boost::this_thread::sleep_for(boost::chrono::nanoseconds(1000000));
			}
			
			taskQueue[i]->push(Task(0, Task::GAME_OVER, true, std::string()));
		}
	}

	void Propagate(std::vector<StringQueuePtr> & stringQueue, StringPtr ptr)
	{
		for (size_t i = 0; i < stringQueue.size(); i++)
		{
			while (stringQueue[i]->write_available() == 0);
			stringQueue[i]->push(ptr);
		}
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, size_t hashFunctions, size_t rounds, size_t threads) :
		vertexSize_(vertexLength)
	{
		std::cout << "Threads = " << threads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << filterSize << std::endl;

		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}

		std::vector<uint64_t> seed(hashFunctions);
		std::generate(seed.begin(), seed.end(), rand);
		size_t edgeLength = vertexLength + 1;
		uint64_t low = 0;
		StringPtr seqEnd(new std::string(1, CHAR_SEQ_END));
		StringPtr seqStart(new std::string(1, CHAR_SEQ_START));
		StringPtr gameOver(new std::string(1, CHAR_GAME_OVER));
		for (size_t round = 0; round < rounds; round++)
		{			
			std::vector<size_t> fastaRecordsSize;
			std::vector<boost::thread> workerThread(threads);
			std::vector<std::unique_ptr<ConcurrentBitVector> > isCandidBit;
			uint64_t high = round == rounds - 1 ? UINT64_MAX : (UINT64_MAX / rounds) * (round + 1);
			time_t mark = time(0);
			{
				std::vector<TaskQueuePtr> taskQueue;
				ConcurrentBitVector bitVector(filterSize);				
				std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
				std::cout << "Counting\tEnumeration\tAggregation" << std::endl;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
					workerThread[i] = boost::thread(CountingWorker, low, high, boost::cref(seed), boost::ref(bitVector), edgeLength, boost::ref(*taskQueue[i]));
				}

				DistributeTasks(fileName, edgeLength, taskQueue, fastaRecordsSize);				
				for (size_t i = 0; i < workerThread.size(); i++)
				{				
					workerThread[i].join();
				}

				for (size_t sz : fastaRecordsSize)
				{
					isCandidBit.push_back(std::unique_ptr<ConcurrentBitVector>(new ConcurrentBitVector(sz)));
				}

				std::cout << time(0) - mark << "\t";
				mark = time(0);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateCheckingWorker, low, high, boost::cref(seed), boost::cref(bitVector), vertexLength, boost::ref(*taskQueue[i]), boost::ref(isCandidBit));
				}

				DistributeTasks(fileName, vertexLength, taskQueue, fastaRecordsSize);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i].join();
				}
			}

			std::cout << time(0) - mark << "\t";
			mark = time(0);
			std::vector<StringQueuePtr> charQueue;
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				charQueue.push_back(StringQueuePtr(new StringQueue(QUEUE_CAPACITY)));
			}

			uint64_t threadLow = low;
			std::vector<size_t> falseCount(workerThread.size(), 0);
			std::vector<std::unique_ptr<VertexSet> > trueBifSet(workerThread.size());
			for (size_t thread = 0; thread < workerThread.size(); thread++)
			{
				trueBifSet[thread] = std::unique_ptr<VertexSet>(new VertexSet(1024, VertexHashFunction(vertexLength), VertexEquality(vertexLength)));
				uint64_t threadHigh = thread == workerThread.size() - 1 ? high : ((high - low) / workerThread.size()) * (thread + 1);				
				workerThread[thread] = boost::thread(AggregationWorker, vertexLength, boost::ref(*charQueue[thread]), threadLow, threadHigh, seed, boost::ref(isCandidBit), boost::ref(*trueBifSet[thread]), boost::ref(falseCount[thread]));
				threadLow = threadHigh + 1;
			}

			StringPtr chunk;
			for (const std::string & nowFileName : fileName)
			{	
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); )
				{										
					char ch;
					Propagate(charQueue, seqStart);
					for (bool go = true; go; )
					{
						if(parser.GetChar(ch))
						{
							if (chunk == 0)
							{
								chunk.reset(new std::string());
							}

							chunk->push_back(ch);
						}
						else
						{
							go = false;
						}

						assert(chunk != 0);
						if (chunk->size() == Task::TASK_SIZE || !go)
						{
							Propagate(charQueue, chunk);
							chunk.reset();
						}
					}
					
					Propagate(charQueue, seqEnd);
				}
			}

			size_t trueCount = 0;
			Propagate(charQueue, gameOver);			
			for (size_t i = 0; i < workerThread.size(); i++)
			{				
				workerThread[i].join();
				trueCount += trueBifSet[i]->size();
				bifurcation_.insert(bifurcation_.end(), trueBifSet[i]->begin(), trueBifSet[i]->end());
				trueBifSet[i]->clear();
			}

			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << trueCount << std::endl;
			std::cout << "FP count = " << std::accumulate(falseCount.begin(), falseCount.end(), size_t(0)) << std::endl;
			std::cout << std::string(80, '-') << std::endl;
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