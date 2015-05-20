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
#include <boost/lockfree/spsc_queue.hpp>
#include "lib/SpookyV2.h"

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

		class VertexLess
		{
		public:
			VertexLess(size_t vertexSize) : vertexSize_(vertexSize)
			{

			}

			bool operator () (const uint64_t & a, const uint64_t & b) const
			{
				return CanonicalKmer(a) < CanonicalKmer(b);
			}

			bool Equal(const uint64_t & a, const uint64_t & b)
			{
				return CanonicalKmer(a) == CanonicalKmer(b);
			}

		private:			
			uint64_t CanonicalKmer(uint64_t kmer) const
			{
				DnaString str(vertexSize_, kmer);
				return std::min(str.GetBody(), str.RevComp().GetBody());
			}

			size_t vertexSize_;
		};

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

		DnaString kmer;			

		const size_t QUEUE_CAPACITY = 16;
		typedef boost::lockfree::spsc_queue<Task> TaskQueue;
		typedef std::unique_ptr<TaskQueue> TaskQueuePtr;

		std::string TempFile(size_t record)
		{
			std::stringstream ss;
			ss << record << ".bin";
			return ss.str();
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

		void DistributeTasks(const std::vector<std::string> & fileName, size_t overlapSize, std::vector<TaskQueuePtr> & taskQueue, std::vector<size_t> & fastaRecordsSize)
		{
			fastaRecordsSize.clear();
			for (size_t file = 0; file < fileName.size(); file++)
			{
				size_t record = 0;
				const std::string & nowFileName = fileName[file];
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

		DnaString MakeRecord(DnaString kmer, char extend, char prev)
		{
			kmer.AppendBack(extend);
			kmer.AppendBack(prev);
			return kmer;
		}

		void ParseRecord(size_t vertexSize, uint64_t kmer, char & extend, char & prev)
		{
			DnaString str(vertexSize, kmer);
			extend = str.GetChar(vertexSize);
			prev = str.GetChar(vertexSize + 1);
		}

		struct Candidate
		{
			char prev;
			char extend;
			DnaString base;

			Candidate(size_t vertexSize, uint64_t record) : base(vertexSize + 2, record)
			{
				extend = base.GetChar(vertexSize);
				prev = base.GetChar(vertexSize + 1);
				base.PopBack();
				base.PopBack();
			}

			Candidate Reverse() const
			{
				Candidate ret(*this);
				ret.base = ret.base.RevComp();
				ret.extend = DnaString::Reverse(prev);
				ret.prev = DnaString::Reverse(extend);
				return ret;
			}
		};

		void DetectTrueBifurcations(std::vector<uint64_t> & candidate, std::vector<uint64_t> & bifurcations, size_t vertexSize)
		{			
			for (size_t i = 0; i < candidate.size(); )
			{
				size_t j = i + 1;
				bool bif = false;				
				Candidate icand(vertexSize, candidate[i]);
				if (icand.base == icand.base.RevComp())
				{
					Candidate rcand = icand.Reverse();
					bif = icand.prev != rcand.prev || icand.extend != rcand.extend;
				}

				for (; j < candidate.size(); j++)
				{
					Candidate jcand(vertexSize, candidate[j]);
					if (jcand.base != icand.base && jcand.base != icand.base.RevComp())
					{
						break;
					}
					else if (!bif)
					{
						if (jcand.base != icand.base)
						{
							jcand = jcand.Reverse();
						}

						bif = icand.prev != jcand.prev || icand.extend != jcand.extend;
					}
				}
				
				if (bif)
				{
					bifurcations.push_back(icand.base.GetBody());
				}

				i = j;
			}
		}
	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName, size_t vertexLength, size_t filterSize, size_t hashFunctions, size_t rounds, size_t threads, size_t aggregationThreads) :
		vertexSize_(vertexLength)
	{		
		std::cout << "Threads = " << threads << std::endl;
		std::cout << "Aggregation threads = " << aggregationThreads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << filterSize << std::endl;

		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}

		std::vector<DnaString> border;
		std::vector<uint64_t> seed(hashFunctions);
		std::generate(seed.begin(), seed.end(), rand);
		size_t edgeLength = vertexLength + 1;
		uint64_t low = 0;
		for (size_t round = 0; round < rounds; round++)
		{
			time_t mark = time(0);
			std::vector<uint64_t> candidate;
			uint64_t high = round == rounds - 1 ? UINT64_MAX : (UINT64_MAX / rounds) * (round + 1);
			{
				std::vector<size_t> fastaRecordsSize;
				std::vector<std::unique_ptr<ConcurrentBitVector> > isCandidBit;
				{
					std::vector<TaskQueuePtr> taskQueue;
					ConcurrentBitVector bitVector(filterSize);				
					std::vector<boost::thread> workerThread(threads);
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
							size_t kmer = 0;						
							DnaString negVertex = posVertex.RevComp();
							if (Within(NormHash(seed, posVertex, negVertex), low, high))
							{
								candidate.push_back(MakeRecord(posVertex, 'A', 'A').GetBody());
								candidate.push_back(MakeRecord(posVertex, 'C', 'A').GetBody());
							}

							for (bool go = true; go; kmer++)
							{
								if (go = parser.GetChar(posExtend))
								{
									if (kmer > 0 && isCandidBit[record]->Get(kmer))
									{
										candidate.push_back(MakeRecord(posVertex, posExtend, posPrev).GetBody());
									}

									posVertex.AppendBack(posExtend);
									negVertex.AppendFront(DnaString::Reverse(posExtend));
									posPrev = posVertex.PopFront();
									negExtend = negVertex.PopBack();
								}
							}

							if (Within(NormHash(seed, posVertex, negVertex), low, high))
							{
								candidate.push_back(MakeRecord(posVertex, 'A', 'A').GetBody());
								candidate.push_back(MakeRecord(posVertex, 'C', 'A').GetBody());
							}
						}
					}
				}
			}

			std::sort(candidate.begin(), candidate.end(), VertexLess(vertexSize_));
			DetectTrueBifurcations(candidate, bifurcation_, vertexSize_);
			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << bifurcation_.size() << std::endl;
			std::cout << "FP count = " << candidate.size() << std::endl;
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