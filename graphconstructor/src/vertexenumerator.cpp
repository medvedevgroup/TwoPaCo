#include <set>
#include <deque>
#include <ctime>
#include <memory>
#include <bitset>
#include <numeric>
#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>

#include <boost/ref.hpp>
#include <boost/lockfree/spsc_queue.hpp>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>

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
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]) % filter.Size();
				filter.SetConcurrently(hvalue);
			}
		}

		bool IsInBloomFilter(const ConcurrentBitVector & filter, const std::vector<uint64_t> & seed, const DnaString & item)
		{
			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t body = item.GetBody();
				uint64_t hvalue = SpookyHash::Hash64(&body, sizeof(body), seed[i]) % filter.Size();
				if (!filter.Get(hvalue))
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
				DnaString stra(vertexSize_, a);
				DnaString strb(vertexSize_, b);
				return stra.GetBody() < strb.GetBody();
			}


		private:			
			size_t vertexSize_;
		};

		size_t CharIndex(char ch)
		{
			return std::find(DnaString::LITERAL.begin(), DnaString::LITERAL.end(), ch) - DnaString::LITERAL.begin();
		}

		struct Task
		{
			bool isFinal;
			uint64_t start;
			std::string str;
#ifdef _DEBUG
			static const size_t TASK_SIZE = 36;
#else
			static const size_t TASK_SIZE = 1 << 18;
#endif			
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(uint64_t start, bool isFinal, std::string && str) :  start(start), isFinal(isFinal), str(std::move(str)) {}
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

		DnaString MakeCanonicalRecord(DnaString posVertex, DnaString negVertex, char posExtend, char posPrev)
		{
			if (posVertex.GetBody() < negVertex.GetBody())
			{
				posVertex.AppendBack(posExtend);
				posVertex.AppendBack(posPrev);
				return posVertex;
			}

			negVertex.AppendBack(DnaString::Reverse(posPrev));
			negVertex.AppendBack(DnaString::Reverse(posExtend));
			return negVertex;
		}

		DnaString Generate(std::string::const_iterator it, size_t size)
		{
			DnaString str;
			for (size_t i = 0; i < size; i++)
			{
				str.AppendBack(*(it + i));
			}

			return str;
		}

		void CandidateCheckingWorker(uint64_t low,
			uint64_t high,
			const std::vector<uint64_t> & seed,
			const ConcurrentBitVector & bitVector,
			size_t vertexLength,
			TaskQueue & taskQueue,
			std::ofstream & outFile,
			boost::mutex & outMutex,
			size_t & total)
		{
			std::vector<size_t> hf(seed.size());
			std::vector<uint64_t> output;
			while (true)
			{
				Task task;
				if (taskQueue.pop(task))
				{
 					output.clear();
					if (task.start == Task::GAME_OVER)
					{
						break;
					}

					if (task.str.size() < vertexLength)
					{
						continue;
					}

					if (task.start == 0)
					{
						DnaString posVertex = Generate(task.str.begin(), vertexLength);
						DnaString negVertex = posVertex.RevComp();
						if (Within(NormHash(seed, posVertex, negVertex), low, high))
						{
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'A', 'A').GetBody());
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'C', 'A').GetBody());
						}
					}

					if (task.isFinal)
					{
						DnaString posVertex = Generate(task.str.end() - vertexLength, vertexLength);
						DnaString negVertex = posVertex.RevComp();
						if (Within(NormHash(seed, posVertex, negVertex), low, high))
						{
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'A', 'A').GetBody());
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'C', 'A').GetBody());
						}
					}

					if (task.str.size() >= vertexLength + 2)
					{
						DnaString posVertex = Generate(task.str.begin() + 1, vertexLength);	
						DnaString negVertex = posVertex.RevComp();
						for (size_t pos = 1;; ++pos)
						{
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							assert(posVertex.RevComp() == negVertex);
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
									if (nextCh == posPrev)
									{
										++inCount;
									}
									else
									{
										DnaString posInEdge = posVertex;
										DnaString negInEdge = negVertex;
										posInEdge.AppendFront(nextCh);
										negInEdge.AppendBack(DnaString::Reverse(nextCh));
										if ((posInEdge.GetBody() <= negInEdge.GetBody() && IsInBloomFilter(bitVector, seed, posInEdge)) ||
											(negInEdge.GetBody() < posInEdge.GetBody() && IsInBloomFilter(bitVector, seed, negInEdge)))
										{
											inCount++;
										}
									}

									if (nextCh == posExtend)
									{
										++outCount;
									}
									else
									{
										DnaString posOutEdge = posVertex;
										DnaString negOutEdge = negVertex;
										posOutEdge.AppendBack(nextCh);
										negOutEdge.AppendFront(DnaString::Reverse(nextCh));
										if ((posOutEdge.GetBody() <= negOutEdge.GetBody() && IsInBloomFilter(bitVector, seed, posOutEdge)) ||
											(negOutEdge.GetBody() < posOutEdge.GetBody() && IsInBloomFilter(bitVector, seed, negOutEdge)))
										{
											outCount++;
										}
									}
								}

								if (inCount > 1 || outCount > 1)
								{
									output.push_back(MakeCanonicalRecord(posVertex, negVertex, posExtend, posPrev).GetBody());
								}
							}

							if (pos + vertexLength + 1 < task.str.size())
							{
								posVertex.AppendBack(posExtend);								
								negVertex.AppendFront(DnaString::Reverse(posExtend));
								posVertex.PopFront();
								negVertex.PopBack();
							}
							else
							{
								break;
							}
						}
						
					}

					if (output.size() > 0)
					{
						boost::lock_guard<boost::mutex> guard(outMutex);
						total += output.size();						
						outFile.write(reinterpret_cast<const char*>(&output[0]), output.size() * sizeof(output[0]));
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
							if (posEdge.GetBody() < negEdge.GetBody())
							{
								PutInBloomFilter(bitVector, seed, posEdge);
							}
							else
							{
								PutInBloomFilter(bitVector, seed, negEdge);
							}
							
						}

						posEdge.PopFront();
						negEdge.PopBack();
					}
				}
			}
		}

		void DistributeTasks(const std::vector<std::string> & fileName, size_t overlapSize, std::vector<TaskQueuePtr> & taskQueue)
		{
			for (size_t file = 0; file < fileName.size(); file++)
			{
				size_t record = 0;
				const std::string & nowFileName = fileName[file];
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
											overlap.assign(buf.end() - overlapSize, buf.end());
										}

										q->push(Task(prev, over, std::move(buf)));
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

				taskQueue[i]->push(Task(Task::GAME_OVER, true, std::string()));
			}
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

		struct TrueBifurcations
		{
			size_t vertexSize;
			uint64_t falsePositives;
			std::vector<uint64_t> * candidate;
			tbb::concurrent_vector<uint64_t> * out;
			TrueBifurcations(std::vector<uint64_t> * candidate, tbb::concurrent_vector<uint64_t> * out, size_t vertexSize) :
				candidate(candidate), out(out), vertexSize(vertexSize), falsePositives(0) {}

			uint64_t operator()(const tbb::blocked_range<size_t> & range, uint64_t init) const
			{
				size_t i = range.begin();
				if (i > 0)
				{
					DnaString istr(vertexSize, (*candidate)[i]);
					if (istr == DnaString(vertexSize, (*candidate)[i - 1]))
					{
						for (++i; i < range.end() && istr == DnaString(vertexSize, (*candidate)[i]); ++i);
					}					
				}

				uint64_t falsePositives = init;
				for (; i < range.end(); )
				{
					size_t j = i + 1;
					bool bif = false;
					Candidate icand(vertexSize, (*candidate)[i]);
					if (icand.base == icand.base.RevComp())
					{
						Candidate rcand = icand.Reverse();
						bif = icand.prev != rcand.prev || icand.extend != rcand.extend;
					}

					for (; j < candidate->size(); j++)
					{
						Candidate jcand(vertexSize, (*candidate)[j]);
						if (jcand.base != icand.base)
						{
							break;
						}
						else if (!bif)
						{
							bif = icand.prev != jcand.prev || icand.extend != jcand.extend;
						}
					}

					if (bif)
					{
						out->push_back(icand.base.GetBody());
					}
					else
					{
						falsePositives++;
					}

					i = j;
				}

				return falsePositives;
			}
		};

	}

	VertexEnumerator::VertexEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t aggregationThreads,
		const std::string & tmpFileName) :
		vertexSize_(vertexLength)
	{		
		std::cout << "Threads = " << threads << std::endl;
		std::cout << "Aggregation threads = " << aggregationThreads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << filterSize << std::endl;
		std::cout << "Files: " << std::endl;
		for (const std::string & fn : fileName)
		{
			std::cout << fn << std::endl;
		}

		std::cout << std::string(80, '-') << std::endl;

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
			time_t mark = time(0);
			size_t totalRecords = 0;
			uint64_t high = round == rounds - 1 ? UINT64_MAX : (UINT64_MAX / rounds) * (round + 1);
			{
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
						workerThread[i] = boost::thread(CountingWorker,
							low,
							high,
							boost::cref(seed),
							boost::ref(bitVector),
							edgeLength,
							boost::ref(*taskQueue[i]));
					}

					DistributeTasks(fileName, edgeLength, taskQueue);				
					for (size_t i = 0; i < workerThread.size(); i++)
					{				
						workerThread[i].join();
					}
					
					std::cout << time(0) - mark << "\t";
					mark = time(0);					
					boost::mutex tmpFileMutex;
					std::ofstream tmpFile(tmpFileName.c_str(), std::ios_base::binary);
					if (!tmpFile)
					{
						throw StreamFastaParser::Exception("Can't open the temporary file");
					}

					for (size_t i = 0; i < workerThread.size(); i++)
					{
						workerThread[i] = boost::thread(CandidateCheckingWorker,
							low,
							high,
							boost::cref(seed),
							boost::cref(bitVector),
							vertexLength,
							boost::ref(*taskQueue[i]),
							boost::ref(tmpFile),
							boost::ref(tmpFileMutex),
							boost::ref(totalRecords));
					}

					DistributeTasks(fileName, vertexLength + 1, taskQueue);
					for (size_t i = 0; i < taskQueue.size(); i++)
					{
						workerThread[i].join();
					}
				}
			
				std::cout << time(0) - mark << "\t";				
			}

			mark = time(0);
			std::vector<uint64_t> candidate(totalRecords);
			std::ifstream tmpFile(tmpFileName.c_str(), std::ios_base::binary);
			if (!tmpFile)
			{
				throw StreamFastaParser::Exception("Can't open the temporary file");
			}

			if (totalRecords > 0)
			{ 
				tmpFile.read(reinterpret_cast<char*>(&candidate[0]), totalRecords * sizeof(candidate[0]));
			}

			tbb::parallel_sort(candidate.begin(), candidate.end(), VertexLess(vertexSize_));
			uint64_t falsePositives = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, candidate.size()),
				uint64_t(0),
				TrueBifurcations(&candidate, &bifurcation_, vertexSize_),
				std::plus<uint64_t>());

			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << bifurcation_.size() << std::endl;
			std::cout << "FP count = " << falsePositives << std::endl;
			std::cout << "Records = " << candidate.size() << std::endl;
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
			tbb::concurrent_vector<uint64_t>::const_iterator it = std::lower_bound(bifurcation_.begin(), bifurcation_.end(), str.GetBody());
			if (it != bifurcation_.end() && *it == str.GetBody())
			{
				return it - bifurcation_.begin();
			}
		}

		return INVALID_VERTEX;
	}
}