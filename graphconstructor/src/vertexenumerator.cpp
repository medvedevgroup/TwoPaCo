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
#include "ngramhashing/cyclichash.h"

#include "vertexenumerator.h"

namespace Sibelia
{
	const size_t VertexEnumerator::INVALID_VERTEX = -1;

	namespace
	{
		uint64_t hashMask;
		typedef CyclicHash<uint64_t> HashFunction;
		typedef std::unique_ptr<HashFunction> HashFunctionPtr;

		bool IsInBloomFilter(const ConcurrentBitVector & filter, const std::vector<uint64_t> & seed, const char * edge, uint64_t hash0, size_t edgeLength)
		{
			if (!filter.Get(hash0))
			{
				return false;
			}

			for (size_t i = 0; i < seed.size(); i++)
			{
				uint64_t hvalue = SpookyHash::Hash64(edge, edgeLength, seed[i]) & hashMask;
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
			Task(uint64_t start, bool isFinal, std::string && str) : start(start), isFinal(isFinal), str(std::move(str)) {}
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

		uint64_t NormHash(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash)
		{
			return std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
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

		std::string RevComp(const std::string & str)
		{
			std::string ret;
			for (auto it = str.rbegin(); it != str.rend(); ++it)
			{
				ret.push_back(DnaString::Reverse(*it));
			}

			return ret;
		}

		void InitializeHashFunctions(const std::vector<HashFunctionPtr> & seed,
			std::vector<HashFunctionPtr> & posEdgeHash,
			std::vector<HashFunctionPtr> & negEdgeHash,
			const std::string & fragment,
			size_t length,
			size_t offset = 0)
		{
			for (size_t i = 0; i < posEdgeHash.size(); i++)
			{
				posEdgeHash[i] = HashFunctionPtr(new HashFunction(*seed[i]));
				negEdgeHash[i] = HashFunctionPtr(new HashFunction(*seed[i]));
				for (auto it = fragment.begin() + offset; it != fragment.begin() + length + offset; ++it)
				{
					posEdgeHash[i]->eat(*it);
				}

				assert(posEdgeHash[i]->hashvalue == posEdgeHash[i]->hash(fragment.substr(offset, length)));
				for (std::string::const_reverse_iterator it(fragment.begin() + length + offset); it != fragment.rend() - offset; ++it)
				{
					char c = DnaString::Reverse(*it);
					negEdgeHash[i]->eat(DnaString::Reverse(*it));
				}

				assert(negEdgeHash[i]->hashvalue == negEdgeHash[i]->hash(RevComp(fragment.substr(offset, length))));
			}
		}

		size_t totalRecords;

		void CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
			const std::vector<HashFunctionPtr> & hashFunction,
			const std::vector<uint64_t> & seed,
			const ConcurrentBitVector & bitVector,
			size_t vertexLength,
			TaskQueue & taskQueue,
			std::ofstream & outFile,
			boost::mutex & outMutex,			
			std::unique_ptr<StreamFastaParser::Exception> & error)
		{
			uint64_t low = bound.first;
			uint64_t high = bound.second;
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

					std::string revTask = DnaString::RevComp(task.str);
					if (task.start == 0)
					{
						DnaString posVertex = Generate(task.str.begin(), vertexLength);
						DnaString negVertex = posVertex.RevComp();
						uint64_t posHash0 = hashFunction[0]->hash(posVertex.ToString());
						uint64_t negHash0 = hashFunction[0]->hash(negVertex.ToString());
						if (Within(std::min(negHash0, posHash0), low, high))
						{
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'A', 'A').GetBody());
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'C', 'A').GetBody());
						}
					}

					if (task.isFinal)
					{
						DnaString posVertex = Generate(task.str.end() - vertexLength, vertexLength);
						DnaString negVertex = posVertex.RevComp();
						uint64_t posHash0 = hashFunction[0]->hash(posVertex.ToString());
						uint64_t negHash0 = hashFunction[0]->hash(negVertex.ToString());
						if (Within(std::min(negHash0, posHash0), low, high))
						{
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'A', 'A').GetBody());
							output.push_back(MakeCanonicalRecord(posVertex, negVertex, 'C', 'A').GetBody());
						}
					}

					if (task.str.size() >= vertexLength + 2)
					{
						std::vector<HashFunctionPtr> posVertexHash(1);
						std::vector<HashFunctionPtr> negVertexHash(1);
						InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
						for (size_t pos = 1;; ++pos)
						{
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							if (Within(std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue), low, high))
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
										uint64_t posHash0 = posVertexHash[0]->hash_prepend(nextCh);
										uint64_t negHash0 = negVertexHash[0]->hash_extend(revNextCh);
										assert(posHash0 == posVertexHash[0]->hash(std::string(1, nextCh) + task.str.substr(pos, vertexLength)));
										assert(negHash0 == negVertexHash[0]->hash(RevComp(std::string(1, nextCh) + task.str.substr(pos, vertexLength))));
										if (posHash0 <= negHash0)
										{
											if (IsInBloomFilter(bitVector, seed, &body[pos], posHash0))
											{
												outCount++;
											}
										}
										else
										{
											char revNextCh = DnaString::Reverse(nextCh);
											if (IsInBloomFilter(bitVector, seed,  , negHash0))
											{
												outCount++;
											}
										}
									}

									if (nextCh == posExtend)
									{
										++outCount;
									}
									else
									{
										uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
										uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
										assert(posHash0 == posVertexHash[0]->hash(task.str.substr(pos, vertexLength) + nextCh));
										assert(negHash0 == negVertexHash[0]->hash(RevComp(task.str.substr(pos, vertexLength) + nextCh)));
										if (posHash0 <= negHash0)
										{
											if (IsInBloomFilter(bitVector, seed, posOutEdge, posHash0))
											{
												outCount++;
											}
										}
										else
										{
											if (IsInBloomFilter(bitVector, seed, negOutEdge, negHash0))
											{
												outCount++;
											}
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
								char negExtend = DnaString::Reverse(posExtend);
								char negPrev = DnaString::Reverse(task.str[pos + vertexLength]);
								for (size_t i = 0; i < posVertexHash.size(); i++)
								{
									posVertexHash[i]->update(posPrev, posExtend);
									negVertexHash[i]->reverse_update(negExtend, negPrev);
									assert(posVertexHash[i]->hashvalue == posVertexHash[0]->hash(task.str.substr(pos, vertexLength)));
									assert(negVertexHash[i]->hashvalue == negVertexHash[0]->hash(RevComp(task.str.substr(pos, vertexLength))));
								}
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
						totalRecords += output.size();
						outFile.write(reinterpret_cast<const char*>(&output[0]), output.size() * sizeof(output[0]));
						if (error != 0)
						{
							return;
						}

						if (!outFile)
						{
							error.reset(new StreamFastaParser::Exception("Can't write to the temporary file"));
							return;
						}
					}
				}
			}
		}

		const uint64_t BINS_COUNT = 1 << 24;
		const uint32_t MAX_COUNTER = UINT32_MAX >> 1;

		void FilterFillerWorker(uint64_t low,
			uint64_t high,
			const std::vector<HashFunctionPtr> & hashFunction,
			const std::vector<uint64_t> & seed,
			ConcurrentBitVector & filter,
			size_t edgeLength,
			TaskQueue & taskQueue)
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
										
					size_t vertexLength = edgeLength - 1;
					size_t revPos = task.str.size() - edgeLength;
					std::vector<HashFunctionPtr> posVertexHash(1);
					std::vector<HashFunctionPtr> negVertexHash(1);
					std::string revTask = DnaString::RevComp(task.str);
					InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength);
					for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
					{
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaString::Reverse(nextCh);
						uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
						uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
						uint64_t fistMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);						
						char prevCh = task.str[pos];
						for (size_t i = 0; i < posVertexHash.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(DnaString::Reverse(nextCh), DnaString::Reverse(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(RevComp(task.str.substr(pos + 1, vertexLength))));
						}

						uint64_t secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						const char * body = posHash0 <= negHash0 ? &task.str[pos] : &revTask[revPos];
						if (Within(fistMinHash0, low, high) || Within(secondMinHash0, low, high))
						{
							filter.SetConcurrently(std::min(posHash0, negHash0));
							for (uint64_t s : seed)
							{
								uint64_t hvalue = SpookyHash::Hash64(body, edgeLength, s) & hashMask;
								filter.SetConcurrently(hvalue);
							}
						}
					}
				}
			}
		}
		/*
		boost::mutex counterMutex;
		std::unordered_set<uint64_t> edgeCounter;
		std::vector<uint32_t> trueBin(BINS_COUNT, 0);*/


		void InitialFilterFillerWorker(uint64_t binSize,
			const std::vector<HashFunctionPtr> & hashFunction,
			ConcurrentBitVector & filter,
			size_t vertexLength,
			TaskQueue & taskQueue,
			std::atomic<uint32_t> * binCounter)
		{
			size_t edgeLength = vertexLength + 1;
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

					size_t vertexLength = edgeLength - 1;
					std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
					std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
					InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength);
					for (size_t pos = 0; pos + edgeLength - 1 < task.str.size(); ++pos)
					{
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaString::Reverse(nextCh);
						uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
						uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);
						uint64_t fistMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						char prevCh = task.str[pos];
						bool wasSet = true;
						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							uint64_t hvalue = posHash0 <= negHash0 ? posVertexHash[i]->hash_extend(nextCh) : negVertexHash[i]->hash_prepend(revNextCh);
							if (!filter.Get(hvalue))
							{
								wasSet = false;
								filter.SetConcurrently(hvalue);
							}
						}
											
						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(DnaString::Reverse(nextCh), DnaString::Reverse(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(RevComp(task.str.substr(pos + 1, vertexLength))));
						}


						uint64_t secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						if (!wasSet)
						{
							uint64_t value[] = { fistMinHash0, secondMinHash0 };
							for (uint64_t v : value)
							{
								uint64_t bin = v / binSize;
								if (binCounter[bin] < MAX_COUNTER)
								{
									binCounter[bin].fetch_add(1);
								}
							}
						}
						/*
						{
							DnaString edge(task.str.substr(pos, edgeLength));
							if (!(posHash0 <= negHash0))
							{
								edge = edge.RevComp();
							}

							boost::lock_guard<boost::mutex> g(counterMutex);
							if (edgeCounter.count(edge.GetBody()) == 0)
							{
								edgeCounter.insert(edge.GetBody());
								uint64_t value[] = { fistMinHash0, secondMinHash0 };
								for (uint64_t v : value)
								{
									uint64_t bin = v / binSize;
									if (trueBin[bin] < MAX_COUNTER)
									{
										trueBin[bin] += 1;
									}
								}
							}
						}*/
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
			std::vector<uint64_t> * out;
			boost::mutex * outMutex;
			TrueBifurcations(std::vector<uint64_t> * candidate, std::vector<uint64_t> * out, boost::mutex * outMutex, size_t vertexSize) :
				candidate(candidate), out(out), vertexSize(vertexSize), outMutex(outMutex), falsePositives(0) {}

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
				for (; i < range.end();)
				{					
					bool bif = false;
					Candidate icand(vertexSize, (*candidate)[i]);
					bool selfRevComp = icand.base == icand.base.RevComp();					
					size_t j = i;
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
							if (selfRevComp)
							{
								bif = bif || icand.prev != DnaString::Reverse(jcand.extend) || icand.extend != DnaString::Reverse(jcand.prev);
							}
						}
					}

					if (bif)
					{
						boost::lock_guard<boost::mutex> guard(*outMutex);
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
		uint64_t realSize = uint64_t(1) << filterSize;
		hashMask = realSize - 1;
		std::cout << "Threads = " << threads << std::endl;
		std::cout << "Aggregation threads = " << aggregationThreads << std::endl;
		std::cout << "Hash functions = " << hashFunctions << std::endl;
		std::cout << "Filter size = " << realSize << std::endl;
		std::cout << "Files: " << std::endl;
		for (const std::string & fn : fileName)
		{
			std::cout << fn << std::endl;
		}

		std::vector<uint64_t> seed(hashFunctions - 1);
		std::generate(seed.begin(), seed.end(), rand);
		std::cout << std::string(80, '-') << std::endl;
		const uint64_t BIN_SIZE = std::max(uint64_t(1), realSize / BINS_COUNT);
		if (vertexLength > 30)
		{
			throw std::runtime_error("The vertex size is too large");
		}

		std::vector<HashFunctionPtr> hashFunction(hashFunctions);
		for (HashFunctionPtr & ptr : hashFunction)
		{
			ptr = HashFunctionPtr(new HashFunction(vertexLength, filterSize));
		}

		size_t edgeLength = vertexLength + 1;
		std::atomic<uint32_t> * binCounter = new std::atomic<uint32_t>[BINS_COUNT];
		{
			std::fill(binCounter, binCounter + BINS_COUNT, 0);
			std::vector<boost::thread> workerThread(threads);
			ConcurrentBitVector bitVector(realSize);
			std::vector<TaskQueuePtr> taskQueue;
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
				workerThread[i] = boost::thread(InitialFilterFillerWorker,					
					BIN_SIZE,
					boost::cref(hashFunction),
					boost::ref(bitVector),
					vertexLength,
					boost::ref(*taskQueue[i]),
					binCounter);
			}

			DistributeTasks(fileName, edgeLength, taskQueue);
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				workerThread[i].join();
			}
		}


		rounds = 1;
		double roundSize = 0;
		for (;; ++rounds)
		{
			roundSize = std::accumulate(binCounter, binCounter + BINS_COUNT, 0.f) / rounds;
			if (realSize / roundSize > 16)
			{
				break;
			}
		}

		/*
		long long diff = 0;
		size_t nonZero = 0;
		std::ofstream diffOut("diff.txt");
		for (size_t i = 0; i < BINS_COUNT; i++)
		{
			long long c1 = binCounter[i];
			long long c2 = trueBin[i];
			diff += std::abs(c1 - c2);
			if (c1 || c2)
			{
				diffOut << c1 << " " << c2 << std::endl;
				nonZero++;
			}
		}

		std::cout << "Diff = " << diff << std::endl;
		std::cout << "Avg diff = " << double(diff) / nonZero << std::endl;*/

		uint64_t low = 0;
		uint64_t high = 0;
		size_t lowBoundary = 0;
		uint64_t totalFpCount = 0;
		for (size_t round = 0; round < rounds; round++)
		{
			totalRecords = 0;
			time_t mark = time(0);			
			uint64_t accumulated = binCounter[lowBoundary];
			for (++lowBoundary; lowBoundary < BINS_COUNT; ++lowBoundary)
			{				
				if (accumulated <= roundSize || round + 1 == rounds)
				{
					accumulated += binCounter[lowBoundary];
				}
				else
				{					
					break;
				}
			}

			std::cout << "Ratio = " << double(realSize) / accumulated << std::endl;
			uint64_t high = lowBoundary * BIN_SIZE;
			{
				std::vector<TaskQueuePtr> taskQueue;
				ConcurrentBitVector bitVector(realSize);
				std::vector<boost::thread> workerThread(threads);
				std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
				std::cout << "Counting\tEnumeration\tAggregation" << std::endl;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue.push_back(TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY)));
					workerThread[i] = boost::thread(FilterFillerWorker,
						low,
						high,
						boost::cref(hashFunction),
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

				std::unique_ptr<StreamFastaParser::Exception> error;
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateCheckingWorker,
						std::make_pair(low, high),						
						boost::cref(hashFunction),
						boost::cref(seed),
						boost::cref(bitVector),
						vertexLength,
						boost::ref(*taskQueue[i]),
						boost::ref(tmpFile),
						boost::ref(tmpFileMutex),
						boost::ref(error));
				}

				if (error != 0)
				{
					throw StreamFastaParser::Exception(*error);
				}

				DistributeTasks(fileName, vertexLength + 1, taskQueue);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i].join();
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
				if (!tmpFile)
				{
					throw StreamFastaParser::Exception("The temporary file is corrupted");
				}
			}

			boost::mutex outMutex;
			tbb::parallel_sort(candidate.begin(), candidate.end(), VertexLess(vertexSize_));
			uint64_t falsePositives = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, candidate.size()),
				uint64_t(0),
				TrueBifurcations(&candidate, &bifurcation_, &outMutex, vertexSize_),
				std::plus<uint64_t>());

			std::cout << time(0) - mark << std::endl;
			std::cout << "Vertex count = " << bifurcation_.size() << std::endl;
			std::cout << "FP count = " << falsePositives << std::endl;
			std::cout << "Records = " << candidate.size() << std::endl;
			std::cout << std::string(80, '-') << std::endl;
			totalFpCount += falsePositives;
			low = high + 1;
		}

		std::cout << "Total FPs = " << totalFpCount << std::endl;
		delete[] binCounter;
		tbb::parallel_sort(bifurcation_.begin(), bifurcation_.end());
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