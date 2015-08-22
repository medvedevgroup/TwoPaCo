#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#define MAX_CAPACITY 5

#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_set>

#include <tpie/tpie.h>
#include <tpie/sort.h>
#include <tpie/dummy_progress.h>
#include <tpie/file_stream.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/concurrent_vector.h>

#include <boost/ref.hpp>
#include <boost/filesystem.hpp> 
#include <boost/lockfree/spsc_queue.hpp>

#include "ngramhashing/cyclichash.h"

#include "streamfastaparser.h"
#include "candidateoccurence.h"
#include "concurrentbitvector.h"

namespace Sibelia
{
	class VertexEnumerator
	{
	public:
		const static size_t INVALID_VERTEX = -1;
		virtual size_t GetVerticesCount() const = 0;
		virtual size_t GetId(const std::string & vertex) const = 0;
		virtual void Dump(std::vector<std::string> & toCopy) const = 0;

		virtual ~VertexEnumerator()
		{

		}
	};

	std::unique_ptr<VertexEnumerator> CreateEnumerator(const std::vector<std::string> & fileName,
		size_t vertexLength,
		size_t filterSize,
		size_t hashFunctions,
		size_t rounds,
		size_t threads,
		size_t aggregationThreads,
		const std::string & tmpFileName,
		const std::string & outFileName);

	template<size_t CAPACITY>
	class VertexEnumeratorImpl : public VertexEnumerator
	{
	private:
		class VertexLess;
		static const size_t BUF_SIZE = 1 << 24;
	public:
		typedef CompressedString<CAPACITY> DnaString;
		typedef CandidateOccurence<CAPACITY> Occurence;

		size_t GetVerticesCount() const
		{
			return bifurcation_.size();
		}

		size_t GetId(const std::string & vertex) const
		{
			DnaString str;
			str.CopyFromString(vertex.begin(), vertexSize_);
			typename std::vector<DnaString>::const_iterator it = std::lower_bound(bifurcation_.begin(), bifurcation_.end(), str, DnaStringLess(vertexSize_));
			if (it != bifurcation_.end() && *it == str)
			{
				return it - bifurcation_.begin();
			}

			return INVALID_VERTEX;
		}

		virtual void Dump(std::vector<std::string> & toCopy) const
		{
			std::set<std::string> ret;
			for (const DnaString & str : bifurcation_)
			{
				std::string pos = str.ToString(vertexSize_);
				std::string neg = DnaChar::ReverseCompliment(pos);
				assert(ret.count(pos) == 0 && ret.count(neg) == 0);
				ret.insert(pos);
				ret.insert(neg);
			}

			toCopy.assign(ret.begin(), ret.end());
		}

		VertexEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t aggregationThreads,
			const std::string & tmpFileName,
			const std::string & outFileName) :
			vertexSize_(vertexLength)
		{
			uint64_t realSize = uint64_t(1) << filterSize;
			std::cout << "Threads = " << threads << std::endl;
			std::cout << "Vertex length = " << vertexLength << std::endl;
			std::cout << "Aggregation threads = " << aggregationThreads << std::endl;
			std::cout << "Hash functions = " << hashFunctions << std::endl;
			std::cout << "Filter size = " << realSize << std::endl;
			std::cout << "Capacity = " << CAPACITY << std::endl;
			std::cout << "Files: " << std::endl;
			for (const std::string & fn : fileName)
			{
				std::cout << fn << std::endl;
			}

			const uint64_t BIN_SIZE = std::max(uint64_t(1), realSize / BINS_COUNT);
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
				roundSize = double(std::accumulate(binCounter, binCounter + BINS_COUNT, size_t(0))) / rounds;
				if (realSize / roundSize >= 8)
				{
					break;
				}
			}

			tpie::tpie_init(tpie::ALL);
			tpie::progress_indicator_null pi;
			tpie::get_memory_manager().set_limit(std::max(uint64_t(1024 * 1024 * 1024), realSize / 8));
			tpie::file_stream<VertexRecord> outFile;
			outFile.open(outFileName.c_str(), tpie::access_write);
			std::cout << "Round size = " << realSize / roundSize << std::endl;
			std::cout << std::string(80, '-') << std::endl;
			uint64_t low = 0;
			uint64_t high = 0;	
			size_t lowBoundary = 0;
			uint64_t totalFpCount = 0;
			uint64_t verticesCount = 0;			
			for (size_t round = 0; round < rounds; round++)
			{
				size_t totalRecords = 0;
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
					tpie::file_stream<Occurence> tmpFile;
					tmpFile.open(tmpFileName.c_str(), tpie::access_write);
					std::unique_ptr<StreamFastaParser::Exception> error;
					for (size_t i = 0; i < workerThread.size(); i++)
					{
						workerThread[i] = boost::thread(CandidateCheckingWorker,
							std::make_pair(low, high),
							boost::cref(hashFunction),
							boost::cref(bitVector),
							vertexLength,
							boost::ref(totalRecords),
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
				tpie::file_stream<Occurence> tmpFile;
				tmpFile.open(tmpFileName.c_str());
				tpie::sort(tmpFile);
				boost::mutex outMutex;
			 	tmpFile.seek(0);
				uint64_t falsePositives = TrueBifurcations(&tmpFile, &bifurcation_, &outMutex, vertexSize_)(totalRecords, 0, verticesCount, outFile);								
				std::cout << time(0) - mark << std::endl;
				std::cout << "Vertex count = " << bifurcation_.size() << std::endl;
				std::cout << "FP count = " << falsePositives << std::endl;
				std::cout << "Records = " << totalRecords << std::endl;
				std::cout << std::string(80, '-') << std::endl;
				totalFpCount += falsePositives;
				low = high + 1;
			}

			outFile.close();
			outFile.open(outFileName.c_str());
			tpie::sort(outFile, pi);
			delete[] binCounter;
			std::cout << "Total FPs = " << totalFpCount << std::endl;
			tbb::parallel_sort(bifurcation_.begin(), bifurcation_.end(), DnaStringLess(vertexSize_));
		}

	private:

		

		static const size_t QUEUE_CAPACITY = 16;
		static const uint64_t BINS_COUNT = 1 << 24;
		static const uint32_t MAX_COUNTER = UINT32_MAX >> 1;

		struct Task
		{
			bool isFinal;
			uint64_t start;
			uint64_t seqId;
			std::string str;
#ifdef _DEBUG
			static const size_t TASK_SIZE = 36;
#else
			static const size_t TASK_SIZE = 1 << 18;
#endif			
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(uint64_t seqId, uint64_t start, bool isFinal, std::string && str) : 
				seqId(seqId), start(start), isFinal(isFinal), str(std::move(str)) {}
		};		

		typedef CyclicHash<uint64_t> HashFunction;
		typedef std::unique_ptr<HashFunction> HashFunctionPtr;
		typedef boost::lockfree::spsc_queue<Task> TaskQueue;
		typedef std::unique_ptr<TaskQueue> TaskQueuePtr;

		enum StrandComparisonResult
		{
			positiveLess,
			negativeLess,
			tie
		};

		static StrandComparisonResult DetermineStrandExtend(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash, char nextCh, char revNextCh)
		{
			for (size_t i = 0; i < posVertexHash.size(); i++)
			{
				uint64_t posHash = posVertexHash[i]->hash_extend(nextCh);
				uint64_t negHash = negVertexHash[i]->hash_prepend(revNextCh);
				if (posHash != negHash)
				{
					return posHash < negHash ? positiveLess : negativeLess;
				}
			}

			return tie;
		}

		static StrandComparisonResult DetermineStrandPrepend(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash, char nextCh, char revNextCh)
		{
			for (size_t i = 0; i < posVertexHash.size(); i++)
			{
				uint64_t posHash = posVertexHash[i]->hash_prepend(nextCh);
				uint64_t negHash = negVertexHash[i]->hash_extend(revNextCh);
				if (posHash != negHash)
				{
					return posHash < negHash ? positiveLess : negativeLess;
				}
			}

			return tie;
		}

		struct VertexRecord
		{
			uint64_t vertexId;
			uint32_t sequenceId;
			uint32_t position;
			VertexRecord() {}
			VertexRecord(uint64_t vertexId, uint32_t sequenceId, uint32_t position) :
				vertexId(vertexId), sequenceId(sequenceId), position(position) {}

			bool operator < (const VertexRecord & record) const
			{
				return std::make_pair(sequenceId, position) < std::make_pair(record.sequenceId, record.position);
			}
		};

		template<class F>
		static bool IsInBloomFilter(const ConcurrentBitVector & filter, std::vector<HashFunctionPtr> & hf, F f, char farg)
		{
			for (size_t i = 0; i < hf.size(); i++)
			{
				uint64_t hvalue = ((*hf[i]).*f)(farg);
				if (!filter.Get(hvalue))
				{
					return false;
				}
			}

			return true;
		}

		static uint64_t NormHash(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash)
		{
			return std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
		}

		static bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}

		static void InitializeHashFunctions(const std::vector<HashFunctionPtr> & seed,
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
					char ch = DnaChar::ReverseChar(*it);
					negEdgeHash[i]->eat(DnaChar::ReverseChar(*it));
				}

				assert(negEdgeHash[i]->hashvalue == negEdgeHash[i]->hash(DnaChar::ReverseCompliment(fragment.substr(offset, length))));
			}
		}

		static void CountEdge(size_t & ordinaryCount, size_t & undefiniteCount, char edgeMark)
		{
			if (DnaChar::IsDefinite(edgeMark))
			{
				++ordinaryCount;
			}
			else
			{
				++undefiniteCount;
			}
		}

		static void CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
			const std::vector<HashFunctionPtr> & hashFunction,
			const ConcurrentBitVector & bitVector,
			size_t vertexLength,
			size_t & totalRecords,
			TaskQueue & taskQueue,
			tpie::file_stream<Occurence> & outFile,
			boost::mutex & outMutex,
			std::unique_ptr<StreamFastaParser::Exception> & error)
		{
			uint64_t low = bound.first;
			uint64_t high = bound.second;
			std::vector<Occurence> output;
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

 					size_t edgeLength = vertexLength + 1;					
					if (task.str.size() >= vertexLength + 2)
					{
						std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
						std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());						
						InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
						size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
						for (size_t pos = 1;; ++pos)
						{
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
							if (Within(std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue), low, high) && definiteCount == vertexLength)
							{
								size_t inCount = 0;
								size_t outCount = 0;
								size_t inUndefCount = 0;
								size_t outUndefCount = 0;

							//	bool x = std::string(task.str.begin() + pos, task.str.begin() + pos + vertexLength) == "AGGTGCAT";
								for (int i = 0; i < DnaChar::LITERAL.size() && inCount < 2 && outCount < 2; i++)
								{
									char nextCh = DnaChar::EXT_LITERAL[i];
									char revNextCh = DnaChar::ReverseChar(nextCh);
									if (nextCh == posPrev)
									{
										CountEdge(inCount, inUndefCount, nextCh);
									}
									else
									{
										StrandComparisonResult result = DetermineStrandPrepend(posVertexHash, negVertexHash, nextCh, revNextCh);
										if (result == positiveLess || result == tie)
										{
											if (IsInBloomFilter(bitVector, posVertexHash, &HashFunction::hash_prepend, nextCh))
											{
												CountEdge(inCount, inUndefCount, nextCh);
											}
										}
										else
										{
											if (IsInBloomFilter(bitVector, negVertexHash, &HashFunction::hash_extend, revNextCh))
											{
												CountEdge(inCount, inUndefCount, nextCh);
											}
										}
									}

									if (nextCh == posExtend)
									{
										CountEdge(outCount, inUndefCount, nextCh);
									}
									else
									{
										StrandComparisonResult result = DetermineStrandExtend(posVertexHash, negVertexHash, nextCh, revNextCh);
										if (result == positiveLess || result == tie)
										{
											if (IsInBloomFilter(bitVector, posVertexHash, &HashFunction::hash_extend, nextCh))
											{
												CountEdge(outCount, outUndefCount, nextCh);
											}
										}
										else 
										{
											if (IsInBloomFilter(bitVector, negVertexHash, &HashFunction::hash_prepend, revNextCh))
											{
												CountEdge(outCount, outUndefCount, nextCh);
											}											
										}
									}
								}

							//	std::cout << std::string(task.str.begin() + pos, task.str.begin() + pos + vertexLength) << std::endl;
							//	std::cout << inCount << ' ' << outCount << std::endl;
								if (inCount > 1 || outCount > 1 || inUndefCount > 0 || outUndefCount > 0)
								{
									output.push_back(Occurence());
									output.back().Set(task.seqId,
										task.start + pos - 1,
										posVertexHash[0]->hashvalue,
										negVertexHash[0]->hashvalue,
										task.str.begin() + pos,
										vertexLength,
										posExtend,
										posPrev);
								}
							}

							if (pos + edgeLength < task.str.size())
							{
								char negExtend = DnaChar::ReverseChar(posExtend);
								char posPrev = task.str[pos];
								char negPrev = DnaChar::ReverseChar(task.str[pos]);
								definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
								for (size_t i = 0; i < hashFunction.size(); i++)
								{
									posVertexHash[i]->update(posPrev, posExtend);
									negVertexHash[i]->reverse_update(negExtend, negPrev);
									assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
									assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
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
						for (const Occurence & str : output)
						{
							if (!outFile.is_writable())
							{
								error.reset(new StreamFastaParser::Exception("Can't write to the temporary file"));
								break;
							}

							outFile.write(str);
						}
						
						output.clear();
						if (error != 0)
						{
							return;
						}						
					}
				}
			}
		}

		static uint64_t FilterFillerWorker(uint64_t low,
			uint64_t high,
			const std::vector<HashFunctionPtr> & hashFunction,
			ConcurrentBitVector & filter,
			size_t edgeLength,
			TaskQueue & taskQueue)
		{
			uint64_t ret = 0;
			std::vector<uint64_t> setup;
			std::vector<uint64_t> hashValue;
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
						hashValue.clear();
						char prevCh = task.str[pos];
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaChar::ReverseChar(nextCh);						
						uint64_t fistMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						StrandComparisonResult result = DetermineStrandExtend(posVertexHash, negVertexHash, nextCh, revNextCh);						

						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							if (result == positiveLess || result == tie)
							{
								hashValue.push_back(posVertexHash[i]->hash_extend(nextCh));
							}
							
							if (result == negativeLess || result == tie)
							{
								hashValue.push_back(negVertexHash[i]->hash_prepend(revNextCh));
							}
						}

						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(revNextCh, DnaChar::ReverseChar(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
						}

						uint64_t secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						if (Within(fistMinHash0, low, high) || Within(secondMinHash0, low, high))
						{
							for (uint64_t value : hashValue)
							{
								setup.push_back(value);
							}
						}
					}
				}

				for (uint64_t hv : setup)
				{
					if (!filter.Get(hv))
					{
						filter.SetConcurrently(hv);
					}
				}

				setup.clear();
			}

			return ret;
		}

		static void InitialFilterFillerWorker(uint64_t binSize,
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
						uint64_t hvalue;
						bool wasSet = true;
						char prevCh = task.str[pos];
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaChar::ReverseChar(nextCh);
						uint64_t firstMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						uint64_t posHash0 = posVertexHash[0]->hash_extend(nextCh);
						uint64_t negHash0 = negVertexHash[0]->hash_prepend(revNextCh);

						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							if (posHash0 < negHash0 || (posHash0 == negHash0 && DnaChar::LessSelfReverseComplement(task.str.begin() + pos, vertexLength)))
							{
								hvalue = posVertexHash[i]->hash_extend(nextCh);
							}
							else
							{
								hvalue = negVertexHash[i]->hash_prepend(revNextCh);
							}

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
							negVertexHash[i]->reverse_update(DnaChar::ReverseChar(nextCh), DnaChar::ReverseChar(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
						}

						uint64_t secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
						if (!wasSet)
						{
							uint64_t value[] = { firstMinHash0, secondMinHash0 };
							for (uint64_t v : value)
							{
								uint64_t bin = v / binSize;
								if (binCounter[bin] < MAX_COUNTER)
								{
									binCounter[bin].fetch_add(1);
								}
							}
						}
					}
				}
			}
		}

		static void DistributeTasks(const std::vector<std::string> & fileName, size_t overlapSize, std::vector<TaskQueuePtr> & taskQueue)
		{
			for (size_t file = 0; file < fileName.size(); file++)
			{
				size_t record = 0;
				const std::string & nowFileName = fileName[file];
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
				{
					char ch;					
					uint64_t prev = 0;
					uint64_t start = 0;
					std::string buf = "N";
					bool over = false;
					do
					{
						over = !parser.GetChar(ch);
						if (!over)
						{
							start++;
							buf.push_back(DnaChar::IsDefinite(ch) ? ch : 'N');
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
										else
										{
											buf.push_back('N');
										}

										q->push(Task(file, prev, over, std::move(buf)));
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

		struct TrueBifurcations
		{
			size_t vertexSize;
			tpie::file_stream<Occurence> * candidate;
			uint64_t falsePositives;			
			std::vector<DnaString> * out;
			boost::mutex * outMutex;
			TrueBifurcations(tpie::file_stream<Occurence> * candidate,
				std::vector<DnaString> * out,
				boost::mutex * outMutex,
				size_t vertexSize) :
				candidate(candidate), out(out), vertexSize(vertexSize), outMutex(outMutex), falsePositives(0) {}

			uint64_t operator()(uint64_t records, uint64_t init, uint64_t & verticesCount, tpie::file_stream<VertexRecord> & outFile) const
			{		
				Occurence base;
				Occurence next;	
				std::vector<Occurence> store;
				uint64_t falsePositives = init;				
				for (uint64_t nowRecord = 0; nowRecord < records;)
				{					
					if (nowRecord == 0)
					{
						++nowRecord;
						base = candidate->read();
					}
					else
					{
						base = next;
					}

					store.clear();
					size_t inUnknownCount = 0;
					size_t outUnknownCount = 0;
					bool bifurcation = false;					
					bool selfReverseCompliment = base.IsSelfReverseCompliment(vertexSize);
					for (next = base; ;)
					{						
						if (!base.EqualBase(next))
						{
							break;
						}

						//bool x = base.GetBase().ToString(vertexSize) == "AGGTGCAT" || base.GetBase().ReverseComplement(vertexSize).ToString(vertexSize) == "AGGTGCAT";
						inUnknownCount += DnaChar::IsDefinite(next.Prev()) ? 0 : 1;
						outUnknownCount += DnaChar::IsDefinite(next.Next()) ? 0 : 1;
						store.push_back(next);
						if (!bifurcation)
						{
							bifurcation = base.Prev() != next.Prev() || base.Next() != next.Next();
							if (selfReverseCompliment)
							{
								inUnknownCount += DnaChar::IsDefinite(next.Next()) ? 0 : 1;
								outUnknownCount += DnaChar::IsDefinite(next.Prev()) ? 0 : 1;
								bifurcation = bifurcation ||
									base.Prev() != DnaChar::ReverseChar(next.Next()) ||
									base.Next() != DnaChar::ReverseChar(next.Prev());
							} 
						}

						if (nowRecord < records)
						{
							++nowRecord;
							next = candidate->read();
						}
						else
						{
							break;
						}
					}

					if (bifurcation || inUnknownCount > 1 || outUnknownCount > 1)
					{
						boost::lock_guard<boost::mutex> guard(*outMutex);
						out->push_back(base.GetBase());
						for (const Occurence & record : store)
						{
							VertexRecord vertex(verticesCount, record.GetSequenceId(), record.GetPosition());
							outFile.write(vertex);
						}

						++verticesCount;
					}
					else
					{
						falsePositives++;
					}
				}

				return falsePositives;
			}
		};

		class DnaStringLess
		{
		public:
			DnaStringLess(size_t vertexSize) : vertexSize_(vertexSize)
			{

			}

			bool operator() (const DnaString & v1, const DnaString & v2) const
			{
				return DnaString::LessPrefix(v1, v2, vertexSize_);
			}

		private:
			size_t vertexSize_;
		};
		

		size_t vertexSize_;
		std::vector<DnaString> bifurcation_;
	};
}

#endif	