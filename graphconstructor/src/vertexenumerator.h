#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#define MAX_CAPACITY 5

#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>

#include <tbb/spin_rw_mutex.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>

#include <boost/ref.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>
#include <boost/dynamic_bitset.hpp>
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
		static const size_t BUF_SIZE = 1 << 24;

		typedef CompressedString<CAPACITY> DnaString;
		typedef CandidateOccurence<CAPACITY> Occurence;

		class OccurenceHash
		{
		public:
			uint64_t operator()(const Occurence & occurence) const
			{
				return occurence.Hash();
			}
		};

		class DnaStringHash
		{
		public:
			uint64_t operator()(const DnaString & dnaString) const
			{
				return dnaString.Hash();
			}
		};

		class OccurenceEquality
		{
		public:
			bool operator()(const Occurence & occurence1, const Occurence & occurence2) const
			{
				return occurence1.EqualBase(occurence2);
			}
		};

		typedef std::unordered_map<DnaString, uint64_t, DnaStringHash> BifurcationMap;
		typedef tbb::concurrent_unordered_multiset<Occurence, OccurenceHash, OccurenceEquality> OccurenceSet;

	public:

		size_t GetVerticesCount() const
		{
			return bifurcationKey_.size();
		}

		size_t GetId(const std::string & vertex) const
		{
			DnaString str;
			str.CopyFromString(vertex.begin(), vertexSize_);
			typename BifurcationMap::const_iterator it = bifurcationKey_.find(str);
			if (it != bifurcationKey_.end())
			{
				return it->second;
			}

			return INVALID_VERTEX;
		}

		virtual void Dump(std::vector<std::string> & toCopy) const
		{/*
			std::set<std::string> ret;
			for (const DnaString & str : bifurcation_)
			{
				std::string pos = str.ToString(vertexSize_);
				std::string neg = DnaChar::ReverseCompliment(pos);
				assert(ret.count(pos) == 0 && ret.count(neg) == 0);
				ret.insert(pos);
				ret.insert(neg);
			}

			toCopy.assign(ret.begin(), ret.end());*/
		}

		VertexEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			size_t aggregationThreads,
			const std::string & tmpDirName,
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
			std::vector<TaskQueuePtr> taskQueue(threads);
			std::vector<boost::thread> workerThread(threads);
			std::atomic<uint32_t> * binCounter = new std::atomic<uint32_t>[BINS_COUNT];
			{
				std::fill(binCounter, binCounter + BINS_COUNT, 0);
				ConcurrentBitVector bitVector(realSize);				
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					taskQueue[i] = TaskQueuePtr(new TaskQueue(QUEUE_CAPACITY));
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

			std::cout << "Round size = " << realSize / roundSize << std::endl;
			std::cout << std::string(80, '-') << std::endl;
			uint64_t low = 0;
			uint64_t high = 0;	
			size_t lowBoundary = 0;
			uint64_t totalFpCount = 0;
			uint64_t verticesCount = 0;
			std::ofstream bifurcationTempWrite((tmpDirName + "/bifurcations.bin").c_str(), ios::binary);
			if (!bifurcationTempWrite)
			{
				throw StreamFastaParser::Exception("Can't create a temp file");
			}

			time_t mark;
			std::unique_ptr<StreamFastaParser::Exception> error;
			for (size_t round = 0; round < rounds; round++)
			{
				mark = time(0);
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

				std::vector<TaskQueuePtr> taskQueue;
				std::vector<boost::thread> workerThread(threads);
				std::cout << "Ratio = " << double(realSize) / accumulated << std::endl;
				uint64_t high = lowBoundary * BIN_SIZE;
				{
					ConcurrentBitVector bitVector(realSize);
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
					std::unique_ptr<StreamFastaParser::Exception> error;
					for (size_t i = 0; i < workerThread.size(); i++)
					{
						workerThread[i] = boost::thread(CandidateCheckingWorker,
							std::make_pair(low, high),
							boost::cref(hashFunction),
							boost::cref(bitVector),
							vertexLength,
							boost::ref(*taskQueue[i]),
							boost::cref(tmpDirName),
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
				tbb::spin_rw_mutex mutex;
				OccurenceSet occurenceSet(1 << 20);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					workerThread[i] = boost::thread(CandidateFilteringWorker,
						boost::cref(hashFunction),
						vertexLength,
						boost::ref(*taskQueue[i]),
						boost::ref(occurenceSet),
						boost::ref(mutex),
						boost::cref(tmpDirName),
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

				size_t falsePositives = 0;
				size_t truePositives = TrueBifurcations(occurenceSet, bifurcationTempWrite, vertexSize_, falsePositives);
				std::cout << time(0) - mark << std::endl;				
				std::cout << "Vertex count = " << truePositives << std::endl;
				std::cout << "FP count = " << falsePositives << std::endl;
				std::cout << "Records = " << occurenceSet.size() << std::endl;
				std::cout << std::string(80, '-') << std::endl;
				totalFpCount += falsePositives;
				verticesCount += truePositives;
				low = high + 1;
			}
			
			delete[] binCounter;
			bifurcationTempWrite.close();
			bifurcationKey_.swap(BifurcationMap(verticesCount));
			std::ifstream bifurcationTempRead((tmpDirName + "/bifurcations.bin").c_str(), ios::binary);
			if (!bifurcationTempRead)
			{
				throw StreamFastaParser::Exception("Can't open the temp file");
			}

			uint64_t bitsPower = 0;
			while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
			{
				++bitsPower;
			}

			hashFunction.clear();
			std::vector<bool> bifurcationFilter(uint64_t(1) << bitsPower);
			hashFunction.resize(double(bifurcationFilter.size()) / verticesCount * 0.7);
			for (HashFunctionPtr & ptr : hashFunction)
			{
				ptr = HashFunctionPtr(new HashFunction(vertexLength, bitsPower));
			}

			DnaString buf;
			std::string stringBuf(vertexLength, ' ');
			for (size_t i = 0; i < verticesCount; i++)
			{
				size_t key = bifurcationKey_.size();
				buf.ReadFromFile(bifurcationTempRead);
				bifurcationKey_[buf] = key;
				buf.ToString(stringBuf, vertexLength);
				for (HashFunctionPtr & ptr : hashFunction)
				{
					bifurcationFilter[ptr->hash(stringBuf)] = true;
				}
			}

			std::cout << "Total FPs = " << totalFpCount << std::endl;
			std::ofstream compressedDbg(outFileName.c_str());
			if (!compressedDbg)
			{
				throw StreamFastaParser::Exception("Can't create the output file");
			}

			mark = time(0);
			std::atomic<uint32_t> currentPiece = 0;
			for (size_t i = 0; i < workerThread.size(); i++)
			{
				workerThread[i] = boost::thread(EdgeConstructionWorker,
					boost::ref(hashFunction),
					vertexLength,
					boost::ref(*taskQueue[i]),
					boost::ref(bifurcationKey_),
					boost::ref(bifurcationFilter),
					boost::ref(compressedDbg),
					boost::ref(currentPiece));
			}

			DistributeTasks(fileName, vertexLength + 1, taskQueue);
			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				workerThread[i].join();
			}

			std::cout << "Edges construction: " << time(0) - mark << std::endl;
			std::cout << std::string(80, '-') << std::endl;
		}

	private:

		static const size_t QUEUE_CAPACITY = 16;
		static const uint64_t BINS_COUNT = 1 << 24;
		static const uint32_t MAX_COUNTER = UINT32_MAX >> 1;

		struct Task
		{			
			bool isFinal;
			uint32_t piece;
			uint64_t start;
			uint64_t seqId;
			std::string str;
#ifdef _DEBUG
			static const size_t TASK_SIZE = 32;
#else
			static const size_t TASK_SIZE = 1 << 18;
#endif					
			static const size_t GAME_OVER = SIZE_MAX;
			Task() {}
			Task(uint64_t seqId, uint64_t start, uint32_t piece, bool isFinal, std::string && str) : 
				seqId(seqId), start(start), piece(piece), isFinal(isFinal), str(std::move(str)) {}
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

		static StrandComparisonResult DetermineStrandPrepend(const std::vector<HashFunctionPtr> & posVertexHash, const std::vector<HashFunctionPtr> & negVertexHash, char prevCh, char revPrevCh)
		{
			for (size_t i = 0; i < posVertexHash.size(); i++)
			{
				uint64_t posHash = posVertexHash[i]->hash_prepend(prevCh);
				uint64_t negHash = negVertexHash[i]->hash_extend(revPrevCh);
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

		static bool IsInBloomFilterExtend(const ConcurrentBitVector & filter, std::vector<HashFunctionPtr> & hf, char farg)
		{
			for (size_t i = 0; i < hf.size(); i++)
			{
				uint64_t hvalue = hf[i]->hash_extend(farg);
				if (!filter.Get(hvalue))
				{
					return false;
				}
			}

			return true;
		}

		static bool IsInBloomFilterPrepend(const ConcurrentBitVector & filter, std::vector<HashFunctionPtr> & hf, char farg)
		{
			for (size_t i = 0; i < hf.size(); i++)
			{
				uint64_t hvalue = hf[i]->hash_prepend(farg);
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

		static std::string CandidateMaskFileName(const std::string & directory, size_t sequence, size_t pos)
		{
			std::stringstream ss;
			ss << directory << "/" << sequence << "_" << pos << ".tmp";
			return ss.str();
		}

		static void CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
			const std::vector<HashFunctionPtr> & hashFunction,
			const ConcurrentBitVector & bitVector,
			size_t vertexLength,
			TaskQueue & taskQueue,
			const std::string & tmpDirectory,
			std::unique_ptr<StreamFastaParser::Exception> & error)
		{
			uint64_t low = bound.first;
			uint64_t high = bound.second;
			typedef uint32_t BITSET_BLOCK_TYPE;
			boost::dynamic_bitset<BITSET_BLOCK_TYPE> candidateMask(Task::TASK_SIZE);
			std::vector<BITSET_BLOCK_TYPE> buf(candidateMask.num_blocks(), 0);
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
					
 					size_t edgeLength = vertexLength + 1;					
					if (task.str.size() >= vertexLength + 2)
					{
						candidateMask.reset();
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
								size_t inCount = DnaChar::IsDefinite(posPrev) ? 0 : 2;
								size_t outCount = DnaChar::IsDefinite(posExtend) ? 0 : 2;
								for (int i = 0; i < DnaChar::LITERAL.size() && inCount < 2 && outCount < 2; i++)
								{
									char nextCh = DnaChar::LITERAL[i];
									char revNextCh = DnaChar::ReverseChar(nextCh);
									if (nextCh == posPrev)
									{
										++inCount;
									}
									else
									{
										StrandComparisonResult result = DetermineStrandPrepend(posVertexHash, negVertexHash, nextCh, revNextCh);
										if (result == positiveLess || result == tie)
										{
											if (IsInBloomFilterPrepend(bitVector, posVertexHash, nextCh))
											{
												++inCount;
											}
										}
										else
										{
											if (IsInBloomFilterExtend(bitVector, negVertexHash, revNextCh))
											{
												++inCount;
											}
										}
									}

									if (nextCh == posExtend)
									{
										++outCount;
									}
									else 
									{
										StrandComparisonResult result = DetermineStrandExtend(posVertexHash, negVertexHash, nextCh, revNextCh);
										if (result == positiveLess || result == tie)
										{
											if (IsInBloomFilterExtend(bitVector, posVertexHash,  nextCh))
											{
												++outCount;
											}
										}
										else 
										{
											if (IsInBloomFilterPrepend(bitVector, negVertexHash, revNextCh))
											{
												++outCount;
											}							
										}
									}
								}

								if (inCount > 1 || outCount > 1)
								{
									candidateMask.set(pos);
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

						std::ofstream candidateMaskFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start).c_str(), std::ios::binary);
						boost::to_block_range(candidateMask, buf.begin());
						candidateMaskFile.write(reinterpret_cast<const char*>(&buf[0]), buf.size() * sizeof(BITSET_BLOCK_TYPE));
						candidateMaskFile.close();
					}
				}
			}
		}

		static void CandidateFilteringWorker(const std::vector<HashFunctionPtr> & hashFunction,
			size_t vertexLength,
			TaskQueue & taskQueue,
			OccurenceSet & occurenceSet,
			tbb::spin_rw_mutex & mutex,
			const std::string & tmpDirectory,
			std::unique_ptr<StreamFastaParser::Exception> & error)
		{
			typedef uint32_t BITSET_BLOCK_TYPE;
			boost::dynamic_bitset<BITSET_BLOCK_TYPE> candidateMask(Task::TASK_SIZE);
			std::vector<BITSET_BLOCK_TYPE> buf(candidateMask.num_blocks(), 0);
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

					std::vector<HashFunctionPtr> posVertexHash(1);
					std::vector<HashFunctionPtr> negVertexHash(1);
					size_t edgeLength = vertexLength + 1;
					if (task.str.size() >= vertexLength + 2)
					{
						InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
						{
							std::ifstream candidateMaskFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start).c_str(), std::ios::binary);
							candidateMaskFile.read(reinterpret_cast<char*>(&buf[0]), buf.size() * sizeof(BITSET_BLOCK_TYPE));
							boost::from_block_range(buf.begin(), buf.end(), candidateMask);
						}
						
						for (size_t pos = 1;; ++pos)
						{
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							if (candidateMask.test(pos))
							{
								Occurence now;
								now.Set(posVertexHash[0]->hashvalue,
									negVertexHash[0]->hashvalue,
									task.str.begin() + pos,
									vertexLength,
									posExtend,
									posPrev,
									false);
								
								size_t count = 0;
								size_t inUnknownCount = now.Prev() == 'N' ? 1 : 0;
								size_t outUnknownCount = now.Next() == 'N' ? 1 : 0;
								bool newBifurcation = false;
								bool alreadyBifurcation = false;
								mutex.lock_read();
								auto range = occurenceSet.equal_range(now);
								for (auto it = range.first; it != range.second; ++it)
								{
									++count;
									inUnknownCount += DnaChar::IsDefinite(it->Prev()) ? 0 : 1;
									outUnknownCount += DnaChar::IsDefinite(it->Next()) ? 0 : 1;
									if (!alreadyBifurcation && it->IsBifurcation())
									{
										alreadyBifurcation = true;
									}

									if (it->Next() != now.Next() || it->Prev() != now.Prev() || inUnknownCount > 1 || outUnknownCount > 1)
									{
										newBifurcation = true;
									}
								}

								if (count == 0)
								{
									occurenceSet.insert(now);
								}
								else
								{
									if ((newBifurcation && !alreadyBifurcation) || (alreadyBifurcation && count > 1))
									{
										now.MakeBifurcation();
										mutex.unlock();
										mutex.lock();
										range = occurenceSet.equal_range(now);
										occurenceSet.unsafe_erase(range.first, range.second);
										occurenceSet.insert(now);
									}
								}

								mutex.unlock();
							}

							if (pos + edgeLength < task.str.size())
							{
								char negExtend = DnaChar::ReverseChar(posExtend);
								char posPrev = task.str[pos];
								char negPrev = DnaChar::ReverseChar(task.str[pos]);
								for (size_t i = 0; i < 1; i++)
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
				}
			}
		}

		struct EdgeResult
		{
			uint32_t seqId;
			uint32_t pieceId;			
			std::vector<uint64_t> pos;			
			std::vector<uint64_t> bifId;
		};

		static void EdgeConstructionWorker(const std::vector<HashFunctionPtr> & hashFunction,
			size_t vertexLength,
			TaskQueue & taskQueue,			
			BifurcationMap & bifurcationKey,
			std::vector<bool> & bifurcationFilter,
			std::ofstream & outFile,
			std::atomic<uint32_t> & currentPiece)
		{
			DnaString bitBuf;
			std::deque<EdgeResult> result;
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

					std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
					std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
					size_t edgeLength = vertexLength + 1;
					if (task.str.size() >= vertexLength + 2)
					{
						result.push_back(EdgeResult());
						result.back().seqId = task.seqId;
						result.back().pieceId = task.piece;
						InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength, 1);
						for (size_t pos = 1;; ++pos)
						{							
							bool posFound = true;
							bool negFound = true;
							char posPrev = task.str[pos - 1];
							char posExtend = task.str[pos + vertexLength];
							for (size_t i = 0; i < posVertexHash.size() && (posFound || negFound); i++)
							{
								if (!bifurcationFilter[posVertexHash[i]->hashvalue])
								{
									posFound = false;
								}

								if (!bifurcationFilter[negVertexHash[i]->hashvalue])
								{
									negFound = false;
								}
							}

							if (posFound)
							{
								posFound = false;
								bitBuf.Clear();
								bitBuf.CopyFromString(task.str.begin() + pos, vertexLength);
								auto it = bifurcationKey.find(bitBuf);
								if (it != bifurcationKey.end())
								{
									posFound = true;
									result.back().pos.push_back(task.start + pos - 1);
									result.back().bifId.push_back(it->second);
								}
								
							}

							if (!posFound && negFound)
							{
								bitBuf.Clear();
								bitBuf.CopyFromReverseString(task.str.begin() + pos, vertexLength);
								auto it = bifurcationKey.find(bitBuf);
								if (it != bifurcationKey.end())
								{
									result.back().pos.push_back(task.start + pos - 1);
									result.back().bifId.push_back(it->second);
								}
							}

							if (pos + edgeLength < task.str.size())
							{
								char negExtend = DnaChar::ReverseChar(posExtend);
								char posPrev = task.str[pos];
								char negPrev = DnaChar::ReverseChar(task.str[pos]);
								for (size_t i = 0; i < 1; i++)
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

						while (result.size() > 0 && result.front().pieceId == currentPiece)
						{
							for (size_t i = 0; i < result.front().pos.size(); i++)
							{
								outFile << result.front().seqId << " " << result.front().pos[i] << " " << result.front().bifId[i] << std::endl;
							}

							++currentPiece;
							result.pop_front();
						}
					}
				}
			}

			while (result.size() > 0 && result.front().pieceId == currentPiece)
			{
				for (size_t i = 0; i < result.front().pos.size(); i++)
				{
					outFile << result.front().pos[i] << " " << result.front().bifId[i] << std::endl;
				}

				++currentPiece;
				result.pop_front();
			}
		}

		static void PutInBloomFilterExtend(const std::vector<HashFunctionPtr> & posVertexHash,
			const std::vector<HashFunctionPtr> & negVertexHash,
			char nextCh,
			char revNextCh,
			std::vector<uint64_t> & hashValue)
		{
			StrandComparisonResult result = DetermineStrandExtend(posVertexHash, negVertexHash, nextCh, revNextCh);
			for (size_t i = 0; i < posVertexHash.size(); i++)
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
		}

		static void PutInBloomFilterPrepend(const std::vector<HashFunctionPtr> & posVertexHash,
			const std::vector<HashFunctionPtr> & negVertexHash,
			char prevCh,
			char revPrevCh,
			std::vector<uint64_t> & hashValue)
		{
			StrandComparisonResult result = DetermineStrandPrepend(posVertexHash, negVertexHash, prevCh, revPrevCh);
			for (size_t i = 0; i < posVertexHash.size(); i++)
			{
				if (result == positiveLess || result == tie)
				{
					hashValue.push_back(posVertexHash[i]->hash_prepend(prevCh));
				}

				if (result == negativeLess || result == tie)
				{
					hashValue.push_back(negVertexHash[i]->hash_extend(revPrevCh));
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
			const char DUMMY_CHAR = DnaChar::LITERAL[0];
			const char REV_DUMMY_CHAR = DnaChar::ReverseChar(DUMMY_CHAR);
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

					uint64_t fistMinHash0;
					uint64_t secondMinHash0;
					size_t vertexLength = edgeLength - 1;
					std::vector<HashFunctionPtr> posVertexHash(hashFunction.size());
					std::vector<HashFunctionPtr> negVertexHash(hashFunction.size());
					size_t definiteCount = std::count_if(task.str.begin(), task.str.begin() + vertexLength, DnaChar::IsDefinite);
					InitializeHashFunctions(hashFunction, posVertexHash, negVertexHash, task.str, vertexLength);
					for (size_t pos = 0; ; ++pos)
					{
						hashValue.clear();
						char prevCh = task.str[pos];
						char nextCh = task.str[pos + edgeLength - 1];
						char revNextCh = DnaChar::ReverseChar(nextCh);
						assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
						if (definiteCount == vertexLength)
						{
							fistMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
							if (DnaChar::IsDefinite(nextCh))
							{
								PutInBloomFilterExtend(posVertexHash, negVertexHash, nextCh, revNextCh, hashValue);
							}
							else
							{
								PutInBloomFilterExtend(posVertexHash, negVertexHash, DUMMY_CHAR, REV_DUMMY_CHAR, hashValue);
								PutInBloomFilterExtend(posVertexHash, negVertexHash, REV_DUMMY_CHAR, DUMMY_CHAR, hashValue);
							}

							if (pos > 0 && !DnaChar::IsDefinite(task.str[pos - 1]))
							{
								PutInBloomFilterPrepend(posVertexHash, negVertexHash, DUMMY_CHAR, REV_DUMMY_CHAR, hashValue);
								PutInBloomFilterPrepend(posVertexHash, negVertexHash, REV_DUMMY_CHAR, DUMMY_CHAR, hashValue);
							}
						}

						for (size_t i = 0; i < hashFunction.size(); i++)
						{
							posVertexHash[i]->update(prevCh, nextCh);
							assert(posVertexHash[i]->hashvalue == posVertexHash[i]->hash(task.str.substr(pos + 1, vertexLength)));
							negVertexHash[i]->reverse_update(revNextCh, DnaChar::ReverseChar(prevCh));
							assert(negVertexHash[i]->hashvalue == negVertexHash[i]->hash(DnaChar::ReverseCompliment(task.str.substr(pos + 1, vertexLength))));
						}

						if (definiteCount == vertexLength)
						{
							secondMinHash0 = std::min(posVertexHash[0]->hashvalue, negVertexHash[0]->hashvalue);
							if (Within(fistMinHash0, low, high) || Within(secondMinHash0, low, high))
							{
								for (uint64_t value : hashValue)
								{
									setup.push_back(value);
								}
							}
						}

						if (pos + vertexLength < task.str.size() - 1)
						{
							definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(prevCh) ? 1 : 0);
						}
						else
						{
							break;
						}
					}
				}

				for (uint64_t hashValue : setup)
				{
					if (!filter.Get(hashValue))
					{
						filter.SetConcurrently(hashValue);
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
			uint32_t pieceCount = 0;
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

										q->push(Task(record, prev, pieceCount++, over, std::move(buf)));
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

				taskQueue[i]->push(Task(0, Task::GAME_OVER, 0, true, std::string()));
			}
		}
		
		uint64_t TrueBifurcations(const OccurenceSet & occurenceSet, std::ofstream & out, size_t vertexSize, size_t & falsePositives) const
		{			
			uint64_t truePositives = falsePositives = 0;
			std::vector<Occurence> store;			
			for (auto it = occurenceSet.begin(); it != occurenceSet.end(); )
			{
				Occurence base = *it;
				size_t inUnknownCount = 0;
				size_t outUnknownCount = 0;
				bool bifurcation = it->IsBifurcation();
				bool selfReverseCompliment = base.IsSelfReverseCompliment(vertexSize);
				auto jt = it;
				for (; jt != occurenceSet.end(); ++jt)
				{
					Occurence next = *jt;
					if (!base.EqualBase(next))
					{
						break;
					}

					inUnknownCount += DnaChar::IsDefinite(next.Prev()) ? 0 : 1;
					outUnknownCount += DnaChar::IsDefinite(next.Next()) ? 0 : 1;					
					if (!bifurcation)
					{
						bifurcation = jt->IsBifurcation() || base.Prev() != next.Prev() || base.Next() != next.Next();
						if (selfReverseCompliment)
						{
							inUnknownCount += DnaChar::IsDefinite(next.Next()) ? 0 : 1;
							outUnknownCount += DnaChar::IsDefinite(next.Prev()) ? 0 : 1;
							bifurcation = bifurcation ||
								base.Prev() != DnaChar::ReverseChar(next.Next()) ||
								base.Next() != DnaChar::ReverseChar(next.Prev());
						}
					}
				}

				if (bifurcation || inUnknownCount > 1 || outUnknownCount > 1)
				{
					++truePositives;
					it->GetBase().WriteToFile(out);
				}
				else
				{
					falsePositives++;
				}

				it = jt;
			}

			return truePositives;
		}

		size_t vertexSize_;
		BifurcationMap bifurcationKey_;
	};
}

#endif	
