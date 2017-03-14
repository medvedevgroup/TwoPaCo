#ifndef _VERTEX_ENUMERATOR_H_
#define _VERTEX_ENUMERATOR_H_

#define MAX_CAPACITY 10

#include <deque>
#include <cstdio>
#include <numeric>
#include <sstream>
#include <unordered_map>

#include <tbb/tbb.h>
#include <tbb/mutex.h>
#include <tbb/compat/thread>
#include <tbb/spin_rw_mutex.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_unordered_set.h>

#include <junctionapi/junctionapi.h>

#include "vertexrollinghash.h"
#include "streamfastaparser.h"
#include "bifurcationstorage.h"
#include "candidateoccurence.h"
#include "concurrentbitvector.h"

namespace TwoPaCo
{
	class VertexEnumerator
	{
	public:
		virtual size_t GetVerticesCount() const = 0;
		virtual int64_t GetId(const std::string & vertex) const = 0;

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
		const std::string & tmpFileName,
		const std::string & outFileName);

	template<size_t CAPACITY>
	class VertexEnumeratorImpl : public VertexEnumerator
	{
	private:
		static const size_t BUF_SIZE = 1 << 24;
		BifurcationStorage<CAPACITY> bifStorage_;
		typedef CompressedString<CAPACITY> DnaString;
		typedef CandidateOccurence<CAPACITY> Occurence;

		class FilterFillerWorker;

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

		typedef tbb::concurrent_unordered_set<Occurence, OccurenceHash, OccurenceEquality> OccurenceSet;

	public:

		int64_t GetId(const std::string & vertex) const
		{
			return bifStorage_.GetId(vertex.begin());
		}

		size_t GetVerticesCount() const
		{
			return bifStorage_.GetDistinctVerticesCount();
		}


		VertexEnumeratorImpl(const std::vector<std::string> & fileName,
			size_t vertexLength,
			size_t filterSize,
			size_t hashFunctions,
			size_t rounds,
			size_t threads,
			const std::string & tmpDirName,
			const std::string & outFileNamePrefix) :
			vertexSize_(vertexLength)
		{
			uint64_t realSize = uint64_t(1) << filterSize;
			std::cout << "Threads = " << threads << std::endl;
			std::cout << "Vertex length = " << vertexLength << std::endl;
			std::cout << "Hash functions = " << hashFunctions << std::endl;
			std::cout << "Filter size = " << realSize << std::endl;
			std::cout << "Capacity = " << CAPACITY << std::endl;
			std::cout << "Files: " << std::endl;
			for (const std::string & fn : fileName)
			{
				std::cout << fn << std::endl;
			}
#ifdef LOGGING
			std::ofstream logFile((tmpDirName + "/log.txt").c_str());
			if (!logFile)
			{
				throw StreamFastaParser::Exception("Can't open the log file");
			}
#else
			std::ostream & logFile = std::cerr;
#endif

			tbb::mutex errorMutex;
			std::unique_ptr<std::runtime_error> error;

			VertexRollingHashSeed hashFunctionSeed(hashFunctions, vertexLength, filterSize);
			size_t edgeLength = vertexLength + 1;
			std::vector<TaskQueuePtr> taskQueue(threads);
			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				taskQueue[i].reset(new TaskQueue());
				taskQueue[i]->set_capacity(QUEUE_CAPACITY);
			}

			std::cout << std::string(80, '-') << std::endl;
			uint64_t low = 0;
			uint64_t high = realSize;
			size_t lowBoundary = 0;
			uint64_t totalFpCount = 0;
			uint64_t verticesCount = 0;
			std::ofstream bifurcationTempWrite((tmpDirName + "/bifurcations.bin").c_str(), ios::binary);
			if (!bifurcationTempWrite)
			{
				throw StreamFastaParser::Exception("Can't create a temp file");
			}

			rounds = 1;
			time_t mark;			
			for (size_t round = 0; round < rounds; round++)
			{
				std::atomic<uint64_t> marks;
				marks = 0;
				mark = time(0);

				{
					ConcurrentBitVector bitVector(realSize);
					std::cout << "Round " << round << ", " << low << ":" << high << std::endl;
					std::cout << "Pass\tFilling\tFiltering" << std::endl << "1\t";
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							FilterFillerWorker worker(low,
								high,
								std::cref(hashFunctionSeed),
								std::ref(bitVector),
								edgeLength,
								std::ref(*taskQueue[i]));
							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, edgeLength, taskQueue, error, errorMutex, logFile);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							workerThread[i]->join();
						}
					}

					std::cout << time(0) - mark << "\t";
					mark = time(0);
					{
						std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
						for (size_t i = 0; i < workerThread.size(); i++)
						{
							CandidateCheckingWorker worker(std::make_pair(low, high),
								hashFunctionSeed,
								bitVector,
								vertexLength,
								*taskQueue[i],
								tmpDirName,
								marks,
								round,
								error,
								errorMutex);

							workerThread[i].reset(new tbb::tbb_thread(worker));
						}

						DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile);
						for (size_t i = 0; i < taskQueue.size(); i++)
						{
							workerThread[i]->join();
						}

						if (error != 0)
						{
							throw *error;
						}
					}

					std::cout << time(0) - mark << "\t" << std::endl;
				}

				mark = time(0);
				tbb::spin_rw_mutex mutex;
				std::cout << "2\t";
				OccurenceSet occurenceSet(1 << 20);
				{
					std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
					for (size_t i = 0; i < workerThread.size(); i++)
					{
						CandidateFinalFilteringWorker worker(hashFunctionSeed,
							vertexLength,
							*taskQueue[i],
							occurenceSet,
							mutex,
							tmpDirName,
							round,
							error,
							errorMutex);

						workerThread[i].reset(new tbb::tbb_thread(worker));
					}

					DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile);
					for (size_t i = 0; i < taskQueue.size(); i++)
					{
						workerThread[i]->join();
					}

					if (error != 0)
					{
						throw std::runtime_error(*error);
					}

					std::cout << time(0) - mark << "\t";
				}

				mark = time(0);
				size_t falsePositives = 0;
				size_t truePositives = TrueBifurcations(occurenceSet, bifurcationTempWrite, vertexSize_, falsePositives);
				std::cout << time(0) - mark << std::endl;
				std::cout << "True junctions count = " << truePositives << std::endl;
				std::cout << "False junctions count = " << falsePositives << std::endl;
				std::cout << "Hash table size = " << occurenceSet.size() << std::endl;
				std::cout << "Candidate marks count = " << marks << std::endl;
				std::cout << std::string(80, '-') << std::endl;
				totalFpCount += falsePositives;
				verticesCount += truePositives;
				low = high + 1;
			}

			mark = time(0);			
			std::string bifurcationTempReadName = (tmpDirName + "/bifurcations.bin");
			bifurcationTempWrite.close();
			{
				std::ifstream bifurcationTempRead(bifurcationTempReadName.c_str(), ios::binary);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't open the temp file");
				}

				bifStorage_.Init(bifurcationTempRead, verticesCount, vertexLength, threads);
			}

			std::remove(bifurcationTempReadName.c_str());
			std::cout << "Reallocating bifurcations time: " << time(0) - mark << std::endl;

			mark = time(0);			
			std::atomic<uint64_t> occurence;
			tbb::mutex currentStubVertexMutex;
			std::atomic<uint64_t> currentPiece;			
			uint64_t currentStubVertexId = verticesCount + 42;
			JunctionPositionWriter posWriter(outFileNamePrefix);
			occurence = currentPiece = 0;
			{
				std::vector<std::unique_ptr<tbb::tbb_thread> > workerThread(threads);
				for (size_t i = 0; i < workerThread.size(); i++)
				{
					EdgeConstructionWorker worker(vertexLength,
						*taskQueue[i],
						bifStorage_,
						posWriter,
						currentPiece,
						occurence,
						currentStubVertexId,
						currentStubVertexMutex,
						tmpDirName,
						rounds,
						error,
						errorMutex);

					workerThread[i].reset(new tbb::tbb_thread(worker));
				}

				DistributeTasks(fileName, vertexLength + 1, taskQueue, error, errorMutex, logFile);
				for (size_t i = 0; i < taskQueue.size(); i++)
				{
					workerThread[i]->join();
				}
			}

			if (error != 0)
			{
				throw std::runtime_error(*error);
			}

			std::cout << "True marks count: " << occurence << std::endl;
			std::cout << "Edges construction time: " << time(0) - mark << std::endl;
			std::cout << std::string(80, '-') << std::endl;
		}

	private:

		static const size_t QUEUE_CAPACITY = 16;
		static const uint64_t BINS_COUNT = 1 << 24;

		static bool IsOutgoingEdgeInBloomFilter(const ConcurrentBitVector & filter, const VertexRollingHash & hf, char farg)
		{
			std::vector<uint64_t> value;
			value.clear();
			hf.GetOutgoingEdgeHash(farg, value);
			for (size_t i = 0; i < value.size(); i++)
			{
				uint64_t hvalue = value[i];
				if (!filter.GetBit(hvalue))
				{
					return false;
				}
			}

			return true;
		}

		static bool IsIngoingEdgeInBloomFilter(const ConcurrentBitVector & filter, const VertexRollingHash & hf, char farg)
		{
			std::vector<uint64_t> value;
			value.clear();
			hf.GetIngoingEdgeHash(farg, value);
			for (size_t i = 0; i < value.size(); i++)
			{
				uint64_t hvalue = value[i];
				if (!filter.GetBit(hvalue))
				{
					return false;
				}
			}

			return true;
		}

		static bool Within(uint64_t hvalue, uint64_t low, uint64_t high)
		{
			return hvalue >= low && hvalue <= high;
		}

		static std::string CandidateMaskFileName(const std::string & directory, size_t sequence, size_t pos, size_t round)
		{
			std::stringstream ss;
			ss << directory << "/" << sequence << "_" << pos << "_" << round << ".tmp";
			return ss.str();
		}

		static void ReportError(tbb::mutex & errorMutex, std::unique_ptr<std::runtime_error> & error, const std::string & msg)
		{
			errorMutex.lock();
			if (error == 0)
			{
				error.reset(new StreamFastaParser::Exception(msg));
			}

			errorMutex.unlock();
		}

		class CandidateCheckingWorker
		{
		public:
			CandidateCheckingWorker(std::pair<uint64_t, uint64_t> bound,
				const VertexRollingHashSeed & hashFunction,
				const ConcurrentBitVector & bitVector,
				size_t vertexLength,
				TaskQueue & taskQueue,
				const std::string & tmpDirectory,
				std::atomic<uint64_t> & marksCount,
				size_t round,
				std::unique_ptr<std::runtime_error> & error,				
				tbb::mutex & errorMutex) : bound(bound), hashFunction(hashFunction), bitVector(bitVector), vertexLength(vertexLength), taskQueue(taskQueue),
				tmpDirectory(tmpDirectory), marksCount(marksCount), error(error), errorMutex(errorMutex), round(round)
			{

			}

			void operator()()
			{
				uint64_t low = bound.first;
				uint64_t high = bound.second;
				ConcurrentBitVector candidateMask(Task::TASK_SIZE);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
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
							candidateMask.Reset();
							VertexRollingHash hash(hashFunction, task.str.begin() + 1, hashFunction.HashFunctionsNumber());
							size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
								if (Within(hash.GetVertexHash(), low, high) && definiteCount == vertexLength)
								{
									size_t inCount = DnaChar::IsDefinite(posPrev) ? 0 : 2;
									size_t outCount = DnaChar::IsDefinite(posExtend) ? 0 : 2;
									for (int i = 0; i < DnaChar::LITERAL.size() && inCount < 2 && outCount < 2; i++)
									{
										char nextCh = DnaChar::LITERAL[i];
										if (nextCh == posPrev || IsIngoingEdgeInBloomFilter(bitVector, hash, nextCh))
										{
											++inCount;
										}

										if (nextCh == posExtend || IsOutgoingEdgeInBloomFilter(bitVector, hash, nextCh))
										{
											++outCount;
										}
									}

									if (inCount > 1 || outCount > 1)
									{
										++marksCount;
										candidateMask.SetBitConcurrently(pos);
									}
								}

								if (pos + edgeLength < task.str.size())
								{
									char posPrev = task.str[pos];
									definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									hash.Update(posPrev, posExtend);
									assert(hash.Assert(task.str.begin() + pos + 1));
								}
								else
								{
									break;
								}
							}

							try
							{
								candidateMask.WriteToFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start, round));
							}
							catch (std::runtime_error & err)
							{
								ReportError(errorMutex, error, err.what());
							}
						}
					}
				}
			}

		private:
			std::pair<uint64_t, uint64_t> bound;
			const VertexRollingHashSeed & hashFunction;
			const ConcurrentBitVector & bitVector;
			size_t vertexLength;
			TaskQueue & taskQueue;
			const std::string & tmpDirectory;
			std::atomic<uint64_t> & marksCount;
			size_t round;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};



		class CandidateFinalFilteringWorker
		{
		public:
			CandidateFinalFilteringWorker(const VertexRollingHashSeed & hashFunction,
				size_t vertexLength,
				TaskQueue & taskQueue,
				OccurenceSet & occurenceSet,
				tbb::spin_rw_mutex & mutex,
				const std::string & tmpDirectory,
				size_t round,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : hashFunction(hashFunction), vertexLength(vertexLength), taskQueue(taskQueue), occurenceSet(occurenceSet),
				mutex(mutex), tmpDirectory(tmpDirectory), round(round), error(error), errorMutex(errorMutex)
			{

			}

			void operator()()
			{
				ConcurrentBitVector candidateMask(Task::TASK_SIZE);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
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
							VertexRollingHash hash(hashFunction, task.str.begin() + 1, 1);
							{
								try
								{
									candidateMask.ReadFromFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start, round), false);
								}
								catch (std::runtime_error & err)
								{
									ReportError(errorMutex, error, err.what());
								}
							}

							for (size_t pos = 1;; ++pos)
							{
								char posPrev = task.str[pos - 1];
								char posExtend = task.str[pos + vertexLength];
								if (candidateMask.GetBit(pos))
								{
									Occurence now;
									bool isBifurcation = false;						
									now.Set(hash.RawPositiveHash(0),
										hash.RawNegativeHash(0),
										task.str.begin() + pos,
										vertexLength,
										posExtend,
										posPrev,
										isBifurcation);
									size_t inUnknownCount = now.Prev() == 'N' ? 1 : 0;
									size_t outUnknownCount = now.Next() == 'N' ? 1 : 0;
									auto ret = occurenceSet.insert(now);
									typename OccurenceSet::iterator it = ret.first;
									if (!ret.second && !it->IsBifurcation())
									{
										inUnknownCount += DnaChar::IsDefinite(it->Prev()) ? 0 : 1;
										outUnknownCount += DnaChar::IsDefinite(it->Next()) ? 0 : 1;
										if (isBifurcation || it->Next() != now.Next() || it->Prev() != now.Prev() || inUnknownCount > 1 || outUnknownCount > 1)
										{
											it->MakeBifurcation();
										}
									}
								}

								if (pos + edgeLength < task.str.size())
								{
									char posPrev = task.str[pos];
									hash.Update(posPrev, posExtend);
									assert(hash.Assert(task.str.begin() + pos + 1));
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

		private:
			const VertexRollingHashSeed & hashFunction;
			size_t vertexLength;
			TaskQueue & taskQueue;
			OccurenceSet & occurenceSet;
			tbb::spin_rw_mutex & mutex;
			const std::string & tmpDirectory;
			size_t round;
			std::unique_ptr<std::runtime_error> & error;
			tbb::mutex & errorMutex;
		};

		struct EdgeResult
		{
			uint32_t pieceId;
			std::vector<JunctionPosition> junction;
		};

		static bool FlushEdgeResults(std::deque<EdgeResult> & result,
			JunctionPositionWriter & writer,
			std::atomic<uint64_t> & currentPiece)
		{
			if (result.size() > 0 && result.front().pieceId == currentPiece)
			{
				for (auto junction : result.front().junction)
				{
					writer.WriteJunction(junction);
				}

				++currentPiece;
				result.pop_front();
				return true;
			}

			return false;
		}

		class EdgeConstructionWorker
		{
		public:
			EdgeConstructionWorker(size_t vertexLength,
				TaskQueue & taskQueue,
				const BifurcationStorage<CAPACITY> & bifStorage,
				JunctionPositionWriter & writer,
				std::atomic<uint64_t> & currentPiece,
				std::atomic<uint64_t> & occurences,
				uint64_t & currentStubVertexId,
				tbb::mutex & currentStubVertexMutex,
				const std::string & tmpDirectory,
				size_t totalRounds,
				std::unique_ptr<std::runtime_error> & error,
				tbb::mutex & errorMutex) : vertexLength(vertexLength), taskQueue(taskQueue), bifStorage(bifStorage),
				writer(writer), currentPiece(currentPiece), occurences(occurences), tmpDirectory(tmpDirectory),
				error(error), errorMutex(errorMutex), currentStubVertexId(currentStubVertexId), currentStubVertexMutex(currentStubVertexMutex), totalRounds(totalRounds)
			{

			}							

			void operator()()
			{
				try
				{
					DnaString bitBuf;
					std::deque<EdgeResult> result;
					ConcurrentBitVector candidateMask(Task::TASK_SIZE);
					while (true)
					{
						Task task;
						if (taskQueue.try_pop(task))
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
								try
								{
									ConcurrentBitVector temporaryMask(Task::TASK_SIZE);
									for (size_t i = 0; i < totalRounds; i++)
									{
										temporaryMask.ReadFromFile(CandidateMaskFileName(tmpDirectory, task.seqId, task.start, i), true);
										candidateMask.MergeOr(temporaryMask);
									}
								}
								catch (std::runtime_error & err)
								{
									ReportError(errorMutex, error, err.what());
								}

								size_t taskstart = task.start;

								EdgeResult currentResult;
								currentResult.pieceId = task.piece;
								size_t definiteCount = std::count_if(task.str.begin() + 1, task.str.begin() + vertexLength + 1, DnaChar::IsDefinite);
								for (size_t pos = 1;; ++pos)
								{
									while (result.size() > 0 && FlushEdgeResults(result, writer, currentPiece));
									int64_t bifId(INVALID_VERTEX);
									assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
									if (definiteCount == vertexLength && candidateMask.GetBit(pos))
									{
										bifId = bifStorage.GetId(task.str.begin() + pos);
										if (bifId != INVALID_VERTEX)
										{
											occurences++;
											currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, bifId));
										}
									}

									if (((task.start == 0 && pos == 1) || (task.isFinal && pos == task.str.size() - vertexLength - 1)) && bifId == INVALID_VERTEX)
									{
										occurences++;
										currentStubVertexMutex.lock();								
										currentResult.junction.push_back(JunctionPosition(task.seqId, task.start + pos - 1, currentStubVertexId++));
										currentStubVertexMutex.unlock();
									}

									if (pos + edgeLength < task.str.size())
									{
										definiteCount += (DnaChar::IsDefinite(task.str[pos + vertexLength]) ? 1 : 0) - (DnaChar::IsDefinite(task.str[pos]) ? 1 : 0);
									}
									else
									{
										break;
									}
								}

								result.push_back(currentResult);
							}
						}
					}

					while (result.size() > 0)
					{
						FlushEdgeResults(result, writer, currentPiece);
					}
				}
				catch (std::runtime_error & e)
				{
					errorMutex.lock();
					error.reset(new std::runtime_error(e));
					errorMutex.unlock();
				}
			}

		private:
			size_t vertexLength;
			TaskQueue & taskQueue;
			uint64_t & currentStubVertexId;
			const BifurcationStorage<CAPACITY> & bifStorage;
			JunctionPositionWriter & writer;
			std::atomic<uint64_t> & currentPiece;
			std::atomic<uint64_t> & occurences;
			const std::string & tmpDirectory;
			std::unique_ptr<std::runtime_error> & error;
			size_t totalRounds;
			tbb::mutex & errorMutex;
			tbb::mutex & currentStubVertexMutex;
		};		
		
		class FilterFillerWorker
		{
		public:
			FilterFillerWorker(uint64_t low,
				uint64_t high,
				const VertexRollingHashSeed & hashFunction,
				ConcurrentBitVector & filter,
				size_t edgeLength,
				TaskQueue & taskQueue) : low(low), high(high), hashFunction(hashFunction), filter(filter), edgeLength(edgeLength), taskQueue(taskQueue)
			{

			}

			void operator()()
			{
				std::vector<uint64_t> setup;
				std::vector<uint64_t> hashValue;
				const char DUMMY_CHAR = DnaChar::LITERAL[0];
				const char REV_DUMMY_CHAR = DnaChar::ReverseChar(DUMMY_CHAR);
				while (true)
				{
					Task task;
					if (taskQueue.try_pop(task))
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
						size_t definiteCount = std::count_if(task.str.begin(), task.str.begin() + vertexLength, DnaChar::IsDefinite);
						VertexRollingHash hash(hashFunction, task.str.begin(), hashFunction.HashFunctionsNumber());						
						for (size_t pos = 0;; ++pos)
						{
							hashValue.clear();
							char prevCh = task.str[pos];
							char nextCh = task.str[pos + edgeLength - 1];
							assert(definiteCount == std::count_if(task.str.begin() + pos, task.str.begin() + pos + vertexLength, DnaChar::IsDefinite));
							if (definiteCount == vertexLength)
							{
								fistMinHash0 = hash.GetVertexHash();
								if (DnaChar::IsDefinite(nextCh))
								{
									hash.GetOutgoingEdgeHash(nextCh, hashValue);
								}
								else
								{
									hash.GetOutgoingEdgeHash(DUMMY_CHAR, hashValue);
									hash.GetOutgoingEdgeHash(REV_DUMMY_CHAR, hashValue);
								}

								if (pos > 0 && !DnaChar::IsDefinite(task.str[pos - 1]))
								{
									hash.GetIngoingEdgeHash(DUMMY_CHAR, hashValue);
									hash.GetIngoingEdgeHash(REV_DUMMY_CHAR, hashValue);
								}
							}

							hash.Update(prevCh, nextCh);
							assert(hash.Assert(task.str.begin() + pos + 1));
							if (definiteCount == vertexLength)
							{
								secondMinHash0 = hash.GetVertexHash();
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
						if (!filter.GetBit(hashValue))
						{
							filter.SetBitConcurrently(hashValue);
						}
					}

					setup.clear();
				}
			}

		private:
			uint64_t low;
			uint64_t high;
			const VertexRollingHashSeed & hashFunction;
			ConcurrentBitVector & filter;
			size_t edgeLength;
			TaskQueue & taskQueue;
		};


		static void DistributeTasks(const std::vector<std::string> & fileName,
			size_t overlapSize,
			std::vector<TaskQueuePtr> & taskQueue,
			std::unique_ptr<std::runtime_error> & error,
			tbb::mutex & errorMutex,
			std::ostream & logFile)
		{
			size_t record = 0;
			size_t nowQueue = 0;
			uint32_t pieceCount = 0;
#ifdef LOGGING
			logFile << "Starting a new stage" << std::endl;
#endif
			for (size_t file = 0; file < fileName.size(); file++)
			{
#ifdef LOGGING
				logFile << "Reading " << fileName[file] << std::endl;
#endif
				const std::string & nowFileName = fileName[file];
				for (StreamFastaParser parser(nowFileName); parser.ReadRecord(); record++)
				{
					{
						errorMutex.lock();
						if (error != 0)
						{
							throw *error;
						}

						errorMutex.unlock();
					}

					std::stringstream ss;
#ifdef LOGGING
					logFile << "Processing sequence " << parser.GetCurrentHeader() << " " << ss.str() << std::endl;
#endif
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
							for (bool found = false; !found; nowQueue = nowQueue + 1 < taskQueue.size() ? nowQueue + 1 : 0)
							{
								TaskQueuePtr & q = taskQueue[nowQueue];
								if (q->capacity() - q->size() > 0)
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
#ifdef LOGGING
									logFile << "Passed chunk " << prev << " to worker " << nowQueue << std::endl;
#endif
									prev = start - overlapSize + 1;
									buf.swap(overlap);
									found = true;
								}
							}

						}

					} while (!over);
				}
			}

			for (size_t i = 0; i < taskQueue.size(); i++)
			{
				TaskQueuePtr & q = taskQueue[nowQueue];
				while (!taskQueue[i]->try_push(Task(0, Task::GAME_OVER, 0, true, std::string())))
				{
					
				}
			}
		}

		uint64_t TrueBifurcations(const OccurenceSet & occurenceSet, std::ofstream & out, size_t vertexSize, size_t & falsePositives) const
		{
			uint64_t truePositives = falsePositives = 0;
			for (auto it = occurenceSet.begin(); it != occurenceSet.end();++it)
			{
				bool bifurcation = it->IsBifurcation();
				if (bifurcation)
				{
					++truePositives;
					it->GetBase().WriteToFile(out);
					if (!out)
					{
						throw StreamFastaParser::Exception("Can't write to a temporary file");
					}
				}
				else
				{
					falsePositives++;
				}
			}

			return truePositives;
		}

		size_t vertexSize_;
	};
}

#endif	
