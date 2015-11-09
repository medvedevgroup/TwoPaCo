#ifndef _JUNCTION_POSITION_API_H_
#define _JUNCTION_POSITION_API_H_

#include <boost/thread.hpp>

namespace Sibelia
{
	struct JunctionPosition
	{
	public:
		JunctionPosition() {}
		JunctionPosition(uint32_t pos, uint64_t bifId) : pos_(pos), bifId_(bifId) {}
		uint32_t GetPos() const
		{
			return pos_;
		}

		uint64_t GetId() const
		{
			return bifId_;
		}

	private:
		uint32_t pos_;
		uint64_t bifId_;		
	};

	class JunctionPositionReader
	{
	public:
	private:
	};
	

	class JunctionPositionWriter
	{
	public:
		JunctionPositionWriter(const std::string & outFileName) : out_(outFileName.c_str(), std::ios::binary)
		{			
			if (!out_)
			{
				throw StreamFastaParser::Exception("Can't create the output file");
			}
		}

		void WriteConcurrently(JunctionPosition pos)
		{
			boost::lock_guard guard(mutex_);

		}

	private:
		std::ostream out_;
		boost::mutex mutex_;
	};
}

#endif