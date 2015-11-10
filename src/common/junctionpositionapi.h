#ifndef _JUNCTION_POSITION_API_H_
#define _JUNCTION_POSITION_API_H_

#include <fstream>
#include <exception>
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
		friend class JunctionPositionWriter;
	};

	class JunctionPositionReader
	{
	public:
		JunctionPositionReader(const std::string & inFileName)
		{

		}

	private:
		std::ifstream in_;
	};
	

	class JunctionPositionWriter
	{
	public:
		JunctionPositionWriter(const std::string & outFileName) : out_(outFileName.c_str(), std::ios::binary)
		{			
			if (!out_)
			{
				throw std::runtime_error("Can't create the output file");
			}
		}

		void WriteSeparator()
		{
			JunctionPosition junction(-1, -1);
			WriteJunction(junction);
		}

		void WriteJunction(JunctionPosition pos)
		{
			out_.write(reinterpret_cast<const char*>(&pos.pos_), sizeof(pos.pos_));
			out_.write(reinterpret_cast<const char*>(&pos.bifId_), sizeof(pos.bifId_));
		}

	private:
		std::ofstream out_;
	};
}

#endif