#ifndef _BIFURCATION_STORAGE_H_
#define _BIFURCATION_STORAGE_H_

#include "common.h"
#include "compressedstring.h"
#include "ngramhashing/cyclichash.h"
#include <fstream>
#include <cmph.h>

namespace Sibelia
{
	typedef CyclicHash<uint64_t> HashFunction;
	typedef std::unique_ptr<HashFunction> HashFunctionPtr;

	template<size_t CAPACITY>
	class BifurcationStorage
	{
	public:
		typedef CompressedString<CAPACITY> DnaString;
		BifurcationStorage(){}

		uint64_t GetDistinctVerticesCount() const
		{
			return bifurcationKey_.size();
		}

		uint64_t GetTotalVerticesCount() const
		{
			return bifurcationKey_.size() * 2;
		}

		void Init(std::istream & bifurcationTempRead, uint64_t verticesCount, uint64_t vertexLength, size_t threads)
		{
			std::ofstream outfile;
			outfile.open("temp.txt", std::ios_base::app);
	
			uint64_t bitsPower = 0;
			vertexLength_ = vertexLength;
			while (verticesCount * 8 >= (uint64_t(1) << bitsPower))
			{
				++bitsPower;
			}

			DnaString buf;
			for (size_t i = 0; i < verticesCount; i++)
			{
				buf.ReadFromFile(bifurcationTempRead);
				if (!bifurcationTempRead)
				{
					throw StreamFastaParser::Exception("Can't read from a temporary file");
				}

				bifurcationKey_.push_back(buf);
			}
			/*----------------------------------------------------------------*/
			for(size_t j = 0; j < verticesCount; j++)
			{
				outfile << bifurcationKey_.at(j).ToString(vertexLength);
				outfile << "\n";
			}
			
			FILE * keys_fd = fopen("temp.txt", "r");

  			if (keys_fd == NULL) 
  			{
  			  fprintf(stderr, "File \"temp.txt\" not found\n");
  			  exit(1);
  			}
	
  			// Source of keys
  			cmph_io_adapter_t *source = cmph_io_nlfile_adapter(keys_fd);
  			cmph_config_t *config = cmph_config_new(source);
  			cmph_config_set_algo(config, CMPH_BDZ);
  			new_cmph_hash = cmph_new(config);
  			cmph_config_destroy(config);
  			/*----------------------------------------------------------------*/
		}

		uint64_t GetId(std::string::const_iterator pos) const
		{
			DnaString bitBuf;
			bool posFound = false;
			bool negFound = false;
			uint64_t ret = INVALID_VERTEX;
			bitBuf.CopyFromString(pos, vertexLength_);
			
			
			/*----------------------------------------------------------------*/
			const char *pos_key = (bitBuf.ToString(vertexLength_)).c_str();
			bitBuf.Clear();
			bitBuf.CopyFromReverseString(pos, vertexLength_);
			const char *neg_key = (bitBuf.ToString(vertexLength_)).c_str();
			
			unsigned int pos_id = cmph_search(new_cmph_hash, pos_key, (cmph_uint32)strlen(pos_key));
			if(pos_id != 0)
			{
				posFound = true;
				ret = pos_id;
			}
			
			if(!posFound)
			{
				unsigned int neg_id = cmph_search(new_cmph_hash, neg_key, (cmph_uint32)strlen(pos_key));
				if(neg_id != 0)
				{
					negFound = true;
					ret = neg_id;
				}
			}
			
			if(negFound && !posFound && !DnaChar::IsSelfReverseCompliment(pos, vertexLength_))
			{
				ret = ret + bifurcationKey_.size();;
			}
			
			
			/*----------------------------------------------------------------
			
			auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
			if (it != bifurcationKey_.end() && *it == bitBuf)
			{
				posFound = true;
				ret = it - bifurcationKey_.begin();
			}

			if (!posFound)
			{
				bitBuf.Clear();
				bitBuf.CopyFromReverseString(pos, vertexLength_);
				auto it = std::lower_bound(bifurcationKey_.begin(), bifurcationKey_.end(), bitBuf, DnaString::Less);
				if (it != bifurcationKey_.end() && *it == bitBuf)
				{
					negFound = true;
					ret = it - bifurcationKey_.begin();
				}
			}

			if (negFound && !posFound && !DnaChar::IsSelfReverseCompliment(pos, vertexLength_))
			{
				ret += bifurcationKey_.size();
			}*/

			return ret;
		}

	private:
		size_t vertexLength_;
		std::vector<DnaString> bifurcationKey_;
		cmph_t *new_cmph_hash;

	};
}

#endif
