#include <sstream>
#include <algorithm>

#include "streamfastaparser.h"

namespace Sibelia
{
	const std::string StreamFastaParser::VALID_CHARS = "ACGTURYKMSWBDHWNX-";

	StreamFastaParser::Exception::Exception(const std::string & msg): std::runtime_error(msg)
	{
	}

	StreamFastaParser::StreamFastaParser(const std::string & fileName): stream_(fileName.c_str())
	{
		if (!stream_ && !stream_.eof())
		{
			throw Exception("Can't open file " + fileName);
		}
	}

	bool StreamFastaParser::ReadRecord()
	{
		if (stream_.eof())
		{
			return false;
		}

		char ch;
		stream_.get(ch);
		if (ch != '>')
		{
			throw Exception("The FASTA header should start with a '>', started with '" + std::string(1, ch) + "'");
		}			

		std::stringstream ss;
		while (stream_)
		{
			stream_.get(ch);			
			if (ch == '\n')
			{
				ss >> currentHeader_;
				break;
			}
			else
			{
				ss << ch;
			}
		}

		return true;
	}

	bool StreamFastaParser::GetChar(char & ch)
	{
		while (stream_)
		{			
			ch = stream_.peek();
			if (stream_.eof())
			{
				break;
			}
			else if (isspace(ch))
			{
				stream_.get(ch);
				continue;
			}
			else if (ch == '>')
			{
				return false;
			}
			else
			{
				if (std::find(VALID_CHARS.begin(), VALID_CHARS.end(), toupper(ch)) == VALID_CHARS.end())
				{
					throw Exception("Found an invalid character '" + std::string(1, ch) + "'");
				}

				stream_.get(ch);
				return true;
			}
		}

		return false;
	}

	std::string StreamFastaParser::GetCurrentHeader() const
	{
		return currentHeader_;
	}

	std::string StreamFastaParser::GetErrorMessage() const
	{
		return errorMessage_;
	}
}
