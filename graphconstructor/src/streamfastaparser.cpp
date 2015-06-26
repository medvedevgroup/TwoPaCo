#include <sstream>
#include <algorithm>

#include "streamfastaparser.h"

namespace Sibelia
{
	//const std::string StreamFastaParser::VALID_CHARS = "ACGTURYKMSWBDHWNX-";
	const std::string StreamFastaParser::VALID_CHARS = "ACGT";

	char StreamFastaParser::MakeUp(char ch)
	{
		switch (ch)
		{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			return UNKNOWN_CHAR;
		}
	}

	char StreamFastaParser::Reverse(char ch)
	{
		switch (ch)
		{
		case 0:
			return 3;
		case 3:
			return 0;
		case 1:
			return 2;
		case 2:
			return 1;
		}

		return UNKNOWN_CHAR;
	}

	char StreamFastaParser::UnMakeUp(char ch)
	{
		switch (ch)
		{
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		default:
			return 'N';
		}
	}

	StreamFastaParser::Exception::Exception(const std::string & msg) : std::runtime_error(msg)
	{
		
	}

	StreamFastaParser::~StreamFastaParser()
	{
		delete [] buffer_;
	}

	StreamFastaParser::StreamFastaParser(const std::string & fileName) : stream_(fileName.c_str()),
		buffer_(new char[BUF_SIZE]), bufferPos_(0), bufferSize_(0)
	{
		if (!stream_ && !stream_.eof())
		{
			throw Exception("Can't open file " + fileName);
		}

		isValid_.assign(256, 0);
		for (char ch : VALID_CHARS)
		{
			isValid_[ch] = 1;
		}
	}

	bool StreamFastaParser::ReadRecord()
	{
		char ch = '\0';
		if (GetCh(ch))
		{
			if (ch != '>')
			{
				throw Exception("The FASTA header should start with a '>', started with '" + std::string(1, ch) + "'");
			}
		}
		else
		{
			return false;
		}

		std::stringstream ss;
		while (GetCh(ch))
		{
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

	bool StreamFastaParser::GetChar(char & ch, bool makeUp)
	{
		while (true)
		{
			if (!Peek(ch))
			{
				break;
			}
			
			if (isspace(ch))
			{
				GetCh(ch);
				continue;
			}
			else if (ch == '>')
			{
				return false;
			}
			else
			{
				if (!isValid_[toupper(ch)])
				{
					throw Exception("Found an invalid character '" + std::string(1, ch) + "'");
				}

				GetCh(ch);
				return true;
			}
		}

		return false;
	}

	bool StreamFastaParser::GetCh(char & ch)
	{
		if (bufferPos_ == bufferSize_)
		{
			if (stream_)
			{				
				stream_.read(buffer_, BUF_SIZE);
				bufferPos_ = 0;
				bufferSize_ = stream_.gcount();
			}
			else
			{
				return false;
			}
		}

		ch = buffer_[bufferPos_++];
		return true;
	}

	bool StreamFastaParser::Peek(char & ch)
	{
		if (bufferPos_ == bufferSize_)
		{
			if (stream_)
			{
				stream_.read(buffer_, BUF_SIZE);
				bufferPos_ = 0;
				bufferSize_ = stream_.gcount();
			}
			else
			{
				return false;
			}
		}

		ch = buffer_[bufferPos_];
		return true;
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
