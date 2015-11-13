#include <junctionpositionapi.h>

int main(int argc, char * argv[])
{
	Sibelia::JunctionPosition pos;
	Sibelia::JunctionPositionReader reader("de_bruijn.bin");
	while (reader.NextJunctionPosition(pos))
	{
		std::cout << pos.GetChr() << ' ' << pos.GetPos() << ' ' << pos.GetId() << std::endl;
	}

	return 0;
}