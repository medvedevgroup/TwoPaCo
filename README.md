This is the code that implements the algorithm described in the paper
"An efficient  algorithm to build the compacted de Bruijn graph from many complete genomes"

To compile the code, you need the following (Linux only):

* CMake 
* A GCC compiler supporting C++11
* Boost library installed on your system
* Intel TBB library installed on your system

Once you've got all the things above, do the following:

* Go to the "build" directory
* Type cmake ../src -DBOOST_ROOT=D1 -DTBB_INCLUDE_DIR=D2 -DTBB_LINK_DIR=D3
Where D1, D2 and D3 are the Boost root directory, the directory with TBB include files
and the directory with TBB library files.
* Type make

To run the graph construction (assuming you're in the "build" dir), type:

	./graphconstructor -q <number_of_hash_functions> -f <filter_size> -k <value_of_k> --tmpdir <directory_for_temporary_files> -o <output_file> -r <number_of_rounds>

You can also type ./graphconstructor --help to get parameter description.
Note that the size of the Bloom filter is actualy 2^{<filter_size>}.
The output file is a binary file that indicates junction positions on the positive strand.
The file consists of pairs of numbers, where the first number is 4 bytes long, and the second is 8 bytes long.
The first number indicates position and the second one indicates the ID of the junction.
The utilities for manipulating the output will be provided a bit later.
The negative strand can be obtained easily since it is symmetric.