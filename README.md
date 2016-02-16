

This is the code that implements the algorithm described in the paper
"An efficient  algorithm to build the compacted de Bruijn graph from many complete genomes"

Disclaimer: this is still a research prototype, the code has not been "officially" released yet.
Things like compilation, installation, output file format, and command line parameters
are subject to change. Though the correctness of the implementation was thoroughly tested,
the testing is still in progress. However, you can use the code to reproduce the experiments
from the paper.

To compile the code, you need the following (Linux only):

* CMake 
* A GCC compiler supporting C++11
* Intel TBB library properly installed on your system. In other words, G++
  should be able to find TBB libs 

Once you've got all the things above, do the following:

* Go to the root directory of the project and create the "build" folder
* Go to the "build" directory
* Type cmake ../src -DCMAKE_BUILD_TYPE=RELEASE
* Type make

This will make two targets: twopaco and graphdump

twopaco
----------------

To run the graph construction (assuming you're in the "build/graphconstructor" dir), type:

	./twopaco -q <number_of_hash_functions> -f <filter_size> -k <value_of_k> --tmpdir <directory_for_temporary_files> -o <output_file> -r <number_of_rounds> <input_files>

You can also type "./twopaco --help" to get parameter description.
Note that the size of the Bloom filter (in bits) is actualy "2^filter_size".
The output file is a binary file that indicates junction positions on the positive strand.
The file consists of pairs of numbers, where the first number is 4 bytes long, and the second is 8 bytes long.
The first number indicates position and the second one indicates the ID of the junction.
Positions appear in the file in the same order they appear in the input genomes
To make the output human readable, use the "graphdump" utility.

A note: this is a prototype version of the graph constructor. This output format is
mostly for testing purposes and not final. The release version will likely use the
GFA format, see: https://github.com/pmelsted/GFA-spec/issues/7

graphdump
---------

This utility turns the binary file into human readable one. Just run

	graphdump <input_file>

It will output a text file to the standard output. Each line will contain a 
triple indicating an occurence of junction:

	<seq_id_i> <pos_i> <junction_id_i>

The first number is the index number of the sequence, the second one is the position,
and the last one is the junction id. The index number of the sequence is just the order
of the sequence in the input file. All positions/orders count from 0.
Positions appear in the file in the same order they appear in the input genomes
This way, one can obtain all multi-edges of the graph with a linear scan, as described in the paper.
The negative strand can be obtained easily since it is symmetric.

One can also output junctions grouped by ids, it is useful for comparison between
different graphs:

	graphdump -g <input_file>

In this format the i-th line line corresponds to the i-th junction and is of format:

	<seq_id_0> <pos_0>; <seq_id_1> <pos_1>; ....

Where each pair "seq_id_i pos_j" corresponds to an occurence of the junction in
sequence "seq_id_i" at position "pos_j". Sequence ids are just the numbers of sequences
in the order they appear in the input. All numbers count from 0.

