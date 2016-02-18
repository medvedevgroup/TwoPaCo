TwoPaCo 0.0.0	

Release date: TBD
=================

Authors
=======
* Ilia Minkin (Pennsylvania State University)
* Son Pham
* Paul Medvedev (Pennsylvania State University)

Introduction
============
It is an implementation of the algorithm described in the paper
"TwoPaCo: An efficient algorithm to build the compacted de Bruijn graph from
many complete genomes".

Test data can be downloaded via the following links:

* 62 E.Coli: https://s3-us-west-2.amazonaws.com/graph.testdata/ecoli.tar.gz

* Seven human: https://s3-us-west-2.amazonaws.com/graph.testdata/human.tar

* Eight primates: https://s3-us-west-2.amazonaws.com/graph.testdata/human.tar

* 93 simulated humans: https://s3-us-west-2.amazonaws.com/graph.testdata/human_sim.tar.gz

This distribution contains two programs:

* twopaco -- a tool for direct construction of the compressed graph from 
multiple complete genomes

* graphdump -- a utility that turns output of twopaco into text format

Disclaimer: this is still a research prototype and the code has not been
"officially" released yet. Things like compilation, installation, output
file format, and commandline parameters are subject to change.

Compilation
===========
To compile the code, you need the following (Linux only):

* CMake 
* A GCC compiler supporting C++11
* Intel TBB library properly installed on your system. In other words, G++
  should be able to find TBB libs 

Once you've got all the things above, do the following:

* Go to the root directory of the project and create the "build" folder
* Go to the "build" directory
* Run cmake ../src
* Run make

This will make two targets: twopaco and graphdump.
Compilation under platforms is possible, portable makefiles are in progress.

The twopaco usage
=============
To construct the graph (assuming you're in the "build/graphconstructor" dir), type:

	./twopaco -f <filter_size> -k <value_of_k> <input_files>

This will constuct the compressed graph for the vertex size of \<value_of_k\> using
2^\<filter_size\> bits in the Bloom filter. The output file is a binary file that
indicates junction positions on the positive strand. The default output file name
is "de_bruiijn. bin". You can read directly using an API (will be documented later)
or make it into a text file using the "graphdump" utility. The program has several
additional parameters, see subsections below. You can also type "./twopaco --help"
to get a short parameter description.

A note: the release version will likely use the GFA format, see:
https://github.com/pmelsted/GFA-spec/issues/7

Here is description of other parameters 

Number of hash functions
------------------------
The number of hash functions used for the Bloom filter. The default is 5. To
change, use:

	-q <number> or --hashfnumber <number>

Number of rounds
----------------
Number of computational rounds. Each rounds processses a separate subset of k-mers
which reduces memory usage. The default is 1. To change, use:

	-r <number> or --rounds <number>

Number of threads
-----------------
twopaco can be run in multiple threads. The default is 1. To change, use:

	-t <number> or --threads <number>

Temporary directory
-------------------
The directory for temporary files. The default is the current working directory.
To change, use (the directory must exist!):

	--tmpdir <path_to_the_directory>

Output file name
----------------
The name of the output file. The default is "de_bruijn.bin". To change, use:

	--o <file_name> or --outfile <file_name>

The graphdump usage
===================
This utility turns the binary file a text one one. You can run:

	graphdump <twopaco_output_file>

It will output a text file to the standard output. Each line will contain a 
triple indicating an occurence of junction:

	<seq_id_i> <pos_i> <positive_junction_id_i> <negative_junction_id_i>

The first number is the index number of the sequence, the second one is the
position, and the last two are the junction id. The index number of the sequence
is just the order of the sequence in the input file. All positions/orders count
from 0. Positions appear in the file in the same order they appear in the input
genomes. The \<positive_junction_id\> is the id of the junction that appears on
the direct strand, while \<negative_junction_id_i> is the id of the junction
that appears on the complementary strand. This way, one can obtain all multi-edges
of the graph with a linear scan, as described in the paper. For example, a sequence
of pairs of junctions ids:

	a_1 b_1
	a_2 b_2
	a_3 b_3

Generates edges a_1 -> a_2, a_2 -> a_3 in the graph corresponding to the positive
strand and edges b_3 -> b_2, b_2 -> b_1 in the graph coming from the reverse 
complementary strand.	

One can also output junctions grouped by ids, it is useful for comparison between
different graphs:

	graphdump -g <twopaco_output_file>

In this format the i-th line line corresponds to the i-th junction and is of format:

	<seq_id_0> <pos_0>; <seq_id_1> <pos_1>; ....

Where each pair "seq_id_i pos_j" corresponds to an occurence of the junction in
sequence "seq_id_i" at position "pos_j". Sequence ids are just the numbers of sequences
in the order they appear in the input. All numbers count from 0.


License
=======
See LICENSE.txt

Contacts
========
Please e-mail your feedback at ivminkin@gmail.com.

You also can report bugs or suggest features using issue tracker at GitHub
https://github.com/medvedevgroup/TwoPaCo
