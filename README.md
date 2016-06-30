TwoPaCo 0.0.0	

Release date: TBD
=================

Authors
=======
* Ilia Minkin (Pennsylvania State University)
* Son Pham (Salk Institute for Biological Studies)
* Paul Medvedev (Pennsylvania State University)

Introduction
============
It is an implementation of the algorithm described in the paper
"TwoPaCo: An efficient algorithm to build the compacted de Bruijn graph from
many complete genomes".

This distribution contains two programs:

* twopaco -- a tool for direct construction of the compressed graph from 
multiple complete genomes
* graphdump -- a utility that turns output of twopaco into text format

Disclaimer: this is still a research prototype and the code has not been
"officially" released yet. Things like compilation, installation, output
file format, and commandline parameters are subject to change.

Test data
=========
Links to the data used for bencharmking in the paper: https://github.com/medvedevgroup/TwoPaCo/blob/master/data.txt

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
Compilation under other platforms is possible, portable makefiles are in progress.

TwoPaCo usage
=============
To construct the graph (assuming you're in the "build/graphconstructor" dir), type:

	./twopaco -f <filter_size> -k <value_of_k> <input_files>

This will constuct the compressed graph for the vertex size of \<value_of_k\> using
2^\<filter_size\> bits in the Bloom filter. The output file is a binary file that
indicates junction positions on the positive strand. The default output file name
is "de_bruiijn. bin". You can it read directly using an API (will be documented later)
or make it into a text file using the "graphdump" utility. The program has several
additional parameters, see subsections below. You can also type "./twopaco --help"
to get a short parameter description. Note that TwoPaCo supports only odd values 
of k for the sake of easier graph manipulation.

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
This utility turns the binary file a text one. There are several output formats
available.

GFA
---
GFA is the most handy option. It **explicitly** represents the graph as a list
non-branching paths, adjacencies between them, and all occurrences of the segments
in the input genomes, i.e. a colored de Bruijn graph. Work on the format specification
is still in progress, the version implemented in TwoPaCo can be found here:

	https://github.com/ilyaminkin/GFA-spec/blob/record_style_occurence/GFA-spec.md

To get GFA output, run:

	graphdummp --gfa -k <value_of_k> -s <input_genomes> <twopaco_output_file>

Although GFA is quite easy to understand and parse, it is quite bulky, especially
for larger graphs.

Junctions List Format
---------------------
In this format the output file only contains positions of junctions in the input
genomes. As described in the paper, you can trivially restore information about
edges from this junctions list. Note that junctions are mapped to genomes, i.e.
one can reconstruct a **colored graph** from it. To get the junctions list, run:

	graphdump <twopaco_output_file>

This command will output a text file to the standard output. Each line will contain a 
triple indicating an occurence of junction:

	<seq_id_i> <pos_i> <junction_id_i> 

The first number is the index number of the sequence, the second one is the
position, and the third one is the junction id. The index number of the sequence
is the order of the sequence in the input file(s). All positions/orders count
from 0. Positions appear in the file in the same order they appear in the input
genomes. The \<junction_id\> is a signed integer, the id of the junction that
appears on the positive strand strand. A positive number indicates "direct" version
of the junction, while a negative one shows the reverse complimentary version of the
same junction. For example +1 and -1 are different versions of the same junction.
This way, one can obtain all multi-edges of the graph with a linear scan, as described
in the paper. For example, a sequence of pairs of junctions ids:

	a_1
	a_2
	a_3

Generates edges a_1 -> a_2, a_2 -> a_3 in the graph corresponding to the positive
strand. To obtain the edges from the positive strand, one has to traverse them 
in the backwards order and negate signs, e.g. for the example above the sequence
will be -a_3 -> -a_2 -> -a_1. One can also output junctions grouped by ids, it is
useful for comparison between different graphs:

	graphdump -g <twopaco_output_file>

In this format the i-th line line corresponds to the i-th junction and is of format:

	<seq_id_0> <pos_0>; <seq_id_1> <pos_1>; ....

Where each pair "seq_id_i pos_j" corresponds to an occurence of the junction in
sequence "seq_id_i" at position "pos_j". Sequence ids are just the numbers of sequences
in the order they appear in the input. All positions count from 0.

Read The Binary File Directly
-----------------------------
This is the most parsimonious option in terms of involved resources.
One can read junctions and/or edges from the output file using a very simple
C++ API. It will be documented slightly later.

License
=======
See LICENSE.txt

Contacts
========
Please e-mail your feedback at ivminkin@gmail.com.

You also can report bugs or suggest features using issue tracker at GitHub
https://github.com/medvedevgroup/TwoPaCo
