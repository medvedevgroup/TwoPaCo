TwoPaCo 0.9.2

Release date: 18th May 2017
=============================

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
* graphdump -- a utility that turns output of twopaco into a text format

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
To construct the graph (assuming you are in dir with "twopaco"), type:

	./twopaco -f <filter_size> -k <value_of_k> <input_files>

This will constuct the compressed graph for the vertex size of \<value_of_k\> using
2^\<filter_size\> bits in the Bloom filter. The output file is a binary that can be
either converted to a text file or read directly using an API (will be available soon).

The filter size -f is a very important parameter that affects both
the memory usage and the speed. TwoPaCo will use at least 2^\<filter_size\> / 8 bytes of
memory, but setting it too low can massively increase the size of
the memory used and slow down the program. We recommend the user to set -f to to the value
so that 2^\<filter_size\> / 8 is the maximum memory in bytes they wish to allocate to
the algorithm. If the memory usage then exceeds the value above, then the number of rounds
should be increased until the memory usage falls below the desired value
(see the section ["Number of rounds"](https://github.com/medvedevgroup/TwoPaCo#number-of-rounds)).

If the memory usage is not a concern, then as a rule of thumb for the fastest speed,
set the parameter -f as large as possible. Here are the recommended settings given
the memory size of a machine:

1) 4GB - 34
2) 8GB - 35
3) 16GB - 36
4) 32GB - 37
5) 64GB - 38
6) 128GB - 39
7) 256GB - 40

For a memory size in between, go up a value, i.e. for 12GB RAM use 36, not 35.
For more details on how the Bloom filter size affects performance, please see the paper.
Below is description of the other parameters.

Number of rounds
----------------
Number of computational rounds. For the fastest performance, use 1 round (the default).
Increasing the number of rounds will decrease the memory usage at the expense of longer
runtime. To set the parameter, use:

	-r <number> or --rounds <number>

K-mer size
----------
This value sets the size of a vertex in the de Bruijn graph. Default is 25, to
change, use:

	-k <number> or --kvalue <number>

Note that:
1) TwoPaCo uses **k** as the size of the vertex and **(k + 1)** as the size of 
the edge 
2) **k** must be odd

The maximum value of K supported by TwoPaCo is determined at the compile time.
To increase the max value of K, increase the value "MAX_CAPACITY" defined in the
header "vertexenumerator.h" and recompile. The value of "MAX_CAPACITY" should be
at least (K + 4) / 32 + 1. Note that increasing the parameter will slow down 
the compilation.

Number of hash functions
------------------------
The number of hash functions used for the Bloom filter. The default is five. To
change, use:

	-q <number> or --hashfnumber <number>

More hash functions increases the running time. At the same time, more hash functions
may decrease the number of false positives and the memory usage.


Number of threads
-----------------
twopaco can be run with multiple threads. The default is 1. To change, use:

	-t <number> or --threads <number>

Temporary directory
-------------------
The directory for temporary files. The default is the current working directory.
To change, use (the directory must exist):

	--tmpdir <path_to_the_directory>

Output file name
----------------
The name of the output file. The default is "de_bruijn.bin". To change, use:

	--o <file_name> or --outfile <file_name>

Running tests
-------------
If the flag is set, TwoPaCo will run a set of internal tests instead of
processing the input file:

	--test

The graphdump usage
===================
This utility turns the binary file a text one. There are several output formats
available. The folder "example" contains an example described in details.

GFF
---
In the next release I will add an option to output coordinates of all occurrences
of the junctions in GFF format.

DOT
---
This format is used for visualization. The resulting DOT file can be converted
into an image using Graphviz package:

	http://www.graphviz.org/

To get the DOT file, use:

	graphdump <twopaco_output_file> -f dot -k <value_of_k>

Note that the graph is a union of graphs built from both strands, with blue edges
coming from the main strand and red ones from reverse one. The labels of the edges
will indicate its position on a chromosome.

GFA
---
GFA is the most handy option. It **explicitly** represents the graph as a list of
edges (non-branching paths in the non-compacted de Bruijn graph) graph and adjacencies
between them. The file also contains all occurrences of the strings spelled by the paths
in the input genomes.

In other words, it describes a colored de Bruijn graph where each path is mapped
to several locations in the input ("colored"). TwoPaCo supports both GFA1 and
GFA2. They are described here:

	https://github.com/GFA-spec/GFA-spec

To get GFA output, run:

	graphdummp <twopaco_output_file> -f gfa[version] -k <value_of_k> -s <input_genomes>

In the resulting file compacted non-branching paths are "segments" with "links"
(GFA1) or "edges" (GFA2) containing them. "Containment" (GFA1) or "Fragment"
(GFA2) records desrcibe the mapping between the non-branching paths in the
graph and the input genomes. For GFA1, each input chromosome is also a "segment"
described in the very beginning of the GFA file. 

GFA1 only: each segment representing an input chromosome has the name of the
corresponding header of the sequence in input FASTA file. In case if there are
duplicate headers, one can add a prefix to segment names:

	"s<number>_" + header of the sequence in input FASTA file

To do so, use the switch:

	--prefix

For an example of GFA output and more detailed explanation, see the "example"
folder.

Junctions List Format
---------------------
In this format the output file only contains positions of junctions in the input
genomes. As described in the paper, you can trivially restore information about
edges from this junctions list. Note that junctions are mapped to genomes, i.e.
one can reconstruct a **colored graph** from it. To get the junctions list, run:

	graphdump -f seq -k <value_of_k>

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
in the paper. For example, a sequence of of junctions ids:

	a_1
	a_2
	a_3

Generates edges a_1 -> a_2, a_2 -> a_3 in the graph corresponding to the positive
strand. To obtain the edges from the positive strand, one has to traverse them 
in the backwards order and negate signs, e.g. for the example above the sequence
will be -a_3 -> -a_2 -> -a_1. One can also output junctions grouped by ids, it is
useful for comparison between different graphs:

	graphdump -f group -k <value_of_k>

In this format the i-th line line corresponds to the i-th junction and is of format:

	<seq_id_0> <pos_0>; <seq_id_1> <pos_1>; ....

Where each pair "seq_id_i pos_j" corresponds to an occurence of the junction in
sequence "seq_id_i" at position "pos_j". Sequence ids are just the numbers of sequences
in the order they appear in the input. All positions count from 0.

Read The Binary File Directly
-----------------------------
This is the most parsimonious option in terms of involved resources.
One can read junctions and/or edges from the output file using a very simple
C++ API. I will add the description in the future release. For now, one can
use the sources of graphdump as a reference, it is relatively straightforward.

License
=======
See LICENSE.txt

Contacts
========
Please e-mail your feedback at ivminkin@gmail.com.

You also can report bugs or suggest features using issue tracker at GitHub
https://github.com/medvedevgroup/TwoPaCo

Citation
================
If you use TwoPaCo, please cite:

	Ilia Minkin, Son Pham, and Paul Medvedev
	"TwoPaCo: An efficient algorithm to build the compacted de Bruijn graph from many complete genomes"
	Bioinformatics, 2016 doi:10.1093/bioinformatics/btw609

This project has been supported in part by NSF awards DBI-1356529, CCF-1439057, IIS-1453527, and IIS-1421908.
