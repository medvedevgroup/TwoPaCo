TwoPaCo 0.0.0	

Release date: TBD
=================

Authors
=======
* Ilia Minkin (Pennsylvania State University)
* Paul Medvedev (Pennsylvania State University)

Introduction
============
It is an implementation of the algorithm described in the paper
"TwoPaCo: An efficient algorithm to build the compacted de Bruijn graph from
many complete genomes".

This distribution contaisn two programs:

* twopaco -- a tool for direct construction of the compressed graph from 
multiple complete genomes

* graphdump -- a utility that turns output of twopaco into text format

Disclaimer: though the correctness of the implementation was thoroughly tested,
this is still a research prototype and the code has not been "officially"
released yet. Things like compilation, installation, output file format, and
commandline parameters are subject to change.

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
* Type cmake ../src -DCMAKE_BUILD_TYPE=RELEASE
* Type make

This will make two targets: twopaco and graphdump.

Compilation under platforms is possible, portable makefiles are in progress.

Usage
=====

To run the graph construction (assuming you're in the "build/graphconstructor" dir), type:

	./twopaco -f <filter_size> -k <value_of_k> <input_files>

This will constuct the compressed graph for the vertex size of \<value_of_k\> using
2^\<filter_size\> bits in the Bloom filter. The output file is a binary file that
indicates junction positions on the positive strand. You can read directly using
an API (will be documented later) or make it into a text file using the "graphdump"
utility. The program has several additional parameters, see subsections below.
You can also type "./twopaco --help" to get a short parameter description.

Note that the size of the Bloom filter (in bits) is actualy "2^filter_size".

The file consists of pairs of numbers, where the first number is 4 bytes long, and the second is 8 bytes long.
The first number indicates position and the second one indicates the ID of the junction.
Positions appear in the file in the same order they appear in the input genomes
To make the output human readable, use the "graphdump" utility.

A note: this is a prototype version of the graph constructor. This output format is
mostly for testing purposes and not final. The release version will likely use the
GFA format, see: https://github.com/pmelsted/GFA-spec/issues/7

License
=======

Contacts
========
Please e-mail your feedback at ivminkin@gmail.com.

You also can report bugs or suggest features using issue tracker at GitHub
https://github.com/medvedevgroup/TwoPaCo
