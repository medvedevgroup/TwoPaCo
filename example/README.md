Introduction
------------
Here is an example of the de Bruijn graph built from a small input with k=11.
Below are commands used to produce the files. Run TwoPaCo:
	
	twopaco -f 20 -k 11 example.fa -o example.dbg

Get DOT file for Graphviz visualisation and render it:

	graphdump -f dot example.dbg -k 11 > example.dot
	dot -Tpng example2.dot > example2.png

Note that the image shows vertices as IDs of the junctions in the graph. To
get the sequences of junctions one may need the list of the junctions in the
order they appear in the input:

	graphdump -f seq example.dbg -k 11 > example.seq

Note that this list is only for the direct strand. To get the list of the junctions
on the reverse strand, reverse the order of the junctions and negate signs. The 
future release will support GFF file for export of those coordinates. To get a
GFA1 file:
	
	graphdump -f gfa1 -k 11 example.dbg -s example.fa > example_gfa1.gfa

or GFA2:

	graphdump -f gfa2 -k 11 example.dbg -s example.fa > example_gfa2.gfa

The resulting GFA1 file can be visualized by Bandage for example, see "example_bandage.png".
One can simply get sequences of all compressed paths in FASTA format:

	graphdump -f fasta -k 11 example.dbg -s example.fa > example_paths.fa