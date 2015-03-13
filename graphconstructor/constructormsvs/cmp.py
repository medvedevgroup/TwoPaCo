def to_line(fn):
	line = set([x for x in open(fn).next().split()])
	return line

ln1 = to_line("hd1.txt")
ln2 = to_line("hd4.txt")

print ln1 - ln2
print ln2 - ln1