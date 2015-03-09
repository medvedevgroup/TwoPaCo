import sys
buf = ""
for line in open(sys.argv[1]):
	line = line.strip()
	if line and line[0] == '>':
		print line
		print buf
		buf = ""
	else:
		buf = buf + line

print buf	