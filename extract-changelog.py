import sys

ver = sys.argv[1][1:]

line = sys.stdin.readline()
while not line.startswith("Version {}".format(ver)):
    line = sys.stdin.readline()

line = sys.stdin.readline()
line = sys.stdin.readline()
line = sys.stdin.readline()

while not line.startswith("Version"):
    print(line, end="")
    line = sys.stdin.readline()
