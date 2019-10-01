import sys

sample_file = sys.argv[1]
t_count_file = sys.argv[2]
seq_len_file = sys.argv[3]

t_count = 0
seq_lengths = []
for line in open(sample_file):
    name, seq = line.strip().split()
    t_count += seq.count('t')
    seq_lengths.append(len(seq))

print(t_count, file=open(t_count_file, 'w'))

file_obj = open(seq_len_file, 'w')
for seq_length in seq_lengths:
    print(seq_length, file=file_obj)
