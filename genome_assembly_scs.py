import itertools

def readFastq(filename):
    seqs = []

    with open(filename) as f:
        lines = f.readlines()

        for i in range(0, len(lines), 4):
            seqs.append(lines[i+1].rstrip())
            if(len(lines[i+1].rstrip()) != 100):
                print(len(lines[i+1].rstrip()))

    return seqs

#calculates the overlap with min_len between two reads
def overlap(x,y,min_len):
    start = 0
    while True:
        start = x.find(y[:min_len], start)
        if start == -1:
            return 0
        if y.startswith(x[start:]):
            return len(x) - start
        start += 1

#calculates the maximum overlap between two sequences in the input array reads
def maximal_overlap(reads, k): 
    reada, readb = None, None 
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a,b, min_len = k)
        if olen > best_olen:
            reada, readb = a,b
            best_olen = olen
    return reada, readb, best_olen

#implements the greedy shortest common superstring algorithm
def greedy_scs (reads, k):
    read_a, read_b, olen = maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = maximal_overlap(reads, k)
    return ''.join(reads)    

reads = readFastq('mystery_virus.fq')
assembled_ss = greedy_scs(reads, k=30)

with open('mystery_virus_genome.txt', 'w') as f:
    f.write(assembled_ss)

num_A = assembled_ss.count('A')
num_C = assembled_ss.count('C')
num_G = assembled_ss.count('G')
num_T = assembled_ss.count('T')

with open('output_stats.txt', 'w') as f1:
    f1.write("Number of A's = ")
    f1.write(str(num_A))
    f1.write('\n')
    f1.write("Number of C's = ")
    f1.write(str(num_C))
    f1.write('\n')
    f1.write("Number of G's = ")
    f1.write(str(num_G))
    f1.write('\n')
    f1.write("Number of T's = ")
    f1.write(str(num_T))
    f1.write('\n')
    f1.write("Total length of mystery virus genome = ")
    f1.write(str(len(assembled_ss)))


