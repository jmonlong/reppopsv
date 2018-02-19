import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Split PacBio reads into smaller subreads.')
parser.add_argument('-pbfa', dest='pbfa', help='a faste file with PacBio reads')
parser.add_argument('-k', dest='k', default=200, type=int, help='the size of the subreads')
args = parser.parse_args()

for pbr in SeqIO.parse(args.pbfa, "fasta"):
    for rs in xrange(len(pbr)/args.k):
        seq = pbr.seq[rs*args.k:(rs+1)*args.k-1]
        print '@' + pbr.id + '--' + str(rs)
        print str(seq)
        print '+'
        print '~' * len(seq)
