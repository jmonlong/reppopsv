import argparse
from pyfaidx import Fasta

# Define arguments
parser = argparse.ArgumentParser(description='Retrieve the reference sequence flanking regions')
parser.add_argument('-bed', dest='bed', help='the file with the regions coordinates')
parser.add_argument('-ref', dest='ref', help='the reference genome fasta')
parser.add_argument('-o', dest='outfile', help='the output fasta file')
parser.add_argument('-addchr', dest='addchr', action='store_true',
                    help='should chr be added to the chr names')
parser.add_argument('-fl', dest='flank', type=int, default=10000,
                    help='the flanks to add')
args = parser.parse_args()

reffa = Fasta(args.ref)
inf = open(args.bed, 'r')
headers = inf.next()
outf = open(args.outfile, 'w')

for region in inf:
    region = region.rstrip('\n').split('\t')
    if args.addchr:
        region[0] = 'chr' + region[0]
    reg_label = region[0] + '_' + region[1] + '_' + region[2]
    reg_start = int(region[1]) - args.flank
    reg_end = int(region[1])
    reg_seq = reffa[region[0]][reg_start:reg_end]
    reg_seq = reg_seq.seq
    outf.write('>' + reg_label + '_flankU\n')
    outf.write(reg_seq + '\n')
    reg_start = int(region[2])
    reg_end = int(region[2]) + args.flank
    reg_seq = reffa[region[0]][reg_start:reg_end]
    reg_seq = reg_seq.seq
    outf.write('>' + reg_label + '_flankD\n')
    outf.write(reg_seq + '\n')
    reg_start = int(region[1]) - args.flank
    reg_end = int(region[2]) + args.flank
    reg_seq = reffa[region[0]][reg_start:reg_end]
    reg_seq = reg_seq.seq
    outf.write('>' + reg_label + '_total\n')
    outf.write(reg_seq + '\n')

# Close connections
inf.close()
outf.close()
