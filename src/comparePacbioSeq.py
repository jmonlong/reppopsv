import argparse
from pyfaidx import Fasta
from Bio import SeqIO
import subprocess
import os
import re


def plotMummer(seq0, seq1, output_prefix, flank):
    ref_file = 'tempRef.fasta'
    contig_file = 'tempContig.fasta'
    ff = open(ref_file, 'w')
    ff.write('>ref\n' + str(seq0) + '\n')
    ff.close()
    ff = open(contig_file, 'w')
    ff.write('>assembly\n' + str(seq1) + '\n')
    ff.close()
    mum_cmd = ['mummer', '-mum', '-b', '-c', ref_file, contig_file]
    dump = open('/dev/null')
    mum_out = subprocess.check_output(mum_cmd, stderr=dump)
    ff = open('mummer.mums', 'w')
    ff.write(mum_out)
    ff.close()
    mum_cmd = ['mummerplot', '-postscript', '-p', output_prefix, 'mummer.mums']
    mum_out = subprocess.check_output(mum_cmd, stderr=dump)
    os.remove('mummer.mums')
    os.remove(ref_file)
    os.remove(contig_file)
    # Add vertical lines
    op = output_prefix.split('_')
    v2 = int(op[2]) - int(op[1]) + flank
    gpfile2 = output_prefix + '_2.gp'
    gpcon2 = open(gpfile2, 'w')
    finish = False
    for line in open(output_prefix + '.gp', 'r'):
        if('set' in line):
            gpcon2.write(line)
        elif (not finish):
            gpcon2.write('set arrow from ', str(flank) + ', graph 0 to ' + flank + ', graph 1 nohead\n')
            gpcon2.write('set arrow from ' + str(v2) + ', graph 0 to ' + str(v2) + ', graph 1 nohead\n')
            finish = True
        if(finish):
            gpcon2.write(line)
    gpcon2.close()
    gnuplot_cmd = ['gnuplot', gpfile2]
    gnuplot_out = subprocess.check_output(gnuplot_cmd, stderr=dump)
    dump.close()


# Define arguments
parser = argparse.ArgumentParser(description='Retrieve sequence of blast hits and compare them')
parser.add_argument('-ass', dest='assfile', help='the pacbio assembled sequences')
parser.add_argument('-ref', dest='reffile', help='the reference genome file')
parser.add_argument('-fl', dest='flank', type=int, help='the flanks size')
parser.add_argument('-outdir', dest='outdir', help='where to put the dotplots')
args = parser.parse_args()

reffa = Fasta(args.reffile)

# Get the sequence of the assembled regions
assseq = []
for record in SeqIO.parse(args.assfile, "fasta"):
    desc = record.description
    midM = float(re.search('midM=([^ ]*)', desc).group(1))
    leftM = float(re.search('leftM=([^ ]*)', desc).group(1))
    rightM = float(re.search('rightM=([^ ]*)', desc).group(1))
    if(midM < .8 and (midM < leftM or midM < rightM)):
        assseq.append(record)

refseq = ''
coord = re.search('pbv-MF-(.*).fa', args.assfile).group(1)
coordS = coord.split('_')
start = int(coordS[1]) - args.flank
end = int(coordS[2]) + args.flank
if(len(assseq) > 0):
    refseq = reffa[coordS[0]][start:end]
else:
    exit()

for record in assseq:
    plotMummer(refseq, str(record.seq), args.outdir + '/' + coord + '-' + record.id, args.flank)
