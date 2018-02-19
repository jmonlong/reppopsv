import argparse
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import subprocess
import os


def plotMummer(seq0, seq1, output_prefix):
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
    v2 = int(op[2]) - int(op[1]) + 50000
    gpfile2 = output_prefix + '_2.gp'
    gpcon2 = open(gpfile2, 'w')
    finish = False
    for line in open(output_prefix + '.gp', 'r'):
        if('set' in line):
            gpcon2.write(line)
        elif (not finish):
            gpcon2.write('set arrow from 50000, graph 0 to 50000, graph 1 nohead\n')
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
parser.add_argument('-ass', dest='assfile', help='the assembly fasta (indexed)')
parser.add_argument('-reg', dest='regfile', help='the sequence of the regions')
parser.add_argument('-blast', dest='blastfile', help='the output from Blast')
parser.add_argument('-blastpairs', dest='blastpairs', default='', help='the output from Blast')
parser.add_argument('-outdir', dest='outdir', help='where to put the dotplots')
args = parser.parse_args()

assfa = Fasta(args.assfile)
blast_con = open(args.blastfile, 'r')
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

# Get the sequence of the regions of interest
regions = {}
for record in SeqIO.parse(args.regfile, "fasta"):
    regions[record.id] = record.seq

# Get the pairs to plot
blastpairs = []
if(args.blastpairs != ''):
    for line in open(args.blastpairs, 'r'):
        line = line.rstrip('\n').split('\t')
        blastpairs.append(line[0] + ' ' + line[1])

for bline in blast_con:
    bline = bline.rstrip('\n').split('\t')
    if('_total' not in bline[0]):
        continue
    if(len(blastpairs) > 0):
        if(bline[0] + ' ' + bline[1] not in blastpairs):
            continue
    regseq = regions[bline[0]]
    reglen = len(regseq)
    qstart = int(bline[6])
    qend = int(bline[7])
    sstart = int(bline[8])
    send = int(bline[9])
    revcomp = False
    if(sstart < send):
        sstart = sstart - qstart
        send = send + reglen - qend
    else:
        sswitch = sstart
        sstart = send - reglen + qend
        send = sswitch + qstart
        revcomp = True
    if(sstart < 0):
        sstart = 0
    graph_title = bline[0] + ' ' + bline[1] + '_' + str(sstart) + '_'
    if(send > len(assfa[bline[1]])):
        send = len(assfa[bline[1]]) - 1
        graph_title = graph_title + str(send) + ' (end)'
    else:
        graph_title = graph_title + str(send)
    print str(sstart) + '-' + str(send)
    assseq = assfa[bline[1]][sstart:send]
    if(revcomp):
        assseq = MutableSeq(str(assseq), generic_dna)
        assseq.reverse_complement()
        assseq = str(assseq)
    plotMummer(regseq, assseq, args.outdir + '/' + bline[0] + '-' + bline[1])

blast_con.close()
