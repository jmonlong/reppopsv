import glob
import math
import re
import os
import argparse
from pyfaidx import Fasta
import subprocess
import pysam
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def mergeRanges(ran):
    flag_rep = []
    for ii in xrange(len(ran)):
        if(ii in flag_rep):
            continue
        for jj in xrange(len(ran)):
            if(ran[ii][0] == ran[jj][0]
               and ran[ii][0] == ran[jj][0]):
                continue
            if(ran[ii][0] <= ran[jj][1]
               and ran[ii][1] >= ran[jj][0]):
                start = min(ran[ii][0], ran[jj][0])
                end = max(ran[ii][1], ran[jj][1])
                ran[ii] = [start, end]
                ran[jj] = [start, end]
                flag_rep.append(max(ii, jj))
    final_ran = []
    for ii in xrange(len(ran)):
        if(ii not in flag_rep):
            final_ran.append(ran[ii])
    return final_ran


def getReads(readids, reads_fn, reads_idx_fn):
    byte_os = []
    for idx_line in open(reads_idx_fn):
        idx_line = idx_line.rstrip('\n').split(':')
        if(idx_line[1].strip('>') in readids):
            byte_os.append(int(idx_line[0]))
    f = open(reads_fn)
    reads = {}
    for bos in byte_os:
        f.seek(bos)
        rid = f.next().rstrip('\n').strip('>')
        reads[rid] = f.next().rstrip('\n')
    f.close()
    return reads


def run_canu(input_file, genomeKb=30, errorRate=0.1, dir_name='test'):
    clean_cmd = ['rm', '-rf', dir_name]
    subprocess.check_output(clean_cmd)
    canu_cmd = ['canu', '-assemble', '-p', dir_name, '-d', dir_name,
                'genomeSize=' + str(genomeKb) + 'k',
                'errorRate=' + str(errorRate),
                'useGrid=false',
                '-pacbio-corrected', input_file]
    dump = open('/dev/null')
    subprocess.check_output(canu_cmd, stderr=dump)
    dump.close()
    # bubbles_fa = dir_name + '/' + dir_name + '.bubbles.fasta'
    contigs_fa = dir_name + '/' + dir_name + '.contigs.fasta'
    seq_desc = []
    for record in SeqIO.parse(contigs_fa, "fasta"):
        seq_desc.append([record.seq, record.description])
    # for record in SeqIO.parse(bubbles_fa, "fasta"):
    #     seq_desc.append([record.seq, record.description])
    subprocess.check_output(clean_cmd)
    return seq_desc


# Define arguments
parser = argparse.ArgumentParser(description='Validate CNV calls using PacBio reads.')
parser.add_argument('-reg', dest='region', help='the region of the CNV')
parser.add_argument('-ref', dest='ref', help='the reference genome fasta')
parser.add_argument('-reads', dest='reads', help='the PacBio reads (fasta)')
parser.add_argument('-readsidx', dest='readsidx',
                    help='the PacBio reads index (produced by `grep "^>" -b)` on the reads).')
parser.add_argument('-subreads', dest='subreads', help='the mapped PacBio subreads')
parser.add_argument('-bwaidx', dest='bwaidx', help='BWA index prefix')
parser.add_argument('-fl', dest='flank', type=int, default=1000, help='the flanks to add')
parser.add_argument('-outdir', dest='outdir', default='', help='the output directory')
args = parser.parse_args()

# Prefix for temp files
tempfn = 'tempPBV-' + str.replace(args.region, ":", "_")

ch_reg = str.split(args.region, ":")[0]
startend = str.split(args.region, ":")[1]
start_reg = int(str.split(startend, "-")[0])
end_reg = int(str.split(startend, "-")[1])

if(end_reg - start_reg > 30000):
    end_reg = start_reg + 30000
    args.region = ch_reg + ':' + str(start_reg) + '-' + str(end_reg)

start_reg += - args.flank
end_reg += args.flank

output_fa = args.outdir + '/pbv-MF-' + args.region.replace(":", "_") + '.fa'
if(os.path.isfile(output_fa)):
    print "Output file already there. Skipping..."
    exit(0)

print 'Region :' + args.region
now = time.time()

# Get subreads mapped to region
bamfile = pysam.AlignmentFile(args.subreads, "rb")
subreads = bamfile.fetch(reference=ch_reg, start=start_reg, end=end_reg)
pb_sel = {}
for read in subreads:
    pb_sel[read.query_name.split('--')[0]] = True

print 'Get ' + str(len(pb_sel)) + " subreads\t(" + str(int(time.time() - now)) + 's)'
now = time.time()

if(len(pb_sel) > 1500):
    print "Too many reads. Skipping..."
    exit(0)

# Get sequence of selected PacBio reads
pb_reads = getReads(pb_sel.keys(), args.reads, args.readsidx)

print 'Retrieve ' + str(len(pb_reads)) + ' PacBio reads\t' + str(int(time.time() - now)) + 's'
now = time.time()

# Write selected reads in Fasta
pbreads_file = tempfn + '-pbreads.fa'
f = open(pbreads_file, 'w')
for pbid in pb_reads.keys():
    f.write('>' + pbid + '\n')
    f.write(pb_reads[pbid])
    f.write("\n")
f.close()

# ----- Exonerate align and metrics

# Write reference in Fasta
reffa = Fasta(args.ref)
ex_ref_start = start_reg-args.flank-1000
ex_ref_end = end_reg+args.flank+1000
ref_region = reffa[str(ch_reg)][ex_ref_start:ex_ref_end]
ref_region = ref_region.seq
locref_file = tempfn + '-ref.fa'
f = open(locref_file, 'w')
f.write('>ref ' + ch_reg + ':' + str(ex_ref_start) + '-' + str(ex_ref_end) + '\n')
f.write(ref_region + '\n')
f.close()

# Prepare read length
read_len = {}
for pbr in pb_reads.keys():
    read_len[pbr] = len(pb_reads[pbr])

# Run Exonerate
exon_cmd = ['exonerate',
            '--model', 'affine:local',
            '--showalignment', 'FALSE', '-n', '1', '--gapextend', '0',
            '--gappedextension', 'FALSE',
            pbreads_file, locref_file]
ex_out = subprocess.check_output(exon_cmd)

print 'Prepare and run Exonerate\t' + str(int(time.time() - now)) + 's'
now = time.time()

# Parse output
read_blocks = {}
break_pos = [[], []]
for exl in ex_out.split('\n'):
    if(exl.find('vulgar') > -1):
        exls = exl.split(' ')
        bl = {'refS': int(exls[6]), 'refE': int(exls[7]), 'refStd': exls[8],
              'readS': int(exls[2]), 'readE': int(exls[3]), 'readStd': exls[4]}
        bl['size'] = max(abs(bl['refS'] - bl['refE']),
                         abs(bl['readS'] - bl['readE']))
        bl['readL'] = read_len[exls[1]]
        if(bl['size'] > 200):
            if(abs(bl['readS']-bl['readE']) < bl['readL'] * .8):
                break_pos[0].append(bl['refS'])
                break_pos[0].append(bl['refE'])
                break_pos[1].append(exls[1])
                break_pos[1].append(exls[1])
            try:
                read_blocks[exls[1]].append(bl)
            except KeyError:
                read_blocks[exls[1]] = [bl]


# Write reads with low-mapping quality in Fasta
pbreadslm_file = tempfn + '-pbreads-lowMap.fa'
reads_written = []
f = open(pbreadslm_file, 'w')
for pbid in break_pos[1]:
    if(pbid not in reads_written):
        f.write('>' + pbid + '\n')
        f.write(pb_reads[pbid])
        f.write("\n")
        reads_written.append(pbid)
f.close()

kb_reg = float(end_reg - start_reg) / 1000
canu_out = run_canu(pbreadslm_file, genomeKb=kb_reg,
                    errorRate=0.1, dir_name=tempfn + '_canu')

print 'Assembly\t' + str(int(time.time() - now)) + 's'
now = time.time()

# Find partially mapped reads breaking at the same location
break_nb = {}
for ii in xrange(len(break_pos[0])):
    try:
        break_nb[break_pos[0][ii]] += 1
    except KeyError:
        break_nb[break_pos[0][ii]] = 1
break_nb_ord = sorted(break_nb.keys(), key=lambda k: -break_nb[k])

# Build consensus for each cluster of reads
cons = {}
cl_ii = 0
while(break_nb[break_nb_ord[cl_ii]] > 2):
    rcl_file = tempfn + '-rcl.fa'
    f = open(rcl_file, 'w')
    nbr = 0
    for rii in xrange(len(break_pos[1])):
        if(break_pos[0][rii] == break_nb_ord[cl_ii]):
            nbr += 1
            pbr = break_pos[1][rii]
            pbseq = pb_reads[pbr]
            if(read_blocks[pbr][0]['refE'] == break_nb_ord[2]):
                readbp = read_blocks[pbr][0]['readE']
            else:
                readbp = read_blocks[pbr][0]['readS']
            pbseq = pbseq[max(0, readbp-400):min(len(pbseq), readbp+400)]
            if(read_blocks[pbr][0]['refStd'] == '-'):
                pbseq = MutableSeq(pbseq, generic_dna)
                pbseq.reverse_complement()
                pbseq = str(pbseq)
            f.write('>' + pbr + '\n')
            f.write(pbseq)
            f.write("\n")
    f.close()
    # Run Clustal
    clo_outfile = tempfn + '-clo-out.fa'
    clo_cmd = ['clustalo', '-i', rcl_file,
               '-o', clo_outfile, '--force']
    clo_out = subprocess.check_output(clo_cmd)
    # Get consensus
    msa_out = []
    for record in SeqIO.parse(clo_outfile, "fasta"):
        msa_out.append(str(record.seq))
    clo_cons = ""
    for cons_i in xrange(len(msa_out[0])):
        nbmsa = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '-': 0}
        for rr in xrange(len(msa_out)):
            nbmsa[msa_out[rr][cons_i]] += 1
            maxmsa = sorted(nbmsa.keys(), key=lambda k: -nbmsa[k])
            maxmsa = maxmsa[0]
        if(maxmsa != '-'):
            clo_cons += str(maxmsa)
    cons[cl_ii] = {'type': 'consensus', 'nbr': nbr,
                   'refPos': break_nb_ord[cl_ii], 'seq': clo_cons,
                   'size': len(clo_cons)}
    cl_ii += 1


print 'Build consensus\t' + str(int(time.time() - now)) + 's'
now = time.time()

# Add contigs to consensus
for seq_desc in canu_out:
    nbr = re.findall('reads=([0-9]+)', seq_desc[1])
    cons_ass = {'type': 'contig', 'seq': str(seq_desc[0]), 'refPos': 'NA',
                'size': len(seq_desc[0]), 'nbr': int(nbr[0])}
    cons[cl_ii] = cons_ass
    cl_ii += 1

if(len(cons.keys()) == 0):
    print 'No consensus or contig. Skipping...'
    # Remove temporary files
    filelist = glob.glob(tempfn + '*')
    for f in filelist:
        os.remove(f)
    exit(0)

# ----- Align assembly/consensus to local region
consensus_file = tempfn + '-cons.fa'
f = open(consensus_file, 'w')
for cc in cons:
    f.write('>' + str(cc) + '\n')
    f.write(cons[cc]['seq'])
    f.write("\n")
f.close()
exon_cmd = ['exonerate',
            '--model', 'affine:local',
            '--showalignment', 'FALSE', '-n', '1', '--gapextend', '0',
            '--gappedextension', 'FALSE',
            consensus_file, locref_file]
ex_out = subprocess.check_output(exon_cmd)

# Parse output
for exl in ex_out.split('\n'):
    if(exl.find('vulgar') > -1):
        exls = exl.split(' ')
        consS = int(exls[2])
        consE = int(exls[3])
        seqC = cons[int(exls[1])]['seq']
        sizeC = cons[int(exls[1])]['size']
        if(consS > sizeC - consE):
            bkPos = consS
        else:
            bkPos = consE
        piece_size = min(bkPos, sizeC-bkPos)
        cons[int(exls[1])]['bkPos'] = bkPos
        cons[int(exls[1])]['left'] = seqC[(bkPos-piece_size):bkPos]
        cons[int(exls[1])]['right'] = seqC[bkPos:(bkPos+piece_size)]
        cons[int(exls[1])]['mid'] = seqC[(bkPos-piece_size/2):(bkPos+piece_size/2)]
        cons[int(exls[1])]['left_match'] = 0
        cons[int(exls[1])]['right_match'] = 0
        cons[int(exls[1])]['left_chr'] = 'NA'
        cons[int(exls[1])]['right_chr'] = 'NA'
        cons[int(exls[1])]['left_pos'] = -1
        cons[int(exls[1])]['right_pos'] = -1
        cons[int(exls[1])]['mid_match'] = 0

# Re-map consensus pieces to genome
conspiece_file = tempfn + '-consPiece.fa'
f = open(conspiece_file, 'w')
for cc in cons:
    if('left' in cons[cc]):
        f.write('>' + str(cc) + '-left\n')
        f.write(cons[cc]['left'])
        f.write("\n")
        f.write('>' + str(cc) + '-right\n')
        f.write(cons[cc]['right'])
        f.write("\n")
        f.write('>' + str(cc) + '-mid\n')
        f.write(cons[cc]['mid'])
        f.write("\n")
f.close()

bwa_cmd = ['bwa', 'mem', '-v', '1', '-a',
           '-E', '0', '-L', '0', args.bwaidx, conspiece_file]
dump = open('/dev/null')
bwa_out = subprocess.check_output(bwa_cmd, stderr=dump)
dump.close()
bwa_out = bwa_out.split('\n')

chrs = [str(c) for c in range(22)] + ['X', 'Y']
for map in bwa_out:
    if(len(map) == 0 or map[0] == '@'):
        continue
    map = map.split('\t')
    readid = map[0]
    readid = readid.split("-")
    cid = int(readid[0])
    ch = map[2]
    pos = int(map[3])
    if(ch not in chrs):
        continue
    rc = int(map[1])
    if(readid[1] == 'left' and math.fmod(rc/16, 2) == 0):
        pos = pos + len(map[9])
    if(readid[1] == 'right' and math.fmod(rc/16, 2) == 1):
        pos = pos + len(map[9])
    if(readid[1] == 'mid'):
        pos = pos + len(map[9])/2
    al_match = re.findall('([0-9]+)M', map[5])
    al_match = [int(alm) for alm in al_match]
    prop_match = float(sum(al_match)) / len(cons[cid]['right'])
    if(cons[cid][readid[1] + '_match'] < prop_match):
        cons[cid][readid[1] + '_match'] = prop_match
        cons[cid][readid[1] + '_chr'] = ch
        cons[cid][readid[1] + '_pos'] = pos


print 'Remap consensus\t' + str(int(time.time() - now)) + 's'
now = time.time()


recs = []
for cid in cons:
    if('left' in cons[cid]):
        desc = 'reads=' + str(cons[cid]['nbr'])
        desc += ' length=' + str(len(cons[cid]['seq']))
        desc += ' class=' + cons[cid]['type']
        desc += ' region=' + args.region
        desc += ' midM=' + str(cons[cid]['mid_match'])
        desc += ' leftM=' + str(cons[cid]['left_match'])
        desc += ' rightM=' + str(cons[cid]['right_match'])
        desc += ' leftChr=' + str(cons[cid]['left_chr'])
        desc += ' rightChr=' + str(cons[cid]['right_chr'])
        desc += ' leftPos=' + str(cons[cid]['left_pos'])
        desc += ' rightPos=' + str(cons[cid]['right_pos'])
        desc += ' refPos=' + str(cons[cid]['refPos'])
        desc += ' bkPos=' + str(cons[cid]['bkPos'])
        recs.append(SeqRecord(MutableSeq(cons[cid]['seq'], generic_dna),
                              id='cons' + str(cid),
                              description=desc))

SeqIO.write(recs, output_fa, "fasta")

# Remove temporary files
filelist = glob.glob(tempfn + '*')
for f in filelist:
    os.remove(f)

print 'Write.output files\t' + str(int(time.time() - now)) + 's'
