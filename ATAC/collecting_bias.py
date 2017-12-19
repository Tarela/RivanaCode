import os,sys
def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
  #      print seq
        nmer_seq[seq] = []
    del allseq
    return nmer_seq

def addbias(rawdict,encodingdict, infile):
    inf = open(infile)
    tmpdict = {}
    for line in inf:
        ll = line.split()
        seq = ll[0]
        raw = ll[1]
        encoding = ll[2]
        tmpdict[seq] = ll[1:]
    inf.close()
    for seqtype in rawdict.keys():
        if tmpdict.has_key(seqtype):
            rawdict[seqtype].append( str(round(float(tmpdict[seqtype][0]),4)) )
            encodingdict[seqtype].append(str(round(float(tmpdict[seqtype][1]),4)) )
        else:
            rawdict[seqtype].append(str(-8))
            encodingdict[seqtype].append(str(-8))
    return
    
            

### 8mer bias
Encodingbias = make_nmer_dict(8)
rawbias = make_nmer_dict(8)

total_sample_list = []

# mouse tissue DNase
folder = '/nv/vol190/zanglab/sh8tv/Project/scATAC/Data/ENCODE_DNase/mouse_tissue/'
suffix = '_8mer_raw_Encoding_bias.txt'
samplelist = ["brainE145","brainE185","cerebellum8w","embryoE115","fatPad8w","gonadalFatPad8w","heart8w","kidney8w","largeIntestine8w","liver8w","lung8w","skeletalMuscle8w","spleen8w","telencephalon8w","thymus8w"]
for sample in samplelist:
    total_sample_list.append("mmDNase_"+sample)
    this_file = folder + sample + suffix
    addbias(rawbias, Encodingbias, this_file)

# human cell line DNase
folder = '/nv/vol190/zanglab/sh8tv/Project/scATAC/Data/ENCODE_DNase/human_cellline/'
suffix = '_DNase_8mer_raw_Encoding_bias.txt'
samplelist = ["H1","K562","GM12878"]
for sample in samplelist:
    total_sample_list.append("hsDNase_"+sample)
    this_file = folder + sample + suffix
    addbias(rawbias, Encodingbias, this_file)



### output
outf1 = open('summary_rawbias_8mer.txt','w')
outf2 = open('summary_encodingbias_8mer.txt','w')

newll = ['seqtype'] + total_sample_list
outf1.write("\t".join(newll)+"\n")
outf2.write("\t".join(newll)+"\n")

for SeqT in sorted(rawbias.keys()):
    newll1 = [SeqT] + rawbias[SeqT]
    newll2 = [SeqT] + Encodingbias[SeqT]
    outf1.write("\t".join(newll1)+"\n")
    outf2.write("\t".join(newll2)+"\n")

outf1.close()
outf2.close()






















