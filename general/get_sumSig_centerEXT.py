#!/usr/bin/env python
'''
Created on XXXX-XX-XX

@author: Tarela
'''
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging,time
import string,numpy
from bx.bbi.bigwig_file import BigWigFile

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

#import subprocess
#def sp(cmd):
#    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
#    ac = a.communicate()
#    return ac
#def bwsig_pattern(bwfile,chrm,start,end,points):
#    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
#    sigPat = sp(cmd)[0].strip().split("\t")
#    return sigPat

def get_signal(inputfile,output,bwfiles,bwfolder,extend):
    signalbw = bwfiles.strip().strip(',').split(',')

    if not bwfolder:
        bwfolder = "./"
    if not bwfolder.endswith('/'):
        bwfolder += '/'

    bwHs = []
    for sb in signalbw:
        if sb.startswith('/'):
            bwHs.append(BigWigFile(open(sb, 'rb')))
        else:
            bwHs.append(BigWigFile(open(bwfolder + sb, 'rb')))

    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        if "_" in ll[0]:
            continue
        #center = (int(ll[1]) + int(ll[2]))/2
        #S = max(0,center - extend)
        #E = center + extend
        C = (int(ll[1]) + int(ll[2]) ) /2
        S = C - extend
        E = C + extend
        for bwHandle in bwHs:
            try:
                signal=(bwHandle.summarize(ll[0],S,E,E-S))
                if type(signal.sum_data) == None:
                    #print 'c1',line
                    addsig = 0
                else:
                    addsig = sum(signal.sum_data)/(E-S)#float(signal.sum_data/signal.valid_count)
            except:
                #print 'c2',line
                addsig = 0#'nan'
               # ll.extend(list(signal.sum_data/signal.valid_count))
            ll.append(addsig)
        outf.write("\t".join(map(str,ll))+"\n")
    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """calculate the average signal within full len of input peak, multiple bw files separated by comma"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="input bed file")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="output bedEX file")
    optparser.add_option("-w","--bws",dest="bws",type="str",default = "",
                         help="bigwig files, comma separate")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",default = "",
                         help="folder of bigwig files")
    optparser.add_option("--ext",dest="ext",type="int",default=1000,
                         help="ext size from center,default=1000")

#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    bws = options.bws
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,bws,options.bwfolder,options.ext)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


