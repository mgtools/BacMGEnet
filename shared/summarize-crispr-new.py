from __future__ import print_function
#!/usr/bin/env python3

#Yuzhen Ye (latest update, July 2021)
#run CRISPR identification tool for all the samples
#Oct 11, 2012 add description of the sequences to crispr ann files

import os
import urllib
import sys
import re
import glob
import operator

def getRepeatSpacer(crtfile, maxmiss, mincopy, minlen, maxlen, seqidprefix=""):
        repeatID, repeatSeq, spacerID, spacerSeq = [], [], [], []
        crisprID, crisprDes, crisprSeq, crisprSeqAnn, crisprCon = [], [], [], [], []
        contig, contigbase, crispr = 0, 0, 0

        infile = open(crtfile, "r")
        crisprnum = 0
        for aline in infile:
                subs = aline.split()
                if aline[:3] == 'SEQ':
                        #seqid = aline[5:-1]
                        subs2 = aline[5:-1].split()
                        seqid = subs2[0]
                        seqdes = " ".join(subs2[1:])
                        crisprnum = 0
                elif aline[:6] == 'CRISPR':
                        repeatlist, repeatpos, repeatcopy, spacerlist, spacerpos, repeatspacer, repeatspacerann = [], [], [], [], [], "", ""
                        add = 0
                        for aline2 in infile:
                                add += 1
                                if add < 3:
                                        continue
                                if aline2[0:3] == '---':
                                        #check if the crispr matches the contraints
                                        break
                                subs2 = aline2.split()
                                repeattmp = subs2[1].replace("-", "")
                                if len(subs2) > 2:
                                        pos, repeat, spacer = int(subs2[0]), subs2[1], subs2[2]
                                        if (repeat[0] != '-') and (repeat[-1] != '-'):
                                                if repeat not in repeatlist:
                                                        repeatlist.append(repeat)
                                                        repeatpos.append(str(pos))
                                                        repeatcopy.append(1)
                                                else:
                                                        whichrepeat = repeatlist.index(repeat)
                                                        repeatpos[whichrepeat] += "+" + str(pos)
                                                        repeatcopy[whichrepeat] += 1
                                        if (spacer not in spacerlist):
                                                spacerlist.append(spacer)
                                                spacerpos.append(pos + len(repeat))
                                        repeatspacer += repeattmp + subs2[2]
                                        #repeatspacerann += "r" + str(len(repeattmp)) + "s" + str(len(subs2[2]))
                                        repeatspacerann += repeattmp.upper() + subs2[2].lower()
                                else:
                                        #July 2012
                                        pos, repeat = int(subs2[0]), subs2[1]
                                        if (repeat[0] != '-') and (repeat[-1] != '-'):
                                                if repeat not in repeatlist:
                                                        repeatlist.append(repeat)
                                                        repeatpos.append(str(pos))
                                                        repeatcopy.append(1)
                                                else:
                                                        whichrepeat = repeatlist.index(repeat)
                                                        repeatpos[whichrepeat] += "+" + str(pos)
                                                        repeatcopy[whichrepeat] += 1
                                        #end of July 2021
                                        repeatspacer += repeattmp
                                        #repeatspacerann += "r" + str(len(repeattmp))
                                        repeatspacerann += repeattmp.upper()
                elif "Repeat-con" in aline:
                        subs = aline.strip().split()    
                        currcons = subs[1]

                        fullrepeat = 0
                        diff = 0
                        valid = 1
                        print("currcons", currcons, "len", len(currcons), "minlen", minlen, "maxlen", maxlen)
                        if (len(currcons) < minlen) or (len(currcons) > maxlen):
                                valid = 0
                        else:
                                for idx in range(len(repeatlist)):
                                        if len(repeatlist[idx]) == len(currcons):
                                                fullrepeat += repeatcopy[idx] 
                                                for pos in range(len(currcons)):
                                                        if repeatlist[idx][pos] != currcons[pos]:
                                                                diff += 1
                                if fullrepeat == 0:
                                    sys.exit(f"Error with {crtfile} file -- zero fullrepeat")
                                print("fullrepeat", fullrepeat, "mincopy", mincopy, "diff", diff, "diff-per-repeat", diff * 1.0 / fullrepeat)
                                if (fullrepeat < mincopy) or (diff > maxmiss * fullrepeat): 
                                        valid = 0
                        print("seqid", seqid, "valid", valid)
                        if valid == 1:
                                crisprnum += 1          
                                allrepeatpos = []
                                for idx in range(len(repeatlist)):
                                        repeatID.append(seqidprefix + seqid + ":c" + str(crisprnum) + ":p" + str(repeatpos[idx])) 
                                        repeatSeq.append(repeatlist[idx])
                                        tmp = repeatpos[idx].split("+")
                                        for one in tmp:
                                            allrepeatpos.append(int(one))
                                for idx in range(len(spacerlist)):
                                        spacerID.append(seqidprefix + seqid + ":c" + str(crisprnum) + ":p" + str(spacerpos[idx])) 
                                        spacerSeq.append(spacerlist[idx])
                                #YY, July 2021 add the start position of the array so it can be matched to the crispr arrays in -sm.gff file
                                #crisprID.append(seqidprefix + seqid + ":c" + str(crisprnum))
                                crisprID.append(f"{seqidprefix}{seqid}:c{crisprnum}:p{min(allrepeatpos)}") ##append the range of the repeats
                                crisprDes.append(seqdes)
                                #repeatspacer = repeatspacer.replace("-", "") 
                                crisprSeq.append(repeatspacer)
                                crisprSeqAnn.append(repeatspacerann)    
                
                                crisprCon.append(currcons)

                elif "Total sequences checked" in aline:
                        contig = subs[-1] 
                elif "Total bases checked" in aline:
                        contigbase = subs[-1] 
                elif "Total sequences with predicted" in aline:
                        crispr = subs[-1] 

        infile.close()

        return (contig, contigbase, crispr, repeatID, repeatSeq, spacerID, spacerSeq, crisprID, crisprDes, crisprSeq, crisprSeqAnn, crisprCon)

def saveSeq(id, seq, linelen, out):
        for idx in range(len(id)):
                out.write(">" + id[idx] + "\n")
                if linelen == 0:
                        out.write(seq[idx] + "\n")
                else:
                        for pos in range(0, len(seq[idx]), linelen):
                                out.write(seq[idx][pos:pos+linelen] + "\n")

def getConsensNR(repeatcon_all_seq, repeatcon_all_id):
        tosort = []
        tot = len(repeatcon_all_seq)
        for idx in range(tot):
                tosort.append((len(repeatcon_all_seq[idx]), repeatcon_all_seq[idx], repeatcon_all_id[idx]))             
        consorted = sorted(tosort, key=operator.itemgetter(0))
        repeatcon_all_seq_sorted = map(operator.itemgetter(1), consorted)
        repeatcon_all_seq_sorted = repeatcon_all_seq_sorted[::-1]
        repeatcon_all_id_sorted = map(operator.itemgetter(2), consorted)
        repeatcon_all_id_sorted = repeatcon_all_id_sorted[::-1]
        repeatcon_valid = [1] * tot

        repeatcon_nr_seq = []
        repeatcon_nr_id = []
        for idx1 in range(tot):
                if repeatcon_valid[idx1] == 0:
                        continue
                repeatcon_nr_seq.append(repeatcon_all_seq_sorted[idx1])
                repeatcon_nr_id.append(repeatcon_all_id_sorted[idx1])
                for idx2 in range(idx1 + 1, tot):
                        if repeatcon_valid[idx2] == 0:
                                continue
                        if repeatcon_all_seq_sorted[idx2] in repeatcon_all_seq_sorted[idx1]: 
                                repeatcon_nr_id[-1] += ":" + repeatcon_all_id_sorted[idx2]
                                repeatcon_valid[idx2] = 0

        return (repeatcon_nr_id, repeatcon_nr_seq)

def main():

        listfile, filename = '', ''
        suffix = "*.crt"
        dir = "./"
        repeatfile, spacerfile, crisprfile, crisprannfile, repeatconfile, seqidprefix = '', '', '', '', '', ''
        repeatcon_nr_file = ''
        maxmiss, mincopy, minlen, maxlen = 10000, 2, 0, 1000
        #default: NO contraints on miss-match & copy number & minimum length of the repeat & maxlen of the repeat

        for idx in range(len(sys.argv)):
                if sys.argv[idx] == "-i" and len(sys.argv) > idx + 1:
                        listfile = sys.argv[idx + 1]    
                elif sys.argv[idx] == "-f" and len(sys.argv) > idx + 1:
                        filename = sys.argv[idx + 1]    
                elif sys.argv[idx] == "-suffix" and len(sys.argv) > idx + 1:
                        suffix = sys.argv[idx + 1]      
                elif sys.argv[idx] == "-dir" and len(sys.argv) > idx + 1:
                        dir = sys.argv[idx + 1] 
                elif ("-maxmis" in sys.argv[idx]) and len(sys.argv) > idx + 1:
                        maxmiss = int(sys.argv[idx + 1])
                elif sys.argv[idx] == "-mincopy" and len(sys.argv) > idx + 1:
                        mincopy = int(sys.argv[idx + 1])
                elif sys.argv[idx] == "-minlen" and len(sys.argv) > idx + 1:
                        minlen = int(sys.argv[idx + 1])
                elif sys.argv[idx] == "-maxlen" and len(sys.argv) > idx + 1:
                        maxlen = int(sys.argv[idx + 1])
                elif sys.argv[idx] == "-repeat" and len(sys.argv) > idx + 1:
                        repeatfile = sys.argv[idx + 1]  
                elif sys.argv[idx] == "-consensus" and len(sys.argv) > idx + 1:
                        repeatconfile = sys.argv[idx + 1]
                elif sys.argv[idx] == "-consensus_nr" and len(sys.argv) > idx + 1:
                        repeatcon_nr_file = sys.argv[idx + 1]
                elif sys.argv[idx] == "-spacer" and len(sys.argv) > idx + 1:
                        spacerfile = sys.argv[idx + 1]  
                elif sys.argv[idx] == "-crispr" and len(sys.argv) > idx + 1:
                        crisprfile = sys.argv[idx + 1]  
                elif sys.argv[idx] == "-crisprann" and len(sys.argv) > idx + 1:
                        crisprannfile = sys.argv[idx + 1]
                elif sys.argv[idx] == "-seqid_prefix" and len(sys.argv) > idx + 1:
                        seqidprefix = sys.argv[idx + 1]   
                idx += 1

        if not (listfile or filename):
                print("summarize-crispr-new.py <-f file> or <-i list-file> <-suffix suffix(default *.crt)> <-dir directory-of-crispr-predictions(default ./) <-repeat repeat-file> <-spacer spacer-file>")
                print("Optional parameters (to include crisprs that meet the following constraints): ")
                print("-maxmis max-mismatches-allowed-per-repeat; default off")
                print("-mincopy min-copy-of-the-intact-repeats; default off")
                print("-minlen min-len-of-the-consensus-of-the-repeats; default off")
                print("-maxlen min-len-of-the-consensus-of-the-repeats; default off")
                print(" example 1: summarize-crispr.py -i /data/hmp/WholeMetagenomic/SRS-info/srs_info_sortbysubject.txt -suffix *.contigs.fa.crt -dir /data/hmp/WholeMetagenomic/CRISPR/")
                print(" example 2: summarize-crispr.py -f test.crt -repeat test-repeat.seq -consensus test-repeatcon.seq -spacer test-spacer.seq -crispr test-crispr.seq -crisprann test-crispr.ann")
                sys.exit()

        filelist, deslist = [], []
        if filename:
                filelist.append(filename)
                deslist.append("")

        if listfile:
                infile = open(listfile, "r")
                for aline in infile:
                        aline = aline.strip()
                        subs = aline.split("\t")
                        srs = subs[0]
                        dtmp = dir + "/" + srs + suffix 
                        print("srs", srs, "suffix", suffix, "dtmp", dtmp)
                        tmps = glob.glob(dtmp)

                        for atmp in tmps:
                                filelist.append(atmp)
                                deslist.append(aline)
                                #only one assembly for each id 
                                break
                infile.close()
        
        if repeatfile:
                repeatout = open(repeatfile, "w")
        if repeatconfile:
                repeatconout = open(repeatconfile, "w")
        if spacerfile:
                spacerout = open(spacerfile, "w")
        if crisprfile:
                crisprout = open(crisprfile, "w")
        if crisprannfile:
                crisprannout = open(crisprannfile, "w")

        repeatcon_all_seq, repeatcon_all_id = [], []
        for idx in range(len(filelist)):
                file, des = filelist[idx], deslist[idx]
                (contig, contigbase, crispr, repeatID, repeatSeq, spacerID, spacerSeq, crisprID, crisprDes, crisprSeq, crisprSeqAnn, crisprCon) = getRepeatSpacer(file, maxmiss, mincopy, minlen, maxlen, seqidprefix)
                print("%s\t%s\t%s\t%s" % (crispr, contig, contigbase, des))
                print("repeat number: ", len(repeatID))
                print("spacer number: ", len(spacerID))
                repeatcon_all_seq.extend(crisprCon)
                repeatcon_all_id.extend(crisprID)
                if repeatfile:
                        saveSeq(repeatID, repeatSeq, 0, repeatout)
                if repeatconfile:
                        saveSeq(crisprID, crisprCon, 0, repeatconout)
                if spacerfile:
                        saveSeq(spacerID, spacerSeq, 0, spacerout)
                if crisprfile:
                        saveSeq(crisprID, crisprSeq, 70, crisprout)
                if crisprannfile:
                        idnew = []
                        for idx in range(len(crisprID)):
                                idnew.append(crisprID[idx] + " " + crisprCon[idx] + " " + crisprDes[idx])
                        saveSeq(idnew, crisprSeqAnn, 70, crisprannout)
        if repeatfile:
                repeatout.close()
        if spacerfile:
                spacerout.close()
        if crisprfile:
                crisprout.close()
        if crisprannfile:
                crisprannout.close()
        if repeatcon_nr_file:
                (repeatcon_nr_id, repeatcon_nr_seq) = getConsensNR(repeatcon_all_seq, repeatcon_all_id)
                file = open(repeatcon_nr_file, "w")
                saveSeq(repeatcon_nr_id, repeatcon_nr_seq, 0, file)
                file.close()
main()
