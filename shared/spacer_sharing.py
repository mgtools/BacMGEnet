#!/usr/bin/env python3

import sys
import re

if len(sys.argv) < 5:
    sys.exit(sys.argv[0] + " spacer-clr-file crt-file outfile1 outfile2")

inf = open(sys.argv[1], "r")
cidx = 0
spacerlist = []
seqlen = {}
for aline in inf:
    if aline[0] == '>':
        cidx += 1
    else:
        subs = aline.strip().split()
        spacerid = subs[2][1:-3]
        subs2 = spacerid.split(":")
        arrayid = ":".join(subs2[:-1])
        spacerpos = int(subs2[-1][1:])
        spacerlist.append([cidx, arrayid, spacerpos, spacerid])
        seqid = subs2[-3]
        if seqid not in seqlen:
            seqlen[seqid] = 0
inf.close()

inf = open(sys.argv[2], "r")
for aline in inf:
    m = re.match("SEQ: (\S+) length (\S+)", aline)
    if m:
        seqlen[m.groups()[0]] = int(m.groups()[1])
inf.close()

#sort spacer according to their sequence id and then position
spacer_sorted = sorted(spacerlist, key = lambda x: (x[1], x[2]))

beg = 0
out1 = open(sys.argv[3], "w")
out2 = open(sys.argv[4], "w")
while beg < len(spacer_sorted):
    for end in range(beg, len(spacer_sorted)):
        if (end == len(spacer_sorted) - 1) or (spacer_sorted[end + 1][1] != spacer_sorted[end][1]):
            break
    if end == beg:
        avelen = 70
    else:
        #number of spacers - 1
        avelen = (spacer_sorted[end][2] - spacer_sorted[beg][2]) / (end - beg)
    if spacer_sorted[beg][2] >= avelen * 1.5:
        tmp = f"{spacer_sorted[beg][1]} end"
    else:
        tmp = f"{spacer_sorted[beg][1]} unk"
    tmp2 = tmp
    for idx in range(beg, end + 1):
        tmp += f" {spacer_sorted[idx][0]}"
        tmp2 += f" {spacer_sorted[idx][0]}[{spacer_sorted[idx][3]}]"
    subs = spacer_sorted[beg][1].split(":")
    seqid = subs[-2] 
    if (seqlen[seqid] - spacer_sorted[end][2]) >= avelen * 2.0:
        tmp += f" end"
        tmp2 += f" end"
    else:
        tmp += f" unk"
        tmp2 += f" unk"
    out1.write(tmp + "\n")
    out2.write(tmp2 + "\n")

    beg = end + 1
out1.close()
out2.close()
