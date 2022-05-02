from __future__ import print_function
#!/usr/bin/env python3
import sys
import os
import subprocess
import argparse
import string

python = "python3"
cmd=sys.argv[0]
ppath=os.path.dirname(sys.argv[0])
if not ppath:
    ppath = "./"
crisprone = ppath + "/../CRISPRone/"
print(f"path {ppath} crisprone {crisprone}")
if os.path.exists(crisprone + "/local") == False:
    sys.exit("CRISPRone directory not set correctly")
cdhit = f"{crisprone}/bin/cd-hit-v4.6.1-2012-08-27/cd-hit-est"

parser = argparse.ArgumentParser()
parser.add_argument("--input", metavar="file", help="a file with the list of genomes", required=True)
parser.add_argument("--prog", metavar="CRISPRone/metaCRT/CRISPRAlign", help="what program to use to predict CRISPR (and cas genes)", default="CRISPRone", required=False)
parser.add_argument("--out", metavar="path", help="specify the path where the outputs will go; use ./ for current folder", required=True)
parser.add_argument("--repeat", metavar="repeat-file", help="if provided, sequences included in the file are to used as the guide repeats; default: off", default="", required=False)
parser.add_argument("--maxmis", metavar="max-mismatch-allowed", help="used when --repeat is given", default=-1, required=False, type=int)
parser.add_argument("--nr", metavar="non-redundant cutoff", help="default 0.85", default=0.85, type=float)
parser.add_argument("--color", metavar="a file specifying colors of nodes in the spacer-graph", required=False)
parser.add_argument("--cleanup", metavar="True/False", help="for debugging only; turn off to keep immediate files", default=True, required=False)
parser.add_argument("--checkmock", metavar="True/False", help="for debugging only: check mock or not, default True", default=True, required=False)
args = parser.parse_args()

if args.prog.lower() not in ["crisprone", "metacrt", "crispralign"]:
    sys.exit("Invalid prog")
if (args.prog.lower() == "crispralign") and (not args.repeat):
    sys.exit("--repeat required for --prog CRISPRAlign")

if not os.path.exists(args.out):
    os.mkdir(args.out)
if not os.path.exists(args.out + "/crispr"):
    os.mkdir(args.out + "/crispr")
if not os.path.exists(args.out + "/spacer_graph"):
    os.mkdir(args.out + "/spacer_graph")

#get the list of sequences (and metadata)
genomelist = []
infile = open(args.input, "r")
#one line per genome
#format: ID,tag,fna,gff
for aline in infile:
    subs = aline.strip().split()
    if len(subs) < 3:
        continue
    genomelist.append(subs)
infile.close()

if len(genomelist) < 1:
    print("no genomes were found; check your input and see if it is the correct format")
    sys.exit(0)

#files are cleaned up
all_crt, all_spacer, tmp_spacer = f"{args.out}/spacer_graph/all.crt", f"{args.out}/spacer_graph/spacer_all.fa", f"{args.out}/spacer_graph/tmp_spacer.fa"
all_repeat, tmp_repeat = f"{args.out}/spacer_graph/repeat_all.fa", f"{args.out}/spacer_graph/tmp_repeat.fa" 
subprocess.check_output("rm -f " + all_crt + " " + all_spacer, shell=True)

#predict CRISPR and collect spacers
spath=f"{ppath}/../shared"
for one in genomelist:
    sid, tag, fna = one[:3]
    print(f"working in progress: {sid} tag {tag} ...")
    if len(one) >= 4:
        gff = one[3]
    if args.prog.lower() == 'crisprone':
        crt_output = f"{args.out}/crispr/{sid}/{sid}.crt"
        cmd = f"{python} {crisprone}/crisprone-local.py --fna {fna} --outbase {args.out}/crispr --outprefix {sid} 1>out.log 2>error.log"
        check_result=subprocess.check_output(cmd, shell=True)
    elif args.prog.lower() == 'metacrt':
        crt_output = f"{args.out}/crispr/{sid}/{sid}.crt"
        cmd = f"java -cp {crisprone}/bin/metaCRT.jar crt {fna} + {crt_output}"
        check_result=subprocess.check_output(cmd, shell=True)
    else:
        print("run CRISPRAlign...")
        crt_output = f"{args.out}/crispr/{sid}.crt"
        if args.maxmis != -1:
            cmd = f"{crisprone}/bin/CrisprAlign -r {args.repeat} -s {fna} -d {args.maxmis} -partial -crt {crt_output} -keepori"
        else:
            cmd = f"{crisprone}/bin/CrisprAlign -r {args.repeat} -s {fna} -p 0.85 -partial -crt {crt_output} -keepori"
        check_result=subprocess.check_output(cmd, shell=True)

    check_result=subprocess.check_output(f"cat {crt_output} >> {all_crt}", shell=True)
    cmd = f"{python} {spath}/summarize-crispr-new.py -f {crt_output} -spacer {tmp_spacer} -repeat {tmp_repeat} -seqid_prefix {sid}:{tag}:"
    check_result=subprocess.check_output(cmd, shell=True)
    check_result=subprocess.check_output(f"cat {tmp_spacer} >> {all_spacer}", shell=True)
    check_result=subprocess.check_output(f"cat {tmp_repeat} >> {all_repeat}", shell=True)

#post-processing
cmd = f"{cdhit} -i {all_spacer} -c 1.0 -o {args.out}/spacer_graph/spacer_nr100.fa -d 200"
check_result=subprocess.check_output(cmd, shell=True)
cmd = f"{cdhit} -i {all_spacer} -c {args.nr} -o {args.out}/spacer_graph/spacer_nr.fa -d 200"
check_result=subprocess.check_output(cmd, shell=True)
cmd = f"{cdhit} -i {all_repeat} -c {args.nr} -o {args.out}/spacer_graph/repeat_nr.fa -d 200"
check_result=subprocess.check_output(cmd, shell=True)
opath = f"{args.out}/spacer_graph"
cmd = f"{python} {spath}/spacer_sharing.py {opath}/spacer_nr.fa.clstr {opath}/all.crt {opath}/module.all {opath}/module.all.ori"
check_result=subprocess.check_output(cmd, shell=True)
cmd = f"{python} {spath}/extractall-sort.py {opath}/module.all {opath}/module.all.ori"
check_result=subprocess.check_output(cmd, shell=True)
cmd = f"{python} {spath}/spacergraph3.py -i {opath}/module.all.ori -o {opath}/spacer_graph.dot --no-collapse -f 40"
cmd2 = f"{python} {spath}/spacergraph3.py -i {opath}/module.all.ori -o {opath}/spacer_graph_compressed.dot -f 40"
if args.color:
    cmd += f" -c {args.color}"
    cmd2 += f" -c {args.color}"
check_result=subprocess.check_output(cmd, shell=True)
check_result=subprocess.check_output(cmd2, shell=True)

print("All completed")
