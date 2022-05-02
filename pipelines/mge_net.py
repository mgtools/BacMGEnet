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
mgenet = ppath + "/../mgenet/"
crisprone = ppath + "/../CRISPRone/"
crispralign = crisprone + "/bin/CrisprAlign"
blastn = crisprone + "/bin/blast+/bin/blastn"
cdhit = f"{crisprone}/bin/cd-hit-v4.6.1-2012-08-27/cd-hit-est"

print(f"path {ppath} mgenet {mgenet}")
if os.path.exists(mgenet + "/build_gff.py") == False:
    sys.exit("mgenet directory not set correctly")
if os.path.exists(crispralign) == False:
    sys.exit(f"CRISPRAlign {crispralign} not found")

parser = argparse.ArgumentParser()
parser.add_argument("--input", metavar="folder", help="folders containing spacer information; if more than one folder, seperate with spaces", required=True)
parser.add_argument("--tag", metavar="string", help="tags to distinguish multiple inputs")
parser.add_argument("--db", metavar="folder", help="folder where the phage/plasmid sequences are located", required=True)
parser.add_argument("--phage", metavar="file", help="blast database of phage sequences (user-provided)")
parser.add_argument("--plasmid", metavar="file", help="blast database of plasmid sequences (user-provided)")
parser.add_argument("--repeat", metavar="file", help="guide repeats for CRISPR identification in phage/plasmid sequences", required=False)
parser.add_argument("--out", metavar="path", help="specify the path where the outputs will go; use ./ for current folder", required=True)
args = parser.parse_args()

if not os.path.exists(args.out):
    os.mkdir(args.out)

subs = args.input.split()
if len(subs) == 1:
    args.spacer = f"{subs[0]}/spacer_graph/spacer_nr100.fa"
    args.array = f"{subs[0]}/spacer_graph/module.all.ori"
else:
    #merge multiple spacers
    print("merge multiple spacer inputs..")
    all_spacer = f"{args.out}/spacer_combined.fa"
    args.spacer = f"{args.out}/spacer_combined_nr100.fa"
    subprocess.check_output("rm -f " + all_spacer, shell=True)
    if args.tag:
        tags = args.tag.split()
    else:
        for idx in range(len(subs)):
            tags.append("inp" + str(idx + 1))
    args.array = ""
    for idx in range(len(subs)):
        thisspacer = f"{subs[idx]}/spacer_graph/spacer_all.fa"
        thisarray = f"{subs[idx]}/spacer_graph/module.all.ori"
        check_result=subprocess.check_output(f"cat {thisspacer} >> {all_spacer}", shell=True)
        if not args.array:
            args.array = tags[idx] + ":" + thisarray  
        else:
            args.array = args.array + " " + tags[idx] + ":" + thisarray  
    #100%-nr 
    cmd = f"{cdhit} -i {all_spacer} -c 1.0 -o {args.spacer} -d 200"
    check_result=subprocess.check_output(cmd, shell=True)
    print("merged spacers saved to file", args.spacer)

    
#search spacer sequences against phage & plasmid databases
phagedb = f"{args.db}/all_phage_seq.fna"
plasmiddb = f"{args.db}/all_plasmid_seq.fna"
m8phage = f"{args.out}/spacer_vs_phages.m8"
m8plasmid = f"{args.out}/spacer_vs_plasmid.m8"
command = f'{blastn} -query {args.spacer} -db {phagedb} -evalue 0.001 -outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"  -num_threads 8 -out {m8phage}'
check_result=subprocess.check_output(command, shell=True)
command = f'{blastn} -query {args.spacer} -db {plasmiddb} -evalue 0.001 -outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"  -num_threads 8 -out {m8plasmid}'
check_result=subprocess.check_output(command, shell=True)

#run CRISPRAlign and filter the similarity search results if repeat is provided
m8phage_sel = m8phage
m8plasmid_sel = m8plasmid
if args.repeat:
    #run crispralign
    command = f"{crispralign} -r {args.repeat} -s {plasmiddb} -partial -d 5 -crt {out}/all_plasmid_seq.db.crt"
    check_result = subprocess.check_output(command, shell=True)
    command = f"{crispralign} -r {args.repeat} -s {phagedb} -partial -d 5 -crt {out}/all_phage_seq.db.crt"
    check_result = subprocess.check_output(command, shell=True)
    #filter
    print(f"Filtering Blast Results: {m8phage}")
    cmd = f"{python} {mgenet}/filter_blast_w_crt.py --phage_crt ${out}/all_phage_seq.db.crt -b ${m8phage}"
    check_result = subprocess.check_output(command, shell=True)
    print(f"Filtering Blast Results: {m8plasmid}")
    cmd = f"{python} {mgenet}/filter_blast_w_crt.py --plasmid_crt ${out}/all_plasmid_seq.db.crt -b ${m8plasmid}"
    check_result = subprocess.check_output(command, shell=True)

    m8phage_sel = "f{out}/spacer_nr_phage.filtered.m8"
    m8plasmid_sel = "f{out}/spacer_nr_plasmids.filtered.m8"

#build gff file
print("generate gff file..")
cmd = f"{python} {mgenet}/build_gff.py --phage {m8phage_sel} --plasmid {m8plasmid_sel} --phagedb {phagedb} --plasmiddb {plasmiddb} --outdir {args.out}"
check_result = subprocess.check_output(cmd, shell=True)
print(f"  result saved to {args.out}/protospacer_w_multihits.gff")

#greedy algorithm
print("apply greedy algorithm to select..")
cmd = f"{python} {mgenet}/mge_nr_greedy.py --gff {args.out}/protospacer_w_multihits.gff --out {args.out}/protospacer_w_multihits-greedy.gff"
check_result = subprocess.check_output(cmd, shell=True)
print(f"  result saved to {args.out}/protospacer_w_multihits-greedy.gff")

#create network file
print("create network..")
cmd = f'{python} {mgenet}/spacer2mge_network_adv.py --gff {args.out}/protospacer_w_multihits-greedy.gff --cdhit {args.spacer}.clstr --ori "{args.array}" --spacer2mge {args.out}/spacer2mge.gml --host2mge {args.out}/host2mge.gml'
check_result = subprocess.check_output(cmd, shell=True)

print("All completed")
