python3 ../pipelines/crispr_ann.py --input meta.txt --out crisprone
#make sure --db points to the folder that contains the MGE database
python3 ../pipelines/mge_net.py --input crisprone/ --db /home/team/mgedb --out mgenet
