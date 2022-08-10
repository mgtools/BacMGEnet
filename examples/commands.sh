#run pipeline 1: call CRISPRone to annotate CRISPR-Cas systems
#meta.txt -- the file contains a list of genomes
#--out -- to specify the output directory
python3 ../pipelines/crispr_ann.py --input meta.txt --out crisprone


#run pipelin 2: identy invaders using spacer matching against MGE database and build bacteria-MGE network
#make sure --db points to the folder that contains the MGE database that you downloaded from our website, or your own customized database
#assume you have phage & plasmid database under /home/team/mgedb, you may call the following command directory
#otherwise, modify the --db setting before you run
python3 ../pipelines/mge_net.py --input crisprone/ --db /home/team/mgedb --out mgenet
