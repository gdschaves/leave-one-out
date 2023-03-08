module load  gcc/6.2.0
module load  xz/5.2.2 
module load  bzip2/1.0.6
module spider bedtools/2.29.0

cd /gpfs/data/applebaum-lab/Gepoliano/leave-one-out/data

bedtools intersect -wa -a /gpfs/data/applebaum-lab/Mohan_sandbox/gencode.v33.annotation.protein_only.gtf \
-b /gpfs/data/applebaum-lab/Gepoliano/leave-one-out/data/HIF1a_N_BE2_merged.bed > \
/gpfs/data/applebaum-lab/Gepoliano/leave-one-out/results/HIF1a_N_BE2_genes.txt