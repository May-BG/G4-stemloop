# G4-stemloop
The repository contains all the code required to run the stem loop test of G-quadruplex and generate results.
## 1 map G4 across seven primate species' T2T genomes


## 2 identify ape shared G4s and build substitution spectrum
### To identify substitutions in chimp:
for i in {1..22} ; do Rscript chimpid_snp_count_sharedG4_noanc.R $i; done \
Rscript chimpid_snp_count_sharedG4_noanc.R "X" \
Rscript chimpid_snp_count_sharedG4_noanc.R "Y" \
Rscript shared_substitution_chimp_hsachr.R \
sbatch intersect_infer_chimp_hsachr.sh \


sbatch intersect_g4_stemloop_subsets.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/human_lost/chimp_id" \
sbatch intersect_g4_stemloop_subsets "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/chimp_specific/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/chimp_specific/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/shared_g4/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/disappear_g4_in_orang/chimp_id" \
sbatch intersect_g4_stemloop_subsets_human.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/disappear_g4_in_orang" \
sbatch intersect_g4_stemloop_subsets_human.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/shared_g4" \
sbatch intersect_g4_stemloop_subsets_human.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/unique_g4" \
sbatch intersect_g4_stemloop_subsets_human.sh "/storage/home/xmz5176/newgroup/data/G4/T2T/autosomal/map_g4_human_complete/disappear_g4" \

## 3 HKA test and visualization

## 4 PhyloFit test
