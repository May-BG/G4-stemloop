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


sbatch intersect_g4_stemloop_subsets.sh "/human_lost/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/chimp_specific/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/chimp_specific/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/shared_g4/chimp_id" \
sbatch intersect_g4_stemloop_subsets.sh "/disappear_g4_in_orang/chimp_id" \
sbatch intersect_g4_stemloop_subsets_human.sh "/disappear_g4_in_orang" \
sbatch intersect_g4_stemloop_subsets_human.sh "/shared_g4" \
sbatch intersect_g4_stemloop_subsets_human.sh "/unique_g4" \
sbatch intersect_g4_stemloop_subsets_human.sh "/disappear_g4" \

Rscript compare_sequence_divergence.R

## 3 HKA test and visualization
Rscript HKA_test_fig3_4.R

## 4 PhyloFit test
