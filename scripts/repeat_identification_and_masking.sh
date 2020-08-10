# Following Avril Coghlan's blogpost (http://avrilomics.blogspot.com/2015/09/ltrharvest.html)

# Download HMM profiles
wget http://gydb.org/gydbModules/collection/collection/db/GyDB_collection.zip
# Unzip program is found in flo environment
sconda /ceph/users/pgonzales/.conda/envs/flo
unzip GyDB_collection.zip
sconda /ceph/users/pgonzales/miniconda3/envs/LTR_retriever
mkdir pfam
# Pfams from tables B1 and B2 from the LTRdigest paper by Steinbiss et al (2009) doi: 10.1093/nar/gkp759
for pfam in PF03732 PF07253 PF00077 PF08284 PF00078 PF07727 PF06817 PF06815 PF00075 PF00552 PF02022 PF00665 PF00098 PF00385 PF01393 PF00692 PF01021 PF03078 PF04094 PF08330 PF04195 PF05380 PF01140 PF02337 PF01141 PF00607 PF02093 PF00692 PF00077 PF00078 PF06817 PF00552 PF02022 PF00665 PF00075 PF00429 PF08791 PF09590 PF03408 PF00516 PF03056 PF00517 PF00098; do
wget http://pfam.xfam.org/family/$pfam/hmm -O pfam/${pfam}.txt
hmmconvert -2 pfam/${pfam}.txt > pfam/${pfam}.hmm
done
# These profiles were empty
rm -f pfam/PF08330.hmm pfam/PF08791.hmm pfam/PF03056.hmm
# COncatenate hmms
mkdir pfam_and_GyDB_profiles
ln -s $PWD/pfam/*.hmm pfam_and_GyDB_profiles
ln -s $PWD/GyDB_collection/profiles/*hmm pfam_and_GyDB_profiles

# LTR harvest
$gt suffixerator -db $assembly -indexname $assemName -tis -suf -lcp -des -ssp -sds -dna
$gt ltrharvest -index $assemName -gff3 ${assemName}.ltrharvest.gff  -tabout no > ${assemName}.ltrharvest.out

# LTR digest
$gt gff3 -sort ${assemName}.ltrharvest.gff > ${assemName}.ltrharvest.sorted.gff
$gt ltrdigest -hmms pfam_and_GyDB_profiles/*hmm -outfileprefix ${assemName}.ltrharvest ${assemName}.ltrharvest.sorted.gff $assemName > ${assemName}.ltrdigest.gff

# Filter predictions that have a hit to one of the TE related domains
$gt select -rule_files /data/blaxter/Otipulae_lr/scripts/filter_protein_match.lua -- < ${assemName}.ltrdigest.gff > ${assemName}.ltrdigest.filteredHit.gff

# Rename the sequences
grep ">" $assembly | sed 's/>//' > tmp.txt
awk -F '\t' 'BEGIN{OFS=FS; count = 0}{print $1, "seq"count; count++}' tmp.txt > rnm.txt
seqkit replace -p '(.+)$' -r '{kv}' -k rnm.txt $assembly > assm.rnm.fa

# Filter to only get coordinates of LTR_retrotransposon feature
awk -F '\t' 'BEGIN{OFS=FS} ($3 ~ /LTR_retrotransposon/){print}' ${assemName}.ltrdigest.filteredHit.gff > ltr.gff
# Extract the sequences
bedtools getfasta -name -fi assm.rnm.fa -bed ltr.gff > ${assemName}.ltrdigest.filteredHit.fa

# Classify the repeats
/ceph/software/repeatmodeler/RepeatModeler.v2.0.1/RepeatClassifier -consensi ${assemName}.ltrdigest.filteredHit.fa > ${assemName}.ltrdigest.filteredHit.classified.fa



# Repeatmodeler
bash /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/repeatmodeler2/oscheius_tipulae.local.v3_3.sh
/ceph/software/repeatmodeler/RepeatModeler.v2.0.1/RepeatClassifier -consensi oscheius_tipulae.local.v3_3.genomic-families.fa

# TransposonPSI
/ceph/software/repeatmodeler/RepeatModeler.v2.0.1/RepeatClassifier -consensi /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/transposonPSI/oscheius_tipulae.local.v3_3.genomic.fa.TPSI.allHits.chains.bestPerLocus.fa 


# USEARCH
# http://avrilomics.blogspot.com/2015/09/creating-merged-repeatmodeler-and.html
cd /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/usearch
ln -s /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/repeatmodeler2/oscheius_tipulae.local.v3_3.genomic-families.fa oscheius_tipulae.local.v3_3.repeatmodeler2.lib
ln -s /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/ltr_finder/oscheius_tipulae.local.v3_3.ltrdigest.filteredHit.fa.classified oscheius_tipulae.local.v3_3.fa.classified.ltrdigest.lib
# Filter transposonPSI sequences shorter than 50 nt
seqkit seq -m 50 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/transposonPSI/oscheius_tipulae.local.v3_3.genomic.fa.TPSI.allHits.chains.bestPerLocus.fa.classified > oscheius_tipulae.local.v3_3.fa.classified.TPSI.lib


# Reformat the fasta files to have one sequence per line
for lib in $(ls *lib); do
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $lib | tail -n +2 > ${lib}2
done

cat *lib2 > merged.fa
usearch7.0.1090_i86linux64 -sortbylength merged.fa --output merged.sorted.fa --log usearch.log
usearch7.0.1090_i86linux64 -cluster_fast merged.sorted.fa --id 0.8 --centroids my_centroids.fa --uc result.uc -consout final.nr.consensus.fa -msaout aligned.fasta --log usearch2.log

# Clean fasta headers
sed 's/centroid=//; s/;seqs=.\+//' final.nr.consensus.fa > final.nr.consensus.fa2
# Check that all headers are shorter than 80 characters
grep ">" final.nr.consensus.fa2 | awk '{print length($0)}' | sort -nr | less

# RepeatMasker
bash /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/repeatmasker/rmask_oscheius_tipulae.local.v3_3.sh
