mask otipu



Functional enrichment - annotate protein domains - kinfin enrichment
separate intron from exon density





grep -P '^\d+_t' /data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/all_longestProteins.faa.tsv | cut -f 1,14 | sed 's/|/, /g' > /data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/oti3.3_longestProteins.faa.topGo.tsv








# einvert parsing
threads=1
condaEnv=/ceph/users/pgonzales/.conda/envs/EMBLmyGFF3
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/einvert/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=einverted_oscheius_tipulae.local.v3_3


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; einverted -gap 12 -threshold 50 -match 3 -mismatch -4 -maxrepeat 2000 -sequence $assembly -outfile ${assemName}.txt -outseq  ${assemName}.fasta; grep -oP \"^\\w+: S\" ${assemName}.txt | sed 's/: S//; s/_/\\t/g' > ${assemName}.tsv"

assemName=testCmd
einverted -gap 12 -threshold 50 -match 3 -mismatch -4 -maxrepeat 2000 -sequence $assembly -outfile ${assemName}.txt -outseq  ${assemName}.fasta
einverted
# grep -oP "\w+: S|^\s*\d+" ${assemName}.txt | sed 's/: S//; s/^\s\+//' | awk 'NR%3{printf "%s\t",$0;next;}1' > ${assemName}.tsv

grep -oP "\w+: S|\d+" oscheius_tipulae.local.v3_3 | sed 's/: S//' | awk 'NR%10{printf "%s\t",$0;next;}1' | awk -F '\t' 'BEGIN{OFS=FS}{print $1, $7, $9, $8-$7, $10-$8, $2}' > oscheius_tipulae.local.v3_3.tsv

grep -oP "\w+: S|\d+" oscheius_tipulae.local.v3_3 | sed 's/: S//' | awk 'NR%10{printf "%s\t",$0;next;}1' | awk -F '\t' 'BEGIN{OFS=FS}{print $1, $7, $8, $2"\n"$1, $10, $9, $2}' | sort -k1,1 -k2,2n > oscheius_tipulae.local.v3_3.bed
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b oscheius_tipulae.local.v3_3.bed > einvert_oscheius_tipulae.local.v3_3.1kw.bed


I: S
2	437 score
3	591	matches
4	895	span
5	66	percent
6	10	gaps
7	2	s1
8	902	e1
9	1802	e2
10	904	s2


threads=1
condaEnv=/ceph/users/pgonzales/.conda/envs/EMBLmyGFF3
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/etandem/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=etandem_oscheius_tipulae.local.v3_3

for chrom in I II III IV V X; do


qsub -cwd -o ${outdir}/${chrom}_${assemName}.log -e ${outdir}/${chrom}_${assemName}.err -N ${chrom}_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; samtools faidx $assembly $chrom > ${chrom}.fa; etandem -sequence ${chrom}.fa -minrepeat 4 -maxrepeat 300 -outfile ${chrom}.etandem -origfile ${chrom}.orig_etandem; sed 's/\\s\\+/\\t/g' ${chrom}.orig_etandem | sed 's/^\\t//' | awk -v chr=$chrom -F '\\t' 'BEGIN{OFS=FS}{print chr, \$0}' > ${chrom}.orig_etandem.tsv"

done

# sed 's/\s\+/\t/g' oscheius_tipulae.local.v3_3.orig_etandem | sed 's/^\t//' | awk -v chr=$chrom -F '\t' 'BEGIN{OFS=FS}{print chr, $0}' > oscheius_tipulae.local.v3_3.orig_etandem.tsv
cat *orig_etandem.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $1, $3, $4, $2}' | sort -k1,1 -k2,2n | bedtools merge -i stdin > etandem_oscheius_tipulae.local.v3_3.orig_etandem.bed
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b etandem_oscheius_tipulae.local.v3_3.orig_etandem.bed > etandem_oscheius_tipulae.local.v3_3.orig_etandem.1kw.bed



# IRF parsing
~/irf305.linux.exe oscheius_tipulae.local.v3_3.fa 3 4 12 80 10 50 30000 300000 -d -h
awk 'BEGIN{FS=" ";OFS="\t"}{if($1 ~ /Sequence/){chrom=$2}; if($1 ~ /^[0-9]/){print chrom, $1, $2, $10"\n"chrom, $4, $5, $10}}' analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.fa.3.4.12.80.10.50.30000.300000.dat | sort -k1,1 -k2,2n > analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.bed

bedtools merge -i analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.bed > analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.merged.bed
bedtools makewindows -g analyses/assemblies/oscheius_tipulae.local.v3_3.genomic.fa.fai -w 1000 | bedtools intersect -a stdin -b analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.merged.bed | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "irf"}'> analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.merged.1kw.bed


# TRF parsing
trf oscheius_tipulae.local.v3_3.genomic.fa 2 7 7 80 10 50 500 -f -d -m
awk 'BEGIN{FS=" ";OFS="\t"}{if($1 ~ /Sequence/){chrom=$2}; if($1 ~ /^[0-9]/){print chrom, $1, $2}}' oscheius_tipulae.local.v3_3.genomic.fa.2.7.7.80.10.50.500.dat | sort -k1,1 -k2,2n | bedtools merge -i stdin >trf_oscheius_tipulae.local.v3_3.bed
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b trf_oscheius_tipulae.local.v3_3.bed | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "trf"}' > trf_oscheius_tipulae.local.v3_3.1kw.bed

# Gene density
sort -k1,1 -k2,2n /data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/dict/oscheius_tipulae.local.v3_3.tsv | bedtools merge -i stdin > oscheius_tipulae.local.v3_3.transcriptCoords.merged.tsv
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b oscheius_tipulae.local.v3_3.transcriptCoords.merged.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "protein_coding"}' > gene_${assemName}.1kw.tsv

# Exon density
zcat /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz | \
 awk -F '\t' 'BEGIN{OFS=FS}($3 == "exon"){print}' | sort -k1,1 -k4,4n | bedtools merge -i stdin > oscheius_tipulae.local.v3_3.codingExons.merged.tsv
 
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b oscheius_tipulae.local.v3_3.codingExons.merged.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "coding_exons"}' > exon_${assemName}.1kw.tsv


# Intron density
zcat /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz | \
 awk -F '\t' 'BEGIN{OFS=FS}($3 == "mRNA"){print}' | sort -k1,1 -k4,4n | bedtools merge -i stdin > ${assemName}.codingGene.merged.tsv
bedtools subtract -a ${assemName}.codingGene.merged.tsv -b ${assemName}.codingExons.merged.tsv > ${assemName}.codingIntrons.merged.tsv
 
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b ${assemName}.codingIntrons.merged.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "coding_introns"}' > intron_${assemName}.1kw.tsv


# C elegans ortholog density
grep Transcript /data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/dict/caenorhabditis_elegans.PRJNA13758.WBPS14.tsv > caenorhabditis_elegans.PRJNA13758.WBPS14.transcriptCoords.tsv

R
library(tidyverse)

rbbh <- read_tsv("/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/rbbh/genemark_vs_celeg.txt",
                 col_names = c("A", "B", "MeanBitscore", "A_bitscore", "B_bitscore", "A_eval", "B_eval", "A_qcov", "B_qcov", "A_pid", "B_pid"),
                 comment = "#") %>%
  mutate(celeg = ifelse(grepl("Transcript", A), A, B),
         otipu = ifelse(grepl("Transcript", A), B, A)) %>%
  select(celeg, otipu)

ce_pos <- read_tsv("caenorhabditis_elegans.PRJNA13758.WBPS14.transcriptCoords.tsv",
                   col_names = c("chr", "start", "end", "strand", "tID"))
ot_gff <- read_tsv("/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz",
                   col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "tID")) %>%
                   mutate(tID = sub(".+=(\\d+_t).*", "\\1", tID)) %>%
                   filter(feature == "exon")

inner_join(ot_gff, rbbh, by = c("tID" = "otipu")) %>%
  mutate(ce_chr = ce_pos$chr[match(celeg, ce_pos$tID)]) %>%
  select(chr, start, end, ce_chr) %>%
  write_tsv("oscheius_tipulae.local.v3_3.ceChr.tsv",
            col_names = F)

exit()

sort -k1,1 -k2,2n ${assemName}.ceChr.tsv | bedtools merge -i stdin -c 4 -o distinct > ${assemName}.ceChr.merged.tsv
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -wb -a stdin -b ${assemName}.ceChr.merged.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $1, $2, $3, $7}' > ceChr_${assemName}.1kw.tsv


# Nigon density

R

library(tidyverse)

ot_nigons <- read_tsv("/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/gene2Nigon20200528_09_55.tsv.gz") %>%
  filter(assembly == "oscheius_tipulae.local.v3_3",
         kinfCut == "0.84") %>%
  select(Orthogroup, nigon, tID)

ot_pos <- read_tsv("/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz",
                   col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "tID")) %>%
                   mutate(tID = sub(".+=(\\d+_t).*", "\\1", tID)) %>%
                   filter(feature == "exon")

inner_join(ot_nigons, ot_pos, by = "tID") %>%
  mutate(nigon = paste0("Ng", nigon)) %>%
  select(chr, start, end, nigon) %>%
  write_tsv("oscheius_tipulae.local.v3_3.nigonCoords.tsv",
            col_names = F)

exit()

sort -k1,1 -k2,2n ${assemName}.nigonCoords.tsv | bedtools merge -i stdin -c 4 -o distinct > ${assemName}.nigonCoords.merged.tsv
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -wb -a stdin -b ${assemName}.nigonCoords.merged.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $1, $2, $3, $7}' > nigon_${assemName}.1kw.tsv


# RepeatMasked densities
sort -k1,1 -k4,4n /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/repeatmasker/oscheius_tipulae.local.v3_3.genomic.fa.out.gff | bedtools merge -i stdin > rmasker_${assemName}.tsv
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b rmasker_${assemName}.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "repeatmasker"}' > rmasker_${assemName}.1kw.tsv
bedtools genomecov -i analyses/annotation/repeats/rmasker_oscheius_tipulae.local.v3_3.1kw.tsv -g analyses/assemblies/oscheius_tipulae.local.v3_3.genomic.fa.fai

# telomere parsing
seqkit locate -p AAGCCT ${assembly} | awk -F '\t' 'BEGIN{OFS=FS}NR>1{print $1, $5, $6}' > ${assemName}.teloRepeats.tsv
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -a stdin -b ${assemName}.teloRepeats.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $0, "telomeres"}' > telomere_${assemName}.1kw.tsv

# In MacBook
for feat in trf gene ceChr rmasker telomere; do
	scp qmaster:~/${feat}_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/
done

cat analyses/annotation/repeats/nigon_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/ceChr_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/rmasker_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/telomere_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/gene_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/trf_oscheius_tipulae.local.v3_3.1kw.bed analyses/annotation/repeats/irf_oscheius_tipulae.local.v3_3.merged.1kw.bed analyses/annotation/repeats/intron_oscheius_tipulae.local.v3_3.1kw.tsv analyses/annotation/repeats/exon_oscheius_tipulae.local.v3_3.1kw.tsv | gzip -c > analyses/annotation/repeats/concat_oscheius_tipulae.local.v3_3.1kw.tsv.gz

###
kat comp
###

threads=10
condaEnv=/ceph/users/pgonzales/miniconda3/envs/kat
outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/kat
pe_reads1=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_1.fastq.gz
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=oscheius_tipulae.local.v3_3


qsub -cwd -o ${outdir}/kat_${assemName}.log -e ${outdir}/kat_${assemName}.err -N kat_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; kat comp -o $assemName -t $threads $pe_reads1 $assembly"

threads=10
condaEnv=/ceph/users/pgonzales/miniconda3/envs/kat
outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/kat
pe_reads1=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_1.fastq.gz
assembly=/data/blaxter/Otipulae_lr/raw/genome/oscheius_tipulae.PRJEB15512.WBPS14.genomic.fa
assemName=oscheius_tipulae.PRJEB15512.WBPS14


qsub -cwd -o ${outdir}/kat_${assemName}.log -e ${outdir}/kat_${assemName}.err -N kat_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; kat comp -o $assemName -t $threads $pe_reads1 $assembly"


###
Rename scaffolds
###
# In MacBook
seqkit replace -p '(\S+) .+$' -r '{kv}' -k analyses/quast/otipu3.3/dict_ot2.0.tsv analyses/assemblies/oscheius_tipulae.PRJEB15512.WBPS14.genomic.fa > analyses/assemblies/oscheius_tipulae.PRJEB15512.WBPS14.genomic.rnm.fa
nucmer -p analyses/dotplot/otipulae/not2nm_vs_ot3.3 analyses/assemblies/freebayes3_otipulaeA_rRNAsized.fasta analyses/assemblies/oscheius_tipulae.PRJEB15512.WBPS14.genomic.rnm.fa
sconda blobtools2
python ~/bin/DotPrep.py --delta analyses/dotplot/otipulae/not2nm_vs_ot3.3.delta --out analyses/dotplot/otipulae/not2nm_vs_ot3.3 --unique-length 500 --overview 999999

assembly=analyses/assemblies/oscheius_tipulae.local.v3_3.genomic.fa
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools nuc -fi ${assembly} -bed stdin > analyses/curation/oscheius_tipulae.local.v3_3.gc.bed






####
Call variants
####


threads=1
condaEnv=/ceph/users/pgonzales/.conda/envs/legacyBlast
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/transposonPSI/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=oscheius_tipulae.local.v3_3.genomic


qsub -cwd -o ${outdir}/psi_${assemName}.log -e ${outdir}/psi_${assemName}.err -N psi_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; perl /data/blaxter/Otipulae_lr/scripts/TransposonPSI_08222010/transposonPSI.pl $assembly nuc"


nextflow -C /data/blaxter/Otipulae_lr/scripts/rbbh.config run /data/blaxter/Otipulae_lr/scripts/rbbh.nf -resume --query "/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/braker/abinitio/braker/augustus.ab_initio.aa" --ref "/data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/caenorhabditis_elegans.PRJNA13758.WBPS14.faa" --out /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/rbbh/braker.abinitio_vs_celeg.txt --chunkSize 500 --minqueryCov 20


# extract proteins
source /ceph/users/pgonzales/miniconda3/bin/activate /ceph/users/pgonzales/.conda/envs/gffmagic
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/agat

annotDir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/
comprAnnotDir=$annotDir/comparedAnnots/raw_link


16,411 # gm_Ce
16,405 # gm_uniP



threads=1
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/agat_gff3
condaEnv=/ceph/users/pgonzales/.conda/envs/gffmagic


for annot in /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/raw_link/augustus_busco.gff /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/raw_link/augustus_ce.gff /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/raw_link/braker_abinitio.gtf /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/raw_link/funannotate_Ce.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/raw_link/funannotate_gm_uniP.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/raw_link/genemark.gtf; do
# Prefix
prFix=$(basename $annot .gff)
prFix2=${prFix%.gtf}
assemName=${prFix2%.gff3}

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/exe_agat.sh; cat ${outdir}/exe_agat.sh; agat_convert_sp_gxf2gxf.pl -g $annot -o $outdir/${assemName}.gff3"

done



protsDir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3
outDir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3
wkdir=/scratch/pgonzales/orthofinder/
threads=56


qsub -cwd -o ${outDir}/orthofinder_t.log -e ${outDir}/orthofinder_t.err -N orthofinder -V -pe smp $threads -b yes "source /ceph/software/conda/bin/activate /ceph/users/pgonzales/.conda/envs/orthofinder; mkdir -p $wkdir; echo localhost: \$HOSTNAME; echo "orthofinder version:" \$(orthofinder --version | head -n 2 | tail -n +2); echo outDir: $outDir; cd $outDir; echo cmnd: ${outDir}/exec_orthofinder.sh; cat ${outDir}/exec_orthofinder.sh; echo date \$(date); orthofinder -t $threads -a $threads -M msa -S blast"




threads=1
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
condaEnv=/ceph/users/pgonzales/.conda/envs/gffmagic

for annot in /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/augustus_busco.longest.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/augustus_ce.longest.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/braker_abinitio.longest.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/funannotate_Ce.longest.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/funannotate_gm_uniP.longest.gff3 /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/genemark.longest.gff3; do

assemName=$(basename $annot .gff3)

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/exe_agat.sh; cat ${outdir}/exe_agat.sh; agat_sp_extract_sequences.pl -g $annot -f $assembly -p --clean_final_stop --clean_internal_stop -o $outdir/${assemName}.faa"

done


threads=12
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/busco
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
condaEnv=/ceph/users/pgonzales/.conda/envs/gffmagic

for gprot in /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/augustus_busco.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/augustus_ce.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/braker_abinitio.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/funannotate_Ce.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/funannotate_gm_uniP.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/genemark.longest.faa; do

assemName=$(basename $gprot .longest.faa)

species=ascaris_suum.local.v2020
scriptFile=$outdir/exe_busco.sh

qsub -cwd -o ${outdir}/run_${assemName}.log -e ${outdir}/run_${assemName}.err -N pbusco_$assemName -V -pe smp $threads -b yes "source /ceph/software/conda/bin/activate /ceph/users/pgonzales/.conda/envs/buscoOnly; cd $outdir; echo localhost: \$HOSTNAME; echo \$(busco --version); echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: $scriptFile; cat $scriptFile; echo date \$(date);  export AUGUSTUS_CONFIG_PATH=/ceph/software/busco_augustus_config_path/config/; cd $outdir; busco --offline -f -c $threads -l /data/blaxter/Otipulae_lr/analyses/qualityMetrics/busco/assemblies/table4/busco_downloads/lineages/nematoda_odb10 -i $gprot --out $assemName --mode prot"

done



outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/rbbh
refProt=/data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/caenorhabditis_elegans.PRJNA13758.WBPS14.faa
cd $outdir

for gprot in /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/augustus_busco.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/augustus_ce.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/braker_abinitio.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/funannotate_Ce.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/funannotate_gm_uniP.longest.faa /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/genemark.longest.faa; do

assemName=$(basename $gprot .longest.faa)

nextflow -C /data/blaxter/Otipulae_lr/scripts/rbbh.config run /data/blaxter/Otipulae_lr/scripts/rbbh.nf --query $gprot --ref $refProt --out $outdir/${assemName}_vs_celeg.txt --chunkSize 500 --minqueryCov 20

rm -rf work

done





cd /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/rbbh
for proteome in /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/braker/abinitio/braker/augustus.ab_initio.aa 


# genemark
threads=6
condaEnv=/ceph/users/pgonzales/.conda/envs/lc_braker2
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/genemark
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=genemark_$(basename $assembly .fa)





qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; perl /data/blaxter/Otipulae_lr/scripts/gmes_linux_64/gmes_petap.pl --cores $threads --ES --sequence $assembly > ${assemName}.gff"



threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/funannotate/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=funannotate_gm_$(basename $assembly .fa)
transcripts=/data/blaxter/Otipulae_lr/analyses/RNA-assembly/mira/otipulae_ERR1880912.fasta
brakerGff=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/agat_gff3/braker_abinitio.gff3

cmd="export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/; export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; funannotate predict --busco_db nematoda --busco_seed_species caenorhabditis --GENEMARK_PATH /data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --cpus $threads --force -i $assembly -o $outdir/$assemName -s \"Oscheius tipulae\" --name OT3 --other_gff ${brakerGff}:2 --weights Augustus:0 snap:1 GlimmerHMM:1 GeneMark:3 --organism other --transcript_evidence $transcripts"


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"





threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/funannotate
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=funannotate_gm_$(basename $assembly .fa)
transcripts=/data/blaxter/Otipulae_lr/analyses/RNA-assembly/mira/otipulae_ERR1880912.fasta
proteins=/data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/caenorhabditis_elegans.PRJNA13758.WBPS14.faa

cmd="export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/; export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; funannotate predict --busco_db nematoda --busco_seed_species caenorhabditis --GENEMARK_PATH /data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --cpus $threads --force -i $assembly -o $outdir/$assemName -s \"Oscheius tipulae\" --name OT3 --augustus_species caenorhabditis --organism other --transcript_evidence $transcripts --protein_evidence $proteins"



qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"


threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/funannotate
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=funannotate_gm_uniP_$(basename $assembly .fa)
transcripts=/data/blaxter/Otipulae_lr/analyses/RNA-assembly/mira/otipulae_ERR1880912.fasta

cmd="export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/; export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; funannotate predict --busco_db nematoda --busco_seed_species caenorhabditis --GENEMARK_PATH /data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --cpus $threads --force -i $assembly -o $outdir/$assemName -s \"Oscheius tipulae\" --name OT3 --augustus_species caenorhabditis --organism other --transcript_evidence $transcripts"



qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"

# Braker protein hints
threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/lc_braker2
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/braker/protHints
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=brakerPHints_$(basename $assembly .fa)
proteins=/data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/caenorhabditis_elegans.PRJNA13758.WBPS14.faa

cmd="export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; braker.pl --GENEMARK_PATH=/data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --genome=$assembly --prot_seq=$proteins --prg=gth --gth2traingenes --trainFromGth --cores=$threads --crf"

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"

# Braker 
threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/lc_braker2
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/braker/abinitio
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=braker_$(basename $assembly .fa)

cmd="export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; braker.pl --GENEMARK_PATH=/data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --genome=$assembly --esmode --cores=$threads --crf"

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"

# funannotate with braker
threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/funannotate/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=funannotate_gm_braker_uniP_$(basename $assembly .fa)
transcripts=/data/blaxter/Otipulae_lr/analyses/RNA-assembly/mira/otipulae_ERR1880912.fasta
auguGff=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/braker/abinitio/braker/augustus.ab_initio.gtf

cmd="export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/; export PERL5LIB=/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/ceph/users/pgonzales/.conda/envs/lc_braker2/lib/site_perl/5.26.2; funannotate predict --busco_db nematoda --busco_seed_species caenorhabditis --GENEMARK_PATH /data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --cpus $threads --force -i $assembly -o $outdir/$assemName -s \"Oscheius tipulae\" --name OT3 --augustus_gff $auguGff --augustus_species caenorhabditis --organism other --transcript_evidence $transcripts"




# Quast compare
/ceph/users/pgonzales/.conda/envs/quast/

# extract proteins with agat
threads=1
condaEnv=/ceph/software/conda/envs/kinfin_env/
outdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/kinfin
orthoFdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/WorkingDirectory/OrthoFinder/Results_May27/
targetFrac=0.84
assemName=kinfin_$targetFrac

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; sed 's/.faa//' $orthoFdir/WorkingDirectory/SpeciesIDs.txt | awk -F \": \" 'BEGIN{print\"#IDX,TAXON\"; OFS = \",\"}{print \$1,\$2}' > config.txt; /ceph/software/kinfin/kinfin/kinfin -c config.txt -g $orthoFdir/Orthogroups/Orthogroups.txt -s $orthoFdir/WorkingDirectory/SequenceIDs.txt -p $orthoFdir/WorkingDirectory/SpeciesIDs.txt --target_fraction $targetFrac -o $outdir/kinfin_$targetFrac"




# repeats
# LTR harvest
cd /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/ltr_finder
gt=/ceph/software/genometools/genometools-1.5.9/bin/gt
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=oscheius_tipulae.local.v3_3


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


# Discard genes that overlap repeats
# Get mRNA coords
zcat /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/oscheius_tipulae.local.v3_3.gff3.gz | awk -F '\t' 'BEGIN{OFS=FS} ($3 ~ /gene/){print}' > protGenes.tsv
bedtools intersect -f 0.4 -wa -u -a protGenes.tsv -b oscheius_tipulae.local.v3_3.genomic.fa.out.gff > re_protGenes.tsv
zcat /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/oscheius_tipulae.local.v3_3.gff3.gz | bedtools subtract -a stdin -b re_protGenes.tsv | gzip -c > /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz

zcat /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz | awk -F '\t' 'BEGIN{OFS=FS} ($3 ~ /mRNA/){print $1, $4, $5, $9}' | sed 's/;Parent.\+//; s/ID=//' | gzip -c > /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.protCodGenes.tsv.gz

curl -T /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/noRe_oscheius_tipulae.local.v3_3.gff3.gz ftp://ftp.ed.ac.uk/edupload/ --user anonymous:a@b.c

curl -T re_protGenes.tsv ftp://ftp.ed.ac.uk/edupload/re_oscheius_tipulae.local.v3_3.tsv --user anonymous:a@b.c


# Parse RM 
perl -I /ceph/software/perl-ensembl/BioPerl-1.6.924/ /data/blaxter/Otipulae_lr/scripts/parseRM.pl -p -i oscheius_tipulae.local.v3_3.genomic.fa.out

# RM out to R

# inverted repeats
cd /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/einvert/

sconda /ceph/users/pgonzales/.conda/envs/EMBLmyGFF3
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=oscheius_tipulae.local.v3_3

ln -s $assembly $(basename $assembly)
einverted -gap 12 -threshold 50 -match 3 -mismatch -4 -sequence $(basename $assembly) -outfile ${assemName}.einverted -outseq ${assemName}.einverted.fa

# tandem repeats
cd /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/etandem/

sconda /ceph/users/pgonzales/.conda/envs/EMBLmyGFF3
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=oscheius_tipulae.local.v3_3

ln -s $assembly $(basename $assembly)
etandem -sequence $(basename $assembly) -outfile ${assemName}.etandem -origfile ${assemName}.orig_etandem -minrepeat 10 -maxrepeat 10



# Annotate proteins

# "Install" interproscan5
cd /data/blaxter/Otipulae_lr/raw/funannotate_db/ 
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.44-79.0/interproscan-5.44-79.0-64-bit.tar.gz
tar -zxf interproscan-5.44-79.0-64-bit.tar.gz


# test dataset

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ceph/users/pgonzales/miniconda3/envs/kat/lib/
$interpro -i tmp.faa -cpu 40 -f xml

download_eggnog_data.py --data_dir /data/blaxter/Otipulae_lr/raw/funannotate_db/

export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/
funannotate annotate --gff /data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/genemark.longest.gff3 --fasta /data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa --iprscan tmp2.faa -o annotFunTest

threads=30
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
interpro=/data/blaxter/Otipulae_lr/raw/funannotate_db/interproscan-5.44-79.0/interproscan.sh
prots=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/genemark.longest.faa
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/functionalAnnot
assemName=interpro5_oscheius_tipulae.local.v3_3


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ceph/users/pgonzales/miniconda3/envs/kat/lib/; $interpro -i $prots -cpu $threads -f xml"

# Kinfin interproscan
cd /data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/
cat /data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/*faa > /data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/all_longestProteins.faa

threads=30
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
interpro=/data/blaxter/Otipulae_lr/raw/funannotate_db/interproscan-5.44-79.0/interproscan.sh
prots=/data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/all_longestProteins.faa
outdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/
assemName=interpro5_longestIsoforms


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp64 $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ceph/users/pgonzales/miniconda3/envs/kat/lib/; $interpro -i $prots -cpu $threads -dp -t p --goterms -appl SignalP-EUK-4.1 -appl Pfam -f TSV"

sconda /ceph/software/conda/envs/kinfin_env/
/ceph/software/kinfin/kinfin/scripts/iprs2table.py -i all_longestProteins.faa.tsv
orthoDir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/WorkingDirectory/OrthoFinder/Results_May27
sed 's/.faa//' ${orthoDir}/WorkingDirectory/SpeciesIDs.txt | awk -F ": " 'BEGIN{print"#IDX,TAXON,OUT"; OFS = ","}{$3=0; if($2 ~ /ascaris/){$3 = 1}; print $0}' > config.txt
/ceph/software/kinfin/kinfin/kinfin -c config.txt -g ${orthoDir}/Orthogroups/Orthogroups.txt -t ${orthoDir}/Species_Tree/SpeciesTree_rooted.txt -s ${orthoDir}/WorkingDirectory/SequenceIDs.txt -p ${orthoDir}/WorkingDirectory/SpeciesIDs.txt --functional_annotation functional_annotation.txt --target_fraction 0.84 -o /data/blaxter/Otipulae_lr/analyses/orthofinder/kinfin/functionalAnnot

#emapper
threads=30
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
prots=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/genemark.longest.faa
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/functionalAnnot
assemName=eggnog_oscheius_tipulae.local.v3_3


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; emapper.py -i $prots --output $outdir/${assemName}.pfam -m diamond --usemem --cpu $threads --data_dir /data/blaxter/Otipulae_lr/raw/funannotate_db/"

# Funannotate annotate

threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
gff=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/genemark.longest.gff3
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/functionalAnnot
iprscan=$outdir/genemark.longest.faa.xml
eggnog=$outdir/eggnog_oscheius_tipulae.local.v3_3.pfam.emapper.annotations
assemName=funannotate_oscheius_tipulae.local.v3_3


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ceph/users/pgonzales/.conda/envs/funannotate/lib/; funannotate annotate --gff $gff --fasta $assembly --iprscan $iprscan --eggnog $eggnog -s Ot3_3 -o $outdir/$assemName --cpus $threads"




source /ceph/software/conda/bin/activate /ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/splitProtsByTelomere
prots=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/genemark.longest.faa
genePos=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/oscheius_tipulae.local.v3_3.protCodGenes.tsv.gz


cd $outdir
for endCoord in $(ls $outdir/rangesDir/*bed); do
teloEnd=$(basename $endCoord .bed)
zcat $genePos | \
 bedtools intersect -a stdin -b $endCoord | cut -f 4 > ${teloEnd}.txt
# Parse ID
cat ${teloEnd}.txt | xargs samtools faidx $prots > protSeqs/${teloEnd}.faa
done
cd protSeqs/
tar -czvf otipulaev3.3.telomericExtension.ProteinSequences.wTEs.tar.gz *faa

# add names and intepro results

threads=1
condaEnv=/ceph/users/pgonzales/.conda/envs/gffmagic
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
gff=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/longest_isoform/genemark.longest.gff3
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/functionalAnnot
iprscan=$outdir/genemark.longest.faa.xml
eggnog=$outdir/eggnog_oscheius_tipulae.local.v3_3.pfam.emapper.annotations
assemName=funannotate_oscheius_tipulae.local.v3_3

cd $outdir
awk -F '\t' 'BEGIN{OFS=FS}{print $1, $2, "100", "100", "0", "0", "0", "100", "0", "100", "0", "1000"}' funannotate_oscheius_tipulae.local.v3_3/annotate_misc/uniprot_eggnog_raw_names.txt > pseudoBlast.tsv
awk -F '\t' 'BEGIN{OFS=FS}{print ">"$2"|"$3"\n"}' funannotate_oscheius_tipulae.local.v3_3/annotate_misc/uniprot_eggnog_raw_names.txt > pseudoSubject.fasta


agat_sp_manage_functional_annotation.pl -f $gff -i genemark.longest.faa.tsv --output tst.gff3
-b pseudoBlast.tsv --db pseudoSubject.fasta --output tst.gff3




# Check chromosome concordance
# Plot Ce chromosome correspondence
# Plot telomeric repeats




threads=12
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/annotation/
prots=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/protSeqs_agat/genemark.longest.faa
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=funannotate_gm_$(basename $assembly .fa)
transcripts=/data/blaxter/Otipulae_lr/analyses/RNA-assembly/mira/otipulae_ERR1880912.fasta
brakerGff=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/agat_gff3/braker_abinitio.gff3








threads=1
condaEnv=/ceph/users/pgonzales/.conda/envs/EMBLmyGFF3
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/repeats/otipu
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_1rx.genomic.fa
assemName=oscheius_tipulae.local.v3_1rx


qsub -cwd -o ${outdir}/einvert_${assemName}.log -e ${outdir}/einvert_${assemName}.err -N einvert_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; einverted -sequence $assembly -outfile $assemName"

rmask_oscheius_tipulae.local.v3_3.sh








/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/ltr_finder/ltrFinder_oscheius_tipulae.local.v3_3.sh
cat ${assemName}.harvest.scn $(basename $assembly).finder.combine.scn > $(basename $assembly .genomic.fa).rawLTR.scn
/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/ltr_finder/ltrRetrieve_oscheius_tipulae.local.v3_3.sh


threads=60
condaEnv=/ceph/users/pgonzales/miniconda3/envs/LTR_retriever
ltrFinder=/data/blaxter/Otipulae_lr/scripts/LTR_FINDER_parallel/LTR_FINDER_parallel
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/ltr_finder
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=ltrFinder_$(basename $assembly .genomic.fa)

cmd="$ltrFinder -seq $assembly -threads $threads -harvest_out -size 1000000 > $outdir/${assemName}.out"

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"




# Plot features



threads=6
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/funannotate
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=funannotate_ceProts_$(basename $assembly .fa)
transcripts=/data/blaxter/Otipulae_lr/analyses/RNA-assembly/mira/otipulae_ERR1880912.fasta
proteins=/data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/caenorhabditis_elegans.PRJNA13758.WBPS14.faa
cmd="export FUNANNOTATE_DB=/data/blaxter/Otipulae_lr/raw/funannotate_db/; funannotate predict --GENEMARK_PATH /data/blaxter/Otipulae_lr/scripts/gmes_linux_64/ --cpus $threads --force -i $assembly -o $outdir/$assemName -s \"Oscheius tipulae\" --name OT3 --augustus_species caenorhabditis --organism other --transcript_evidence $transcripts --protein_evidence $proteins"



qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; $cmd"




threads=40
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/repeats/repeatmodeler/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=$(basename $assembly .fa)



qsub -cwd -o ${outdir}/RM_${assemName}.log -e ${outdir}/RM_${assemName}.err -N RM_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate /ceph/software/conda/envs/repeat_env; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; ln -s $assembly .;  /ceph/software/repeatmodeler/RepeatModeler.v2.0.1/BuildDatabase -name $assemName -engine ncbi ${assemName}.fa; /ceph/software/repeatmodeler/RepeatModeler.v2.0.1/RepeatModeler -pa $threads -database $assemName; rm -f ${assemName}.fa" 






threads=60
RAMmem=60
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/assemblies/freebayes/
pe_reads1=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_1.fastq.gz
pe_reads2=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_2.fastq.gz
assembly=/data/blaxter/Otipulae_lr/analyses/assemblies/gap5/O.tipAlanCuratedPreCorrection.fa
assemName=freebayes_otipulaeA
wkdir=/scratch/pgonzales/snippy/$assemName


qsub -cwd -o ${outdir}/snippy_${assemName}.log -e ${outdir}/snippy_${assemName}.err -N snippy_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; perl -I /ceph/software/perl-ensembl/BioPerl-1.6.924/ /ceph/users/pgonzales/miniconda3/envs/callvars/bin/snippy --cpus $threads --ram $RAMmem --outdir $wkdir --ref $assembly --R1 $pe_reads1 --R2 $pe_reads2 && cp $wkdir/snps.consensus.fa $outdir/${assemName}.fa && cp $wkdir/snps.log $outdir/${assemName}.snps.log"


threads=30
RAMmem=90
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/assemblies/freebayes/
pe_reads1=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_1.fastq.gz
pe_reads2=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_2.fastq.gz
assembly=/data/blaxter/Otipulae_lr/analyses/assemblies/freebayes/freebayes3_racon_otipulaeA.fa
assemName=freebayes4_racon_otipulaeA
wkdir=/scratch/pgonzales/snippy/$assemName


qsub -cwd -o ${outdir}/snippy_${assemName}.log -e ${outdir}/snippy_${assemName}.err -N snippy_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; perl -I /ceph/software/perl-ensembl/BioPerl-1.6.924/ /ceph/users/pgonzales/miniconda3/envs/callvars/bin/snippy --cpus $threads --ram $RAMmem --outdir $wkdir --ref $assembly --R1 $pe_reads1 --R2 $pe_reads2 && cp $wkdir/snps.consensus.fa $outdir/${assemName}.fa && cp $wkdir/snps.log $outdir/${assemName}.snps.log"







assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
illumiBam=/scratch/pgonzales/snippy/freebayes4_otipulaeA/snps.bam
promeBam=/scratch/pgonzales/minimap2/oti3.2/oti3.2.bam
intrstRegions=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/igv-reports/freebayes3_otipulaeA.bed
assemName=freebayes3_otipulaeA

qsub -cwd -o ${outdir}/snippy_${assemName}.log -e ${outdir}/snippy_${assemName}.err -N snippy_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; perl -I /ceph/software/perl-ensembl/BioPerl-1.6.924/ /ceph/users/pgonzales/miniconda3/envs/callvars/bin/snippy --cpus $threads --ram $RAMmem --outdir $wkdir --ref $assembly --R1 $pe_reads1 --R2 $pe_reads2 && cp $wkdir/snps.consensus.fa $outdir/${assemName}.fa && cp $wkdir/snps.log $outdir/${assemName}.snps.log"

awk 'BEGIN{print"#IDX,TAXON"; FS = ": "; OFS = ","}{print $1,$2}' /data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/WorkingDirectory/SpeciesIDs.txt

awk 'BEGIN{FS = ": "; OFS = ","}{print $0}' /data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/WorkingDirectory/SpeciesIDs.txt



orthofinder -t $threads -a $threads -M msa -S blast -b $wkdir -f $protsDir

protsDir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/ascaris_suum.local.v2020/
outDir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3
wkdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/WorkingDirectory/
threads=64


qsub -cwd -o ${outDir}/orthofinder_wAsuum.log -e ${outDir}/orthofinder_wAsuum.err -N orthofinder_wAsuum -V -pe smp64 $threads -b yes "source /ceph/software/conda/bin/activate /ceph/users/pgonzales/.conda/envs/orthofinder; mkdir -p $wkdir; echo localhost: \$HOSTNAME; echo "orthofinder version:" \$(orthofinder --version | head -n 2 | tail -n +2); echo outDir: $outDir; cd $outDir; echo cmnd: ${outDir}/exec_orthofinder.sh; cat ${outDir}/exec_orthofinder.sh; echo date \$(date); orthofinder -t $threads -a $threads -M msa -S blast -b $wkdir -f $protsDir"



threads=1
condaEnv=/ceph/software/conda/envs/kinfin_env/
outdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/kinfin
orthoFdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/
spIds=$orthoFdir/WorkingDirectory/SpeciesIDs.txt
orthos=$orthoFdir/Orthogroups/Orthogroups.txt
targetFrac=0.84
assemName=kinfin_$targetFrac

for 

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; sed 's/.faa//' $orthoFdir/WorkingDirectory/SpeciesIDs.txt | awk -F \": \" 'BEGIN{print\"#IDX,TAXON\"; OFS = \",\"}{print \$1,\$2}' > config.txt; /ceph/software/kinfin/kinfin/kinfin -c config.txt -g $orthoFdir/Orthogroups/Orthogroups.txt -s $orthoFdir/WorkingDirectory/SequenceIDs.txt -p $orthoFdir/WorkingDirectory/SpeciesIDs.txt --target_fraction $targetFrac -o $outdir/kinfin_$targetFrac"




threads=1
condaEnv=/ceph/software/conda/envs/kinfin_env/
outdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/kinfin
orthoFdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/
spIds=$orthoFdir/WorkingDirectory/SpeciesIDs.txt
orthos=$orthoFdir/Orthogroups/Orthogroups.txt
targetFrac=0.84
assemName=kinfin_$targetFrac

qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; sed 's/.faa//' $orthoFdir/WorkingDirectory/SpeciesIDs.txt | awk -F \": \" 'BEGIN{print\"#IDX,TAXON\"; OFS = \",\"}{print \$1,\$2}' > config.txt; /ceph/software/kinfin/kinfin/kinfin -c config.txt -g $orthoFdir/Orthogroups/Orthogroups.txt -s $orthoFdir/WorkingDirectory/SequenceIDs.txt -p $orthoFdir/WorkingDirectory/SpeciesIDs.txt --target_fraction $targetFrac -o $outdir/kinfin_$targetFrac"


threads=1
condaEnv=/ceph/software/conda/envs/kinfin_env/
outdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/kinfin
orthoFdir=/data/blaxter/Otipulae_lr/analyses/orthofinder/otipu3.3/OrthoFinder/Results_May26_1/
spIds=$orthoFdir/WorkingDirectory/SpeciesIDs.txt
orthos=$orthoFdir/Orthogroups/Orthogroups.txt
targetFrac=0.84
assemName=kinfin_$targetFrac


orthoFile=$orthoFdir/



threads=1
host=c6
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/igv-reports/
assembly=/data/blaxter/Otipulae_lr/analyses/assemblies/freebayes/freebayes3_otipulaeA.fa
illumiBam=/scratch/pgonzales/snippy/freebayes4_otipulaeA/snps.bam
promeBam=/scratch/pgonzales/minimap2/oti3.2/oti3.2.bam
intrstRegions=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/igv-reports/freebayes3_otipulaeA.bed
assemName=freebayes3_otipulaeA



qsub -cwd -o ${outdir}/igvRep_${assemName}.log -e ${outdir}/igvRep_${assemName}.err -N igvRep_${assemName} -V -pe smp $threads -l h=$host -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; create_report $intrstRegions $assembly --tracks $illumiBam $promeBam --output ${outdir}/${assemName}.html"


outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/comparedAnnots/rbbh/
refProt=/data/blaxter/Otipulae_lr/analyses/orthofinder/longestIsoforms/caenorhabditis_elegans.PRJNA13758.WBPS14.faa
gprot=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.3/proteinCoding/benchmark/funannotate/funannotate_gm_oscheius_tipulae.local.v3_3.genomic/predict_results/Oscheius_tipulae.proteins.fa
assemName=funannotate_gm_braker_Ce
nextflow -C /data/blaxter/Otipulae_lr/scripts/rbbh.config run /data/blaxter/Otipulae_lr/scripts/rbbh.nf --query $gprot --ref $refProt --out $outdir/${assemName}_vs_celeg.txt --chunkSize 500 --minqueryCov 20
rm -rf work 



outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/minimap2
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
readset=/data/blaxter/Otipulae_lr/analyses/decontamination/kraken2/dbs/kraken_clean_reads.fasta
assemName=oti3.3
threads=64
wkDir=/scratch/pgonzales/minimap2/$assemName
qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N minimap_$assemName -V -pe smp64 $threads -b yes "source /ceph/software/conda/bin/activate /ceph/users/pgonzales/miniconda3/envs/redbean; echo localhost: \$HOSTNAME; echo "busco version:" \$(busco --version); echo outDir: $outdir; cd $outdir; echo conda_env: \$CONDA_PREFIX; echo \"### packages\"; conda list --export; echo \"###\"; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; echo date \$(date); mkdir -p $wkDir; cd $wkDir; minimap2 --MD -L -ax map-ont -t $threads $assembly $readset | samtools sort -@ 3 - | samtools view -b > ${assemName}.bam; samtools index -@ $threads ${assemName}.bam"


threads=1
host=bigbang
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/minimap2/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=depth_oscheius_tipulae.local.v3_3
bam=/scratch/pgonzales/minimap2/oti3.3/oti3.3.bam


qsub -hold_jid 47609 -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -l h=$host -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; samtools depth -Q 10 -d 0 $bam > $outdir/${assemName}_promethion.cov"


threads=1
host=bigwig
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/bwa/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=freebayes_oscheius_tipulae.local.v3_3
bam=/scratch/pgonzales/snippy/freebayes_oscheius_tipulae.local.v3_3/snps.bam


qsub -hold_jid 47610 -cwd -o ${outdir}/depth_${assemName}.log -e ${outdir}/depth_${assemName}.err -N depth_${assemName} -V -pe smp64 $threads -l h=$host -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; samtools depth -d 0 $bam > $outdir/${assemName}_ERR1656469.cov"



threads=64
RAMmem=150
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/assemblies/freebayes/
pe_reads1=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_1.fastq.gz
pe_reads2=/data/blaxter/Otipulae_lr/analyses/assemblies/trimmd_Illumina/ERR1656469_2.fastq.gz
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=freebayes_oscheius_tipulae.local.v3_3
wkdir=/scratch/pgonzales/snippy/$assemName


qsub -cwd -o ${outdir}/snippy_${assemName}.log -e ${outdir}/snippy_${assemName}.err -N snippy_${assemName} -V -pe smp64 $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; perl -I /ceph/software/perl-ensembl/BioPerl-1.6.924/ /ceph/users/pgonzales/miniconda3/envs/callvars/bin/snippy --cpus $threads --ram $RAMmem --outdir $wkdir --ref $assembly --R1 $pe_reads1 --R2 $pe_reads2 && cp $wkdir/snps.consensus.fa $outdir/${assemName}.fa && cp $wkdir/snps.log $outdir/${assemName}.snps.log"


threads=1
host=bigwig
condaEnv=/ceph/users/pgonzales/miniconda3/envs/callvars
outdir=/data/blaxter/Otipulae_lr/analyses/qualityMetrics/mappingStats/bwa/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_3.genomic.fa
assemName=freebayes_oscheius_tipulae.local.v3_3
bam=/scratch/pgonzales/snippy/freebayes_oscheius_tipulae.local.v3_3/snps.bam


qsub -hold_jid 47610 -cwd -o ${outdir}/depth_${assemName}.log -e ${outdir}/depth_${assemName}.err -N depth_${assemName} -V -pe smp64 $threads -l h=$host -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; samtools depth -d 0 $bam > $outdir/${assemName}_ERR1656469.cov"





threads=30
condaEnv=/ceph/users/pgonzales/.conda/envs/funannotate
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/otipu3.2/ncRNA
assembly=/data/blaxter/Otipulae_lr/analyses/assemblies/freebayes/freebayes3_racon_otipulaeA.fa
assemName=infernal_fb3_otipulaeA
wkdir=/scratch/pgonzales/infernal/fb3_V/


qsub -cwd -o ${outdir}/${assemName}.log -e ${outdir}/${assemName}.err -N ${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; mkdir -p $wkdir; samtools faidx $assembly V > $wkdir/${assemName}_V.fa; cmscan --cut_nc --rfam --cpu $threads /data/blaxter/Otipulae_lr/raw/Rfam.cm $wkdir/${assemName}_V.fa"




stripchart(decrease ~ treatment,
           main = "stripchart(OrchardSprays)",
           vertical = TRUE, log = "y", data = OrchardSprays, method = "stack", pch = "_")



###
TransposonPSI
###

threads=1
condaEnv=/ceph/users/pgonzales/.conda/envs/legacyBlast
outdir=/data/blaxter/Otipulae_lr/analyses/annotation/repeats/otipu/tranposonPSI
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_1rx.genomic.fa
assemName=oscheius_tipulae.local.v3_1rx


qsub -cwd -o ${outdir}/psi_${assemName}.log -e ${outdir}/psi_${assemName}.err -N psi_${assemName} -V -pe smp $threads -b yes "cd $outdir; source /ceph/software/conda/bin/activate $condaEnv; echo localhost: \$HOSTNAME; echo outDir: $outdir; echo conda_env: \$CONDA_PREFIX; echo cmnd: ${outdir}/${assemName}.sh; cat ${outdir}/${assemName}.sh; perl /data/blaxter/Otipulae_lr/scripts/TransposonPSI_08222010/transposonPSI.pl $assembly nuc"




sconda medaka
cd /data/blaxter/Otipulae_lr/analyses/repeatmodeler/

bash /data/blaxter/Otipulae_lr/analyses/repeatmodeler/oscheius_tipulae.local.v3_1rx.sh

threads=40
outdir=/data/blaxter/Otipulae_lr/analyses/repeatmodeler/
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_1rx.genomic.fa
assemName=oscheius_tipulae.local.v3_1rx

cat oscheius_tipulae.local.v3_1rx-families.fa <(/ceph/users/pgonzales/.conda/envs/ncRNA/bin/perl /ceph/software/repeatmasker/RepeatMasker-4.1.0/util/queryRepeatDatabase.pl -species nematoda) > allNema_otipu_v3_1rx.fa



threads=32
rmaskExeDir=/ceph/software/repeatmasker/RepeatMasker-4.1.0
outdir=/data/blaxter/Otipulae_lr/analyses/repeatmasker/
repLib=/data/blaxter/Otipulae_lr/analyses/repeatmodeler/allNema_otipu_v3_1rx.fa
assembly=/data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_1rx.genomic.fa
assemName=oscheius_tipulae.local.v3_1rx


bash ../repeatmasker/rmask_otipu_v3_1rx.sh

../repeatmasker/oscheius_tipulae.local.v3_1rx.genomic.fa.out.gff


outDir=/data/blaxter/Otipulae_lr/analyses/comparGenomics/repRegions
windSize=50000

tail -n +4 ../repeatmasker/oscheius_tipulae.local.v3_1rx.genomic.fa.out.gff | sort -k1,1 -k4,4n > sorted_otipu.gff
bedtools merge -i sorted_otipu.gff > merged_otipu.bed

awk -F '\t' '!/^ *#/{sum1+=$2; sum2+=$3} END{print sum2 - sum1}' merged_otipu.bed
awk -F '\t' '{sum1+=$2} END{print sum1}' /data/blaxter/Otipulae_lr/raw/alternativeGenomeSources/assembly/oscheius_tipulae.local.v3_1rx.genomic.fa.fai

for tsvF in $(ls spDict/*); do 
outFile=$(echo bedNemaPos/$(basename $tsvF) | sed 's/tsv/bed/')
awk -F '\t' '{print $1, $2, $3, $5}' $tsvF > $outFile
done


bedtools makewindows -g /data/blaxter/Otipulae_lr/raw/genome/caenorhabditis_elegans.PRJNA13758.WBPS12.genomic_softmasked.fa.fai -w $windSize > genome.bed
bedtools coverage -a genome.bed -b merged_celeg.bed > $outDir/caenorhabditis_elegans.PRJNA13758.WBPS12.genomic_softmasked_softDens.bed