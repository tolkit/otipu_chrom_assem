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

# Start an R terminal
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
# Close terminal
exit()

sort -k1,1 -k2,2n ${assemName}.ceChr.tsv | bedtools merge -i stdin -c 4 -o distinct > ${assemName}.ceChr.merged.tsv
bedtools makewindows -g ${assembly}.fai -w 1000 | bedtools intersect -wb -a stdin -b ${assemName}.ceChr.merged.tsv | awk -F '\t' 'BEGIN{OFS=FS}{print $1, $2, $3, $7}' > ceChr_${assemName}.1kw.tsv


# Nigon density

# Start an R terminal
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
# Close terminal
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
