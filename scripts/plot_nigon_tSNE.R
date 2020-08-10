library(Rtsne)
library(tidyverse)
library(reshape2)
library(cluster)

# Perplexity values ro explore
tsPerplex = c(10,30,50,100,200,500,1000)

# Load previous Nigon definition
nigonDict <- read_tsv("analyses/orthoFinder/gene2Nigon20200528_09_55.tsv.gz")
# Load all gene positions
gPos <- read_tsv("analyses/annotation/transcript_positions.tsv.gz")
# Load identifiers of 1-to-1 orthogroups
kinFile <- paste0("analyses/kinfin/all/all.all.cluster_1to1s.txt")
# Load Orthofinder result
orthos <- read_tsv("analyses/orthoFinder/Results_May27/Orthogroups/Orthogroups.tsv") %>%
  melt(id.vars = "Orthogroup") %>%
  as_tibble() %>%
  # Edit gene IDs from "_" to ":". 
  # Assuming orthofinder formated these IDS
  mutate(value = ifelse(!variable %in% c("auanema_rhodensis.local.CHRR_integrated",
                                         "oscheius_tipulae.local.v3_3",
                                         "ascaris_suum.local.v2020"),
                        sub("_", ":", value), value))

# List of SCO orthorgroups
one2ones <- read_tsv(kinFile)
ogroups <- filter(orthos, Orthogroup %in% one2ones$`#cluster_id`)

# True single copy orthologs
allSps <- c("auanema_rhodensis.local.CHRR_integrated", "brugia_malayi.PRJNA10729.WBPS14", "caenorhabditis_elegans.PRJNA13758.WBPS14", "haemonchus_contortus.PRJEB506.WBPS14", "onchocerca_volvulus.PRJEB513.WBPS14", "oscheius_tipulae.local.v3_3", "pristionchus_pacificus.PRJNA12644.WBPS14", "steinernema_carpocapsae.GCA_000757645.3_ASM75764v3", "strongyloides_ratti.PRJEB125.WBPS14")
one2ones <- filter(orthos, variable %in% allSps) %>%
  filter(!is.na(value)) %>%
  filter(!grepl(",", value)) %>%
  group_by(Orthogroup) %>%
  summarise(tDif = n_distinct(variable)) %>%
  ungroup() %>% 
  filter(tDif == length(allSps)) %>%
  pull(Orthogroup)
ogroups <- filter(orthos, Orthogroup %in% one2ones)


# split multigroup 
# I will discard those found in multiple chromosomes
# because I cannot say anything about them
# regarding their movement between chromosomes
# they are neither evidence for nor against
# macrosynteny

uniOgroup <- filter(ogroups, !grepl(",", value))
multiOgroup <- filter(ogroups, grepl(",", value))

lmultiOgroup <- map2(multiOgroup$value, 1:nrow(multiOgroup), function(x, y){
  spt_gs <- as.character(unlist(str_split(x, ", ")))
  cbind(multiOgroup[rep(y, length(spt_gs)), -3],
        value = spt_gs,
        stringsAsFactors = F)
}) %>%
  bind_rows() %>%
  as_tibble()

allogroups <- bind_rows(uniOgroup, lmultiOgroup) %>%
  # Assign each gene to its corresponding chromosome
  left_join(gPos, by=c("variable" = "assembly", "value"="tID")) %>%
  select(-stPos, -ePos, -strnd, -qsp)

# Filter out species orthogroups that
# have more than one chromosome
multiChroms <- group_by(allogroups, Orthogroup, variable) %>%
  summarise(dChroms = n_distinct(scaffold)) %>%
  ungroup() %>% # count(dChroms)
  filter(dChroms > 1)

uniChroms <- anti_join(allogroups, multiChroms,
                       by = c("Orthogroup", "variable")) %>%
  group_by(Orthogroup) %>%
  # Leave only one gene per species as 
  # chromosome representative
  filter(!duplicated(variable)) %>%
  ungroup()




# Convert to wide format
ogrp2Chroms <- group_by(uniChroms, Orthogroup, variable) %>%
  summarise(scaffold = paste(unique(scaffold), collapse = ",")) %>%
  ungroup %>%
  dcast(Orthogroup ~ variable, value.var = "scaffold", fill = NA) %>%
  as.matrix()

ogrp2Chroms[ogrp2Chroms == "NA"] <- NA

# Remove C. briggsae, C. nigoni, C. remanei and C. inopinata
nonRednd <- as_tibble(ogrp2Chroms) %>%
  select(Orthogroup, auanema_rhodensis.local.CHRR_integrated,
         brugia_malayi.PRJNA10729.WBPS14,
         caenorhabditis_elegans.PRJNA13758.WBPS14,
         haemonchus_contortus.PRJEB506.WBPS14,
         onchocerca_volvulus.PRJEB513.WBPS14,
         oscheius_tipulae.local.v3_3,
         pristionchus_pacificus.PRJNA12644.WBPS14,
         steinernema_carpocapsae.GCA_000757645.3_ASM75764v3,
         strongyloides_ratti.PRJEB125.WBPS14) %>%
  column_to_rownames("Orthogroup")
nonRednd <- apply(nonRednd, 2, factor) %>%
  as.data.frame()


gower.dist <- 1-as.matrix(daisy(nonRednd, metric = c("gower")))
colnames(gower.dist) <- rownames(nonRednd)
rownames(gower.dist) <- rownames(nonRednd)

gower.dist[is.na(gower.dist)] <- 0



uniNigDict <- filter(nigonDict, !duplicated(Orthogroup))

# Colors for Nigons
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "-" = "#000000")

for (perplx in tsPerplex) {
  tsne <- Rtsne(gower.dist, dims = 2,
                perplexity=perplx,
                verbose=TRUE, max_iter = 1000,
                check_duplicates = FALSE)
  etsne <-data.frame(cbind(tsne$Y), Orthogroup = rownames(gower.dist))

  # Plot t-SNE
  mutate(etsne, nigon = uniNigDict$nigon[match(etsne$Orthogroup, uniNigDict$Orthogroup)],
         Nigon = ifelse(is.na(nigon), "-", nigon),
         Nigon = factor(Nigon, levels = c(LETTERS[1:5], "N", "X", "-"))) %>%
    arrange(desc(Nigon)) %>%
    ggplot(aes(colour=Nigon, y=X1, x=X2)) + 
    geom_point(size = 0.8) +
    scale_colour_manual(values = cols) +
    theme_minimal() +
    ggtitle(paste0("tSNE, perplexity = ", perplx))
  
  ggsave(paste0("report/figures/extras/gower_tSNE_",
                perplx, ".png"), width = 5, height = 4)
  
}