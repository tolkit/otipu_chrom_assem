# Nigonize orthologs

library(reshape2)
library(cluster)
library(tidyverse)


# Prior from Tandonnet et al. (2017) model.
# Used ONLY for giving Nigon unit names
priorNigon <- tibble(chromosome = c("ce_I", "ce_II", "ce_III", "ce_IV", "ce_V", "ot_V", "pp_ChrX"),
                     nigon = c(LETTERS[1:5], "N", "X"))


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

# Cluster Gower distances using CLARA
fit <- clara(gower.dist, 7, sampsize = 0.1 * nrow(nonRednd))

ortoWclust <- cbind(nonRednd,
                    cluster = fit$cluster[match(names(fit$cluster), rownames(gower.dist))],
                    Orthogroup = rownames(gower.dist))


chromClust_dict <- melt(ortoWclust, id.vars = c("Orthogroup", "cluster"), value.name = "chromosome") %>%
  group_by(variable, chromosome) %>%
  count(cluster) %>%
  filter(n > max(n)*0.2) %>%
  ungroup %>%
  select(-n) %>%
  reshape2::dcast(chromosome ~ variable,
                  value.var = "cluster",
                  fun.aggregate=function(x) paste0(unique(x), collapse = ",")) %>%
  melt(id.vars = "chromosome", value.name = "cluster") %>%
  filter(cluster != "", !is.na(chromosome)) %>%
  as_tibble()


clustogroups <- inner_join(uniChroms, chromClust_dict,
                           by = c("variable", "scaffold" = "chromosome"))

ortoAndClustSupp <- group_by(clustogroups, Orthogroup) %>%
  summarise(cl1 = sum(str_detect(cluster, "1")),
            cl2 = sum(str_detect(cluster, "2")),
            cl3 = sum(str_detect(cluster, "3")),
            cl4 = sum(str_detect(cluster, "4")),
            cl5 = sum(str_detect(cluster, "5")),
            cl6 = sum(str_detect(cluster, "6")),
            cl7 = sum(str_detect(cluster, "7"))) %>%
  ungroup() %>%
  melt(id.vars = "Orthogroup")

cntClust <-  group_by(ortoAndClustSupp, Orthogroup) %>%
  top_n(1, value) %>%
  ungroup() %>%
  rename(oGroupClust = variable,
         oGroupClustN = value) %>%
  mutate(oGroupClust = sub("cl", "", oGroupClust))


# Remove orthorgoups that had a tie in the most
# supported cluster.
cntClust <- filter(cntClust, !Orthogroup %in% cntClust$Orthogroup[duplicated(cntClust$Orthogroup)])


# How many are supported
inner_join(clustogroups, cntClust,
           by = "Orthogroup") %>%
  group_by(Orthogroup) %>%
  arrange(desc(oGroupClustN)) %>%
  slice(1) %>%
  ungroup() %>%
  count(oGroupClustN)


# Create dict with only highly consistent orthorgoups
consOrts <- inner_join(clustogroups, cntClust,
                       by = "Orthogroup") %>%
  # More than 7 assemblies supporting the majority cluster  
  filter(oGroupClustN > 7) %>%
  pull(Orthogroup) %>% unique


# Link scaffolds to most abundant cluster
consChromClust <- filter(ortoWclust, Orthogroup %in% consOrts) %>%
  melt(id.vars = c("Orthogroup", "cluster"), value.name = "chromosome") %>%
  group_by(variable, chromosome) %>%
  count(cluster) %>%
  filter(n > max(n)*0.2) %>%
  ungroup %>%
  select(-n) %>%
  reshape2::dcast(chromosome ~ variable,
                  value.var = "cluster",
                  fun.aggregate=function(x) paste0(unique(x), collapse = ",")) %>%
  melt(id.vars = "chromosome", value.name = "cluster") %>%
  filter(cluster != "", !is.na(chromosome)) %>%
  as_tibble()

# Translate Nigons to clusters
nigon2clust <- mutate(priorNigon,
                      cluster = consChromClust$cluster[match(chromosome, consChromClust$chromosome)])

# Translate clusters to Nigons
gene2Nigon <- filter(ortoWclust, Orthogroup %in% consOrts) %>%
  mutate(nigon = nigon2clust$nigon[match(cluster, nigon2clust$cluster)]) %>%
  select(Orthogroup, nigon) %>%
  right_join(uniChroms, by = "Orthogroup") %>%
  as_tibble() %>%
  filter(!is.na(nigon)) %>%
  inner_join(gPos, by = c("value" = "tID",
                          "variable" = "assembly",
                          "scaffold")) %>%
  select(Orthogroup, nigon, assembly = variable,
         tID = value, scaffold, stPos)


write_tsv(gene2Nigon, paste0("nalyses/orthoFinder/gene2Nigon",
                             format(Sys.time(), "%Y%m%d_%H_%M"),
                             ".tsv.gz"))

