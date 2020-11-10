library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(scales)
library(gtools)
library(ggplot2)
library(cowplot)


nigonDict <- read_tsv("analyses/orthoFinder/gene2Nigon20200528_09_55.tsv.gz") %>%
  select(-kinfCut) %>%
  mutate(assembly = ifelse(assembly == "caenorhabditis_sp34.PRJDB5687.WBPS14",
                           "caenorhabditis_inopinata.PRJDB5687.WBPS14",
                           assembly))

### Add Bursaphelenchus
bursaphelenchus_nigonDict <- read_tsv("analyses/orthoFinder/bursaphelenchus_gene2Nigon20201109_14_00.tsv.gz")

### Add Heterorhabditis bacteriophora
hbacteriophoraDict <- read_tsv("analyses/orthoFinder/hbacteriophora_gene2Nigon20201109_14_00.tsv.gz")

### Add diploscapter coronatus
diplosdictNX <- read_tsv("analyses/orthoFinder/dcoronatusNX_gene2Nigon20201109_14_00.tsv.gz")


### Concatenate
nigonDict_allSp <- bind_rows(nigonDict, hbacteriophoraDict,
                             diplosdictNX, bursaphelenchus_nigonDict)


write_tsv(nigonDict_allSp, "analyses/orthoFinder/gene2Nigon20201109_14_00.tsv.gz")


cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Plot function
colrNigon <- function(species, windwSize){
  spName <- sub("\\..+", "", species)
  spName <- sub("_", " ", spName)
  spName <- paste(toupper(substr(spName, 1, 1)), substr(spName, 2, nchar(spName)), sep="")
  
  filter(nigonDict_allSp, assembly == species) %>%
    mutate(scaffold = sub("[cC]hr", "", sub("^.._", "", scaffold)),
           scaffold = sub("TMUE_", "", scaffold),
           Nigon = nigon) %>%
    group_by(scaffold) %>%
    mutate(nGenes = n(),
           mxGpos = max(stPos)) %>%
    ungroup() %>%
    filter(nGenes > 15, mxGpos > windwSize * 2) %>%
    group_by(scaffold) %>%
    mutate(ints = as.numeric(as.character(cut(stPos,
                                              breaks = seq(0, max(stPos), windwSize),
                                              labels = seq(windwSize, max(stPos), windwSize)))),
           ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
    count(ints, Nigon) %>%
    ungroup() %>%
    mutate(scaffold_f = factor(scaffold,
                               levels = mixedsort(unique(scaffold)))) %>%
    ggplot(aes(fill=Nigon, y=n, x=ints-windwSize)) + 
    facet_grid(scaffold_f ~ ., switch = "y") +
    geom_bar(position="stack", stat="identity") +
    theme_minimal() +
    scale_y_continuous(position = "right") +
    scale_x_continuous(labels = label_number_si()) +
    scale_fill_manual(values = cols) +
    ggtitle(paste0("*", spName, "*")) +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          panel.border = element_blank(),
          plot.title = ggtext::element_markdown(),
          legend.position = "none")
}


# plot 9 species per page when they have less than 8 scaffolds

## first page

spVec <- grep("enorhabditis|conto|rhod|paci|stron", unique(nigonDict_allSp$assembly),
     value = T)

myplots <- list()
for(species in spVec){
  myplots[[species]] <- colrNigon(species, 5e5)
}

p1 <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]],
                myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]],
                myplots[[9]],
          labels = "AUTO", label_size = 12, nrow = 3)

ggsave("p1_t.pdf", p1, width = 21, height = 29.7, units = "cm")

## second page, with legend
exampPlot <- filter(nigonDict_allSp, assembly == "caenorhabditis_elegans.PRJNA13758.WBPS14") %>%
  mutate(scaffold = sub("[cC]hr", "", sub("^.._", "", scaffold)),
         scaffold = sub("TMUE_", "", scaffold),
         Nigon = nigon) %>%
  group_by(scaffold) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > 15, mxGpos > 5e5 * 2) %>%
  group_by(scaffold) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), 5e5),
                                            labels = seq(5e5, max(stPos), 5e5)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + 5e5, ints)) %>%
  count(ints, Nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(scaffold,
                             levels = mixedsort(unique(scaffold)))) %>%
  ggplot(aes(fill=Nigon, y=n, x=ints-5e5)) + 
  facet_grid(scaffold_f ~ ., switch = "y") +
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = cols)

sxP <- get_legend(exampPlot)

spVec <- grep("steine|bursa|brug|oncho", unique(nigonDict_allSp$assembly),
              value = T)

myplots <- list()
for(species in spVec){
  myplots[[species]] <- colrNigon(species, 5e5)
}

p2 <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]],
                myplots[[5]], sxP,
                labels =  LETTERS[10:15], label_size = 12, nrow = 2)

ggsave("p2_t.pdf", p2, width = 21, height = 29.7, units = "cm")




## third page

# First get Ascaris
spVec <- grep("scari", unique(nigonDict_allSp$assembly),
              value = T)

myplots <- list()
for(species in spVec){
  myplots[[species]] <- colrNigon(species, 5e5)
}

spVec <- grep("hete|diplo", unique(nigonDict_allSp$assembly),
              value = T)

for(species in spVec){
  myplots[[species]] <- colrNigon(species, 1e5)
}

p3 <- plot_grid(myplots[[1]], myplots[[2]], myplots[[3]], 
                labels = LETTERS[16:19], label_size = 12, nrow = 1)

ggsave("p3_t.pdf", p3, width = 21, height = 29.7, units = "cm")


# "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o test.pdf *t.pdf





