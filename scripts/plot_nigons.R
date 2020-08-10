# Plot stacked barplots per window
# Input is chr, start, refChr
library(tidyverse)
library(scales)
library(gtools)


nigonDict <- read_tsv("analyses/orthoFinder/gene2Nigon20200528_09_55.tsv.gz") %>%
    select(-kinfCut)

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")


# Plot Nigons per assembly
for(species in unique(nigonDict$assembly)){
  plotHeight <- 6
  if(species == "meloidogyne_hapla.PRJNA29083.WBPS14")
    next
  if(species == "ascaris_suum.local.v2020.genomic")
    plotHeight <- 20
  spName <- sub("\\..+", "", species)
  spName <- sub("_", " ", spName)
  spName <- paste(toupper(substr(spName, 1, 1)), substr(spName, 2, nchar(spName)), sep="")
  # Bin size
  windwSize <- 5e5
  
  filter(nigonDict, assembly == species) %>%
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
    geom_bar(position="stack", stat="identity"#,
             # width = windwSize, colour="black",
             # size=0.2
             ) +
    # coord_cartesian(xlim = c(0,11170816+5e4)) +
    theme_minimal() +
    scale_y_continuous(position = "right") +
    scale_x_continuous(labels = label_number_si()) +
    scale_fill_manual(values = cols) +
    ggtitle(paste0("*", spName, "*")) +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          panel.border = element_blank(),
          plot.title = ggtext::element_markdown())
  
  ggsave(paste0("report/figures/extras/chrom2Nigon/", species,
                ".png"), width = 5, height = plotHeight)
}

