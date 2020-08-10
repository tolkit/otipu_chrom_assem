library(reshape2)
library(tidyverse)
library(scales)
library(gtools)
library(ggpubr)


# Otipulae features
otiFeats <- read_tsv("analyses/annotation/repeats/concat_oscheius_tipulae.local.v3_3.1kw.tsv.gz",
                     col_names = c("chr", "start", "end", "feature")) %>%
  mutate(feature = case_when(
    feature == "trf" ~ "Tandem repeats",
    feature == "irf" ~ "Inverted repeats",
    feature == "repeatmasker" ~ "All RE",
    feature == "coding_introns" ~ "Coding introns",
    feature == "coding_exons" ~ "Coding exons",
    TRUE ~ as.character(feature)
  ))

####
# Single plot
####

cols <- c("ce_I" = "#bebebe", "ce_II" = "#cda662",
          "ce_III" = "#77ba6d", "ce_IV" = "#00c1b1",
          "ce_V" = "#6fb1e1", "ce_X" = "#D795DC",
          "NgA" = "#af0e2b", "NgB" = "#e4501e",
          "NgC" = "#4caae5", "NgD" = "#f3ac2c",
          "NgE" = "#57b741", "NgN" = "#8880be",
          "NgX" = "#81008b", "other" = "#706E6E")

ceChroms <- c("I", "II",
              "III", "IV", "V", "X")

featLevels <- c("Tandem repeats", "Inverted repeats",
                "All RE", "Coding introns",
                "Coding exons",
                "C. elegans chr", "Nigon")

windwSize <- 5e5
mutate(otiFeats, ort = ifelse(feature %in% ceChroms,
                              paste0("ce_", feature),
                              ifelse(grepl("Ng", feature),
                                     feature, "other")),
       feature = ifelse(feature %in% ceChroms, "C. elegans chr",
                        ifelse(grepl("Ng", feature), "Nigon", 
                               feature))) %>%
  filter(feature %in% featLevels) %>%
  mutate(feature = factor(feature, levels = featLevels)) %>%
  group_by(chr) %>%
  mutate(ints = as.numeric(as.character(cut(start,
                                            breaks = seq(0, max(start), windwSize),
                                            labels = seq(windwSize, max(start), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  ungroup() %>%
  group_by(chr, feature, ort, ints) %>%
  summarise(rdens = sum(end - start) / windwSize) %>%
  ungroup() %>%
  ggplot(aes(x=ints, y=rdens, fill = ort)) + 
  facet_grid(chr ~ feature, scales = "free", switch = "y") +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  coord_flip() +
  # xlab(paste0("Feature density per ", windwSize/1e3, "kbs")) +
  # ylab(my_x_title) +
  theme_minimal() +
  scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, NA)) +
  scale_x_continuous(labels = label_number_si(), position = "top") +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(nrow = 2,
                             title = element_blank())) +
  theme(axis.title.y=element_blank(),
        panel.border = element_blank(),
        legend.position="bottom",
        panel.spacing.x = unit(1.1, "lines"),
        panel.spacing.y = unit(0.5, "lines")
  )


ggsave("report/figures/extras/otipu3.3_Fig2.png", width = 9, height = 10.5)
