library(reshape2)
library(tidyverse)
library(scales)
library(gtools)
library(ggpubr)

# Function to calculate neighboring distances
calcNeighDist <- function(stCoord, eCoord){
  # combine start and end
  # in a single column to be able to use the diff function
  gSpOdds <- ifelse(1:length(stCoord)%%2, stCoord, eCoord)
  gSpEven <- ifelse(!1:length(stCoord)%%2, stCoord, eCoord)
  neigDifOdds <- c(0, diff(gSpOdds))
  neigDifEven <- c(0, diff(gSpEven))
  neigDif <- ifelse(1:length(stCoord)%%2, neigDifOdds, neigDifEven)
  return(neigDif)
}

# Read input files
otiDict <- read_tsv("analyses/orthoFinder/longMar06/otipu_celeChr.dict.tsv")
rbbh <- read_tsv("analyses/orthoFinder/longMar06/rbbh_otipu_vs_celeg.tsv.gz")
ceOperons <- read_tsv("analyses/orthoFinder/celeNigons/cele_gID_gName_operon.tsv",
                      col_names = c("gID", "gName", "oID"))
geneDict <- read_csv("analyses/orthoFinder/caenorhabditis_elegans.PRJNA13758.WBPS14.csv.gz",
                     col_names = c("tID", "gID", "gName", "chr", "type")) %>%
  mutate(gID = sub("gene=Gene:", "", gID)) %>%
  left_join(ceOperons, by = "gID") %>%
  select(gID, tID, oID)



# Discard genes which are part of an operon
uniOprnGenes <- filter(allGpos, is.na(oID)) %>%
  select(-oID, -gID)

# Calculate distances between pair of genes in C. elegans and
# the distance between the corresponding orthologous pair
# in O. tipulae
pairCe_Ot_dist <- group_by(uniOprnGenes, chr.ce) %>%
  arrange(stPos.ce) %>%
  mutate(gRank = 1:n(),
         neigDif.ce = calcNeighDist(stPos.ce, stPos.ce),
         neigDif.ot = calcNeighDist(stPos.ot, stPos.ot)) %>%
  ungroup()

### Plot correspondence between chr I ###
ceOtplot <- filter(pairCe_Ot_dist, chr.ce == chr.ot) %>%
mutate(
       # arms positios from Table S6 of
       # 10.1186/gb-2010-11-12-r120
       inArms = ifelse(stPos.ce < 4856294, "left", ifelse(stPos.ce > 10127632, "right", "center"))
) %>%
  select(-celeg, -chr.ce, -neigDif.ce, -neigDif.ot) %>%
  melt(id.vars = c("otipu", "chr.ot", "gRank", "inArms")) %>%
  mutate(spPos = factor(sub("stPos.(.+)", "\\1", variable),
                        levels = c("ot", "ce")),
         variable = sub("(stPos).+", "\\1", variable) , value = ifelse(spPos == "ce", 15072434/2 - (15072434 - value), 8900181/2 - value)
  ) %>% filter(variable == "stPos", chr.ot == "I") %>%
  slice(sample(n())) %>%
  ggplot(aes(x = value, y = spPos,
             group = otipu, colour = gRank)) + 
  geom_line(alpha = 0.3) +
  scale_x_continuous(labels = label_number_si()) +
  scale_y_discrete(breaks = c("ot", "ce"),
                   labels = c("O. tipulae", "C. elegans")) +
  scale_color_gradientn(colours = rainbow(5)) +
  coord_cartesian(expand = 0) +
  # xlab("Distance from center of chromosome I") +
  xlab("") + 
  ylab("") +
  theme_light() +
  theme(plot.margin = margin(1, 0.5, 0.5, 0, "cm"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


otCeplot <- filter(pairCe_Ot_dist, chr.ce == chr.ot, chr.ot == "I") %>%
  mutate(
    gRank = rank(-stPos.ot),
    inArms = ifelse(stPos.ce < 4856294, "left", ifelse(stPos.ce > 10127632, "right", "center"))) %>%
  select(-celeg, -chr.ce, -neigDif.ce, -neigDif.ot) %>%
  melt(id.vars = c("otipu", "chr.ot", "gRank", "inArms")) %>%
  mutate(spPos = factor(sub("stPos.(.+)", "\\1", variable),
                        levels = c("ot", "ce")),
         gRank = as.integer(as.character(gRank)),
         variable = sub("(stPos).+", "\\1", variable) , value = ifelse(spPos == "ce", 15072434/2 - (15072434 - value), 8900181/2 - value)
  ) %>% filter(variable == "stPos", chr.ot == "I") %>%
  sample_frac(1) %>%
  ggplot(aes(x = value, y = spPos,
             group = otipu, colour = gRank)) + 
  geom_line(alpha = 0.3) +
  scale_x_continuous(labels = label_number_si()) +
  scale_y_discrete(breaks = c("ce", "ot"),
                   labels = c("C. elegans", "O. tipulae"),
                   limits = c("ce", "ot")) +
  scale_color_gradientn(colours = rainbow(5)) +
  coord_cartesian(expand = 0) +
  ylab("") +
  xlab("Distance from center of chromosome I") +
  theme_light() +
  theme(plot.margin = margin(0, 0.5, 0.5, 0, "cm"),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


ggarrange(ceOtplot, otCeplot,
          ncol = 1, nrow = 2)

ggsave(paste0("report/figures/celegans_vs_otipu_chrI.png"), width = 9, height = 6)

### Plot correlation of distance per chromosome  ###

lenCutoff <- 5e4
otCutoff <- 11e6
my_x_title <- expression(paste("Distance between adjacent genes in ", italic("C. elegans")))
my_y_title <- expression(paste("Distance between corresponding rbbh in ", italic("O. tipulae")))
filter(pairCe_Ot_dist, neigDif.ce < lenCutoff, abs(neigDif.ot) < otCutoff) %>%
  ggplot(aes(x=neigDif.ce, y=abs(neigDif.ot))) +
  facet_wrap(chr.ot ~ ., ncol = 3) +
  scale_x_continuous(labels = label_number_si()) +
  scale_y_continuous(labels = label_number_si()) +
  geom_point(alpha = 0.3) +
  geom_density_2d() +
  geom_text(data = count(pairCe_Ot_dist, chr.ot) %>%
              mutate(n = paste0("n=", n)),
            aes(x=lenCutoff*0.65, y=otCutoff*0.8, label=n), 
            colour="black", inherit.aes=FALSE,
            hjust = 0,
            parse=FALSE, size = 4) +
  geom_text(data = filter(pairCe_Ot_dist, neigDif.ce > lenCutoff | abs(neigDif.ot) > otCutoff) %>%
              count(chr.ot) %>%
              mutate(n = paste0("m=", n)),
            aes(x=lenCutoff*0.65, y=otCutoff*0.7, label=n), 
            colour="black", inherit.aes=FALSE,
            hjust = 0,
            parse=FALSE, size = 4) +
  geom_text(data = group_by(pairCe_Ot_dist, chr.ot) %>%
              summarise(r=cor(neigDif.ce, abs(neigDif.ot), method = "spearman")) %>%
              mutate(r = paste0("R=", round(r, 2))),
            aes(x=lenCutoff * 0.65, y=otCutoff*0.9, label=r), 
            colour="black", inherit.aes=FALSE,
            hjust = 0,
            parse=FALSE, size = 4) +
  ylab(my_y_title) + 
  xlab(my_x_title) +
  theme_light() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill="white"),
        strip.text.x = element_text(color = "black", face = "bold"))

ggsave(paste0("report/figures/otipuNeighborDistance_vs_celegansNeighborDistance.png"), width = 8, height = 5)
