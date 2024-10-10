#==load packages==
library(stringr)
library(Biostrings)
library(ggplot2)
library(tidyverse)
library(readxl)
library(lubridate)
library(ggpubr)
library(rworldmap)
library(ape)
library(ggtree)
library(phangorn)
library(phytools)
library("scales")
library(data.table)
library(lubridate)
library(treeio)
library(ggsci)

#==define color==
color <- pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)
color1 <- pal_npg("nrc", alpha =1)(10)
show_col(color1)
show_col(color)

#==read MCC tree (H3N2)==
mcc_tree <- read.beast("../data/tree_data/seq_h3n2_p1_geo_MCC.tre") 
mcc_tree@phylo$edge.length[mcc_tree@phylo$edge.length < 0] <- 0
n <- length(mcc_tree@data$geography[mcc_tree@data$geography == "Temperate_region+Southeast_Asia"])
set.seed(2000)
mcc_tree@data$geography[mcc_tree@data$geography == "Temperate_region+Southeast_Asia"] <- sample(c("Temperate_region","Southeast_Asia"), n, replace = T)

#==read metadata and persistent lineages (H3N2)==
h3n2_even_meta <- read.delim("../data/others/date_h3n2_p1.tsv") %>% 
  mutate(region = sapply(str_split(Name, "\\|"), function(x) x[4])) %>%
  mutate(Location = sapply(str_split(Name, "\\|"), function(x) x[5]))
h3n2_seq_type <- read.csv("../data/others/seq_h3n2_p1_seq_type.csv")
h3n2_even_meta$region[h3n2_even_meta$Name %in% h3n2_seq_type$tip[h3n2_seq_type$type == "Persistent lineages"]] <- "Southeastern Asia\n(Persistent lineages)"
h3n2_even_meta$region[h3n2_even_meta$Name %in% h3n2_seq_type$tip[h3n2_seq_type$type == "Non-persistent lineages"]] <- "Southeastern Asia\n(Non-persistent lineages)"
h3n2_even_meta$region[h3n2_even_meta$region == "Temperate_region"] <- "Temperate regions"
names(h3n2_even_meta)[1] <- "taxa"
factor(h3n2_even_meta$region, levels = unique(h3n2_even_meta$region)[c(3,2,1)]) -> h3n2_even_meta$region

#==plot (H3N2)==
ggtree(mcc_tree, mrsd = as.Date("2023-12-13"), aes(color = geography), size = 0.5) %<+% h3n2_even_meta[,c(1,3)] + geom_tippoint(aes(fill = region),size=1.2, color='transparent',shape=21, stroke=0.1)+
  annotate("rect", xmin = decimal_date(as.Date("2009-07-01")), xmax = decimal_date(as.Date("2010-06-30")),
           ymin = -30, ymax = 4400,  alpha = 0.3, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = decimal_date(as.Date("2020-07-01")), xmax = decimal_date(as.Date("2021-06-30")),
           ymin = -30, ymax = 4400, alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  theme_tree2()+
  scale_x_continuous(breaks = seq(2004,2024,2),limits = c(2004,2024), expand = c(0.015,0))+
  scale_y_continuous(limits = c(-50, 4400), expand = c(0,0))+
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_line())+
  scale_color_manual("Geographic regions\nalong the branch", 
                     labels = c("Southeastern Asia",
                                "Temperate regions"), values = color1[c(1,4)])+
  scale_fill_manual("Regions of virus collection", values = c(color1[c(4)],color[c(5,6)])) -> p2

#==read MCC tree (BV)==
mcc_tree <- read.beast("../data/tree_data/seq_bv_p1_geo_MCC.tre") 
mcc_tree@phylo$edge.length[mcc_tree@phylo$edge.length < 0] <- 0
n <- length(mcc_tree@data$geography[mcc_tree@data$geography == "Temperate_region+Southeast_Asia"])
set.seed(2000)
mcc_tree@data$geography[mcc_tree@data$geography == "Temperate_region+Southeast_Asia"] <- sample(c("Temperate_region","Southeast_Asia"), n, replace = T)

#==read metadata and persistent lineages (BV)==
bv_even_meta <- read.delim("../data/others/date_bv_p1.tsv") %>% 
  mutate(region = sapply(str_split(Name, "\\|"), function(x) x[4])) %>%
  mutate(Location = sapply(str_split(Name, "\\|"), function(x) x[5]))
bv_seq_type <- read.csv("../data/others/seq_bv_p1_seq_type.csv")
bv_even_meta$region[bv_even_meta$Name %in% bv_seq_type$tip[bv_seq_type$type == "Persistent lineages"]] <- "Southeastern Asia (Persistent lineages)"
bv_even_meta$region[bv_even_meta$Name %in% bv_seq_type$tip[bv_seq_type$type == "Non-persistent lineages"]] <- "Southeastern Asia (Non-persistent lineages)"
bv_even_meta$region[bv_even_meta$region == "Temperate_region"] <- "Temperate regions"
names(bv_even_meta)[1] <- "taxa"
factor(bv_even_meta$region, levels = unique(bv_even_meta$region)[c(3,2,1)]) -> bv_even_meta$region

#==plot (BV)==
ggtree(mcc_tree, mrsd = as.Date("2023-12-15"), aes(color = geography), size = 0.5)  %<+% bv_even_meta[,c(1,3)] + geom_tippoint(aes(fill = region),size=1.2, color='transparent',shape=21, stroke=0.1)+
  annotate("rect", xmin = decimal_date(as.Date("2009-07-01")), xmax = decimal_date(as.Date("2010-06-30")),
           ymin = -30, ymax = 4300, alpha = 0.3, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = decimal_date(as.Date("2020-07-01")), xmax = decimal_date(as.Date("2021-06-30")),
           ymin = -30, ymax = 4300, alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  theme_tree2()+
  scale_x_continuous(breaks = seq(2004,2024,2),limits = c(2004,2024), expand = c(0.015,0))+
  scale_y_continuous(limits = c(-50, 4300), expand = c(0,0))+
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = c(0.3,0.7),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    # legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_line())+
  scale_color_manual("Geographic regions along the branch", 
                     labels = c("Southeastern Asia",
                                "Temperate regions"), values = color1[c(1,4)])+
  guides(color = guide_legend(order = 2),
         fill = guide_legend(order = 1))+
  scale_fill_manual("Regions of virus collection", values = c(color1[c(4)],color[c(5,6)]))-> p3

#==proportion of persistent lineages==
#1. H3N2
h3n2_seq_type$Date <- as.Date(h3n2_seq_type$Date)
cut(h3n2_seq_type$Date, breaks = c(as.Date("2006-12-31"), as.Date("2007-06-30"),
                                   as.Date("2008-06-30"), as.Date("2009-06-30"),
                                   as.Date("2010-06-30"), as.Date("2011-06-30"),
                                   as.Date("2012-06-30"), as.Date("2013-06-30"),
                                   as.Date("2014-06-30"), as.Date("2015-06-30"),
                                   as.Date("2016-06-30"), as.Date("2017-06-30"),
                                   as.Date("2018-06-30"), as.Date("2019-06-30"),
                                   as.Date("2020-06-30"), as.Date("2021-06-30"),
                                   as.Date("2022-06-30"), as.Date("2023-06-30"),
                                   as.Date("2024-06-30")),
    labels = c("2006/2007","2007/2008","2008/2009","2009/2010","2010/2011","2011/2012","2012/2013","2013/2014","2014/2015",
               "2015/2016","2016/2017","2017/2018","2018/2019","2019/2020","2020/2021","2021/2022","2022/2023","2023/2024")) -> h3n2_seq_type$season

h3n2_dis <- h3n2_seq_type %>% group_by(season,type) %>% summarise(n = n()) %>%
  ungroup() %>% group_by(season) %>% mutate(n_total = sum(n))

factor(h3n2_dis$type, levels = unique(h3n2_dis$type)[2:1]) -> h3n2_dis$type
ggplot(h3n2_dis)+
  annotate("rect", xmin = 3.5, xmax = 4.5,
           ymin = -5, ymax = 310,alpha = 0.3, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = 14.5, xmax = 15.5,
           ymin = -5, ymax = 310,alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  geom_bar(aes(x = season, y = n, fill = type), stat = "identity", width = 0.6)+
  geom_line(h3n2_dis[h3n2_dis$type == "Persistent lineages",], mapping = aes(x = season, y = n*300/n_total, group = 1), color = "grey50")+
  geom_point(h3n2_dis[h3n2_dis$type == "Persistent lineages",], mapping = aes(x = season, y = n*300/n_total), color = "grey50")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%b %Y", expand = c(0.01,0))+
  scale_x_discrete( expand = c(0.03,0))+
  theme_bw()+
  scale_y_continuous(limits = c(-5, 310), expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/300,
                                         name="Proportion of persistent lineages (%)",
                                         # breaks = seq(0,0.15,0.05),
                                         labels = seq(0,100,25)
                     ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank())+
  theme(legend.position = "none",legend.background = element_rect(fill = "transparent", color = "transparent"))+
  scale_fill_manual("Types of lineage\n(Southeastern Asia)", values = color[c(5,6)])+
  labs(x = "Influenza seasons", y = "Number of A/H3N2 sequences",  tag = "d")-> p4

#2. B/V
bv_seq_type$Date <- as.Date(bv_seq_type$Date)
cut(bv_seq_type$Date, breaks = c(as.Date("2006-12-31"), as.Date("2007-06-30"),
                                 as.Date("2008-06-30"), as.Date("2009-06-30"),
                                 as.Date("2010-06-30"), as.Date("2011-06-30"),
                                 as.Date("2012-06-30"), as.Date("2013-06-30"),
                                 as.Date("2014-06-30"), as.Date("2015-06-30"),
                                 as.Date("2016-06-30"), as.Date("2017-06-30"),
                                 as.Date("2018-06-30"), as.Date("2019-06-30"),
                                 as.Date("2020-06-30"), as.Date("2021-06-30"),
                                 as.Date("2022-06-30"), as.Date("2023-06-30"),
                                 as.Date("2024-06-30")),
    labels = c("2006/2007","2007/2008","2008/2009","2009/2010","2010/2011","2011/2012","2012/2013","2013/2014","2014/2015",
               "2015/2016","2016/2017","2017/2018","2018/2019","2019/2020","2020/2021","2021/2022","2022/2023","2023/2024")) -> bv_seq_type$season

bv_dis <- bv_seq_type %>% group_by(season,type) %>% summarise(n = n()) %>%
  ungroup() %>% group_by(season) %>% mutate(n_total = sum(n))
factor(bv_dis$type, levels = unique(bv_dis$type)[2:1]) -> bv_dis$type

ggplot(bv_dis)+
  annotate("rect", xmin = 3.5, xmax = 4.5,
           ymin = -5, ymax = 310,alpha = 0.3, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = 14.5, xmax = 15.5,
           ymin = -5, ymax = 310,alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  geom_bar(aes(x = season, y = n, fill = type), stat = "identity", width = 0.6)+
  geom_line(bv_dis[bv_dis$type == "Non-persistent lineages",], mapping = aes(x = season, y = (1-n/n_total)*300, group = 1), color = "grey50")+
  geom_point(bv_dis[bv_dis$type == "Non-persistent lineages",], mapping = aes(x = season, y = (1-n/n_total)*300), color = "grey50")+
  
  # scale_x_date(date_breaks = "1 year", date_labels = "%b %Y", expand = c(0.01,0))+
  scale_x_discrete(expand = c(0.03,0))+
  theme_bw()+
  scale_y_continuous(limits = c(-5, 310), expand = c(0,0),
                     sec.axis = sec_axis(trans=~.* 1/300,
                                         name="Proportion of persistent lineages (%)",
                                         # breaks = seq(0,0.15,0.05),
                                         labels = seq(0,100,25)
                     ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  theme(legend.position = c(0.3,0.845),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        panel.grid = element_blank(),
        # legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  scale_fill_manual("Types of lineage (Southeastern Asia)", values = color[c(5,6)])+
  labs(x = "Influenza seasons", y = "Number of B/Victoria sequences", tag = "e")-> p5

(p4|p5)+plot_layout(guides = "collect")
