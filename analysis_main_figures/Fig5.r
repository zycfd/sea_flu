#==load packages==
library(stringr)
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
library(sf)
library(rnaturalearth)
library(patchwork)
library(grid)

#==define color==
colors <- c(pal_npg("nrc", alpha =1)(10)[c(1:7,9:10)],"darkred","#FADDA9","grey90")
colors1 <- c(pal_aaas("default", alpha = 0.3)(10)[10:1])
colors2 <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
             "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
             "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
             "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
             "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
             "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
colors3 <- pal_lancet(palette = c("lanonc"), alpha = 1)(9)
show_col(colors3)

#==define order of sub-location==
levels = c("North China", "Central China", "South China", "Hong Kong", "Hong Kong, China","Macao","Macao, China", "Taiwan","Taiwan, China",
           "Japan", "South Korea", "Bangladesh",  "Bhutan","Nepal", "India","Sri Lanka" ,"Myanmar",
           "Laos",    "Philippines", "Thailand","Cambodia", "Vietnam",
           "Malaysia-Brunei", "Singapore", "Indonesia-East Timor")
levels1 = c("NorthChina", "CentralChina", "SouthChina", "HongkongChina","MacaoChina", "TaiwanChina",
            "Japan", "SouthKorea", "Bangladesh",  "Bhutan","Nepal", "India","SriLanka" ,"Myanmar",
            "Laos",    "Philippines", "Thailand","Cambodia", "Vietnam",
            "MalaysiaBrunei", "Singapore", "IndonesiaEastTimor")

#==read data (h3n2)==
h3n2_jump <- rbind(read.delim("../phylogeographic_analysis/jump_SEA/h3n2_p1_sea_Jumps_0_2.5M.txt", sep = ","),
                   read.delim("../phylogeographic_analysis/jump_SEA/h3n2_p1_sea_Jumps_2.5_5M.txt", sep = ","))
h3n2_jump$date <- decimal2Date(decimal_date(as.Date("2023-11-27")) - h3n2_jump$time)
h3n2_jump$period[h3n2_jump$date < as.Date("2010-07-01") & h3n2_jump$date >= as.Date("2009-07-01")] <- "Period 1"
h3n2_jump$period[h3n2_jump$date < as.Date("2020-07-01") & h3n2_jump$date >= as.Date("2010-07-01")] <- "Period 2"
h3n2_jump$period[h3n2_jump$date < as.Date("2021-07-01") & h3n2_jump$date >= as.Date("2020-07-01")] <- "Period 3"
h3n2_jump1 <- h3n2_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())

length(unique(h3n2_jump1$treeId))
length(unique(h3n2_jump1$startLocation))
length(unique(h3n2_jump1$endLocation))

h3n2_tmp <- rbind(data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation), 91), each = 22)) %>% mutate(period = "Period 1"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation), 91), each = 22)) %>% mutate(period = "Period 2"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "Period 3"))
h3n2_jump2 <- left_join(h3n2_tmp, h3n2_jump1)
nrow(h3n2_jump2[!is.na(h3n2_jump2$n),])
h3n2_jump2$n[is.na(h3n2_jump2$n)] <- 0

h3n2_jump2 <- h3n2_jump2 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n),
            sd = mean(n))

h3n2_jump2$startLocation[h3n2_jump2$startLocation == "CentralChina"] <- "Central China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "HongkongChina"] <- "Hong Kong, China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "IndonesiaEastTimor"] <- "Indonesia-East Timor"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "MacaoChina"] <- "Macao, China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "MalaysiaBrunei"] <- "Malaysia-Brunei"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "NorthChina"] <- "North China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthChina"] <- "South China"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SouthKorea"] <- "South Korea"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "SriLanka"] <- "Sri Lanka"
h3n2_jump2$startLocation[h3n2_jump2$startLocation == "TaiwanChina"] <- "Taiwan, China"

h3n2_jump2$endLocation[h3n2_jump2$endLocation == "CentralChina"] <- "Central China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "HongkongChina"] <- "Hong Kong, China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "IndonesiaEastTimor"] <- "Indonesia-East Timor"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "MacaoChina"] <- "Macao, China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "MalaysiaBrunei"] <- "Malaysia-Brunei"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "NorthChina"] <- "North China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthChina"] <- "South China"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SouthKorea"] <- "South Korea"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "SriLanka"] <- "Sri Lanka"
h3n2_jump2$endLocation[h3n2_jump2$endLocation == "TaiwanChina"] <- "Taiwan, China"

h3n2_jump2$num_season[h3n2_jump2$period == "Period 1"] <- 1
h3n2_jump2$num_season[h3n2_jump2$period == "Period 2"] <- 10
h3n2_jump2$num_season[h3n2_jump2$period == "Period 3"] <- 1
h3n2_jump2$mean_per_senson <- h3n2_jump2$mean / h3n2_jump2$num_season

range(h3n2_jump2$mean_per_senson)
cut(h3n2_jump2$mean_per_senson, breaks = c(0, 0.5, 1, 1.5, 2, 3.5),right = T,
    labels = c("(0, 0.5]","(0.5, 1]","(1, 1.5]","(1.5, 2]","(2, 3.5]")) -> h3n2_jump2$mean_per_senson1
table(h3n2_jump2$mean_per_senson1)

#==read map data==
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  mutate(adm0_a3 = ifelse(adm0_a3 == "SDS", "SSD" ,adm0_a3)) %>% 
  # filter(!sovereignt %in% c("China", "Taiwan")) %>%
  filter(admin != "Antarctica")
world <- world[,c(4,5,10,11,169)]
worldmap2 <- world[(world$adm0_a3 %in%
                      c("BGD","BTN","IND","NPL","LKA","CHN","HKG","MAC","TWN","JPN","KOR",
                        "BRN", "KHM", "IDN", "LAO", "MYS", "MMR", "PHL", "SGP", "THA", "TLS", "VNM")),]
worldmap2 <- st_transform(worldmap2, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

map <- worldmap2
map$region_final1 <- NA
map$region_final1[map$region_final %in% c("Japan","South Korea", "Taiwan", "China")] <- "East Asia"
map$region_final1[map$region_final %in% c("Bangladesh","Bhutan","India","Nepal","Sri Lanka")] <- "South Asia"
map$region_final1[is.na(map$region_final1)] <- "Southeast Asia"

#==read geo data==
tmp <- readRDS("../data/geo_data/lati_long.rds")
tmp$geometry <- as.character(tmp$geometry)
tmp$lati <- str_replace_all(tmp$geometry,"c|\\(|\\)","")
tmp$lati <- as.numeric(str_trim(sapply(str_split(tmp$lati,","),function(x) x[2])))
tmp$long <- str_replace_all(tmp$geometry,"c|\\(|\\)","")
tmp$long <- as.numeric(str_trim(sapply(str_split(tmp$long ,","),function(x) x[1])))

tmp$long[tmp$region_final == "Japan"] <- 140
tmp$long[tmp$region_final == "South China"] <- 117
tmp$long[tmp$region_final == "Vietnam"] <- 107
tmp$long[tmp$region_final == "Malaysia-Brunei"] <- 115
tmp$long[tmp$region_final == "Indonesia-East Timor"] <- 120

h3n2_jump3 <- h3n2_jump2 %>%
  left_join(tmp, by = c("startLocation" = "region_final"))%>%
  left_join(tmp, by = c("endLocation" = "region_final"))

h3n2_jump4 <- left_join(h3n2_jump3[h3n2_jump3$period %in% c("Period 1", "Period 3"),], 
                        h3n2_jump3[h3n2_jump3$period %in% c("Period 2"),] %>% ungroup() %>% select("startLocation", "endLocation", "mean_per_senson") %>% 
                          rename("mean_per_senson_interpandemic" = mean_per_senson) )

#==map (panel a-c)==
#Map files in China at province-resolution can be downloaded from https://github.com/GaryBikini/ChinaAdminDivisonSHP
ggplot() + 
  geom_sf(data = map, aes(fill = region_final1) ,color = "black", size = 0.02)+
  geom_curve(data = h3n2_jump3[h3n2_jump3$period == "Period 2" & h3n2_jump3$mean_per_senson >= 0.5,],
             aes(x = as.double(long.x), 
                 y = as.double(lati.x), 
                 xend = as.double(long.y), 
                 yend = as.double(lati.y),
                 linewidth = mean_per_senson
             ),
             color = "grey50",
             arrow = arrow(angle = 40, length = unit(0.1, "inches"),
                           ends = "last", type = "open"),
             alpha = 0.5,
             curvature = 0.35)+
  # guides(fill = F, linewidth = F)+
  # guides(fill = F)+
  guides(linewidth = guide_legend(nrow = 1, title.position = "top"),
         fill = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = data.frame(tmp), aes(x = long, y = lati),shape = 21, size = 3,fill = "grey50")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  # scale_x_continuous(limits=c(-170,170))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Viral movement between each pair of sub-locations per season", range = c(0.5,3.5), 
                             # breaks = seq(1,5,1),
                             # labels = c("(0.1, 0.5]","(0.5, 1.0]","(1.0, 2.0]","(2.0, 5.0]","(5.0, 10.0]"), 
                             limits = c(0.5,3.05)
  )+
  scale_fill_manual("Sub-regions within Southeastern Asia", values = colors3[c(4,5,6)])+
  labs(subtitle = "a. During interpandemic period")-> p1

ggplot() +
  geom_sf(data = map, aes(fill = region_final1) ,color = "black", size = 0.02)+
  geom_curve(data = h3n2_jump3[h3n2_jump3$period == "Period 1" & h3n2_jump3$mean_per_senson >= 0.5,],
             aes(x = as.double(long.x),
                 y = as.double(lati.x),
                 xend = as.double(long.y),
                 yend = as.double(lati.y),
                 linewidth = mean_per_senson
             ),
             color = "grey50",
             arrow = arrow(angle = 40, length = unit(0.1, "inches"),
                           ends = "last", type = "open"),
             alpha = 0.5,
             curvature = 0.35)+
  # guides(fill = F, linewidth = F)+
  # guides(fill = F)+
  guides(linewidth = guide_legend(nrow = 1, title.position = "top"),
         fill = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = data.frame(tmp), aes(x = long, y = lati),shape = 21, size = 3,fill = "grey50")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  # scale_x_continuous(limits=c(-170,170))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Viral movement between each pair of sub-locations per season", range = c(0.5,3.5),
                             # breaks = seq(1,5,1),
                             # labels = c("(0.1, 0.5]","(0.5, 1.0]","(1.0, 2.0]","(2.0, 5.0]","(5.0, 10.0]"),
                             limits = c(0.5,3.05)
  )+
  scale_fill_manual("Sub-regions within Southeastern Asia", values = colors3[c(4,5,6)])+
  labs(subtitle = "b. During A/H1N1 pandemic season")-> p2

ggplot() +
  geom_sf(data = map, aes(fill = region_final1) ,color = "black", size = 0.02)+
  geom_curve(data = h3n2_jump3[h3n2_jump3$period == "Period 3" & h3n2_jump3$mean_per_senson >= 0.5,],
             aes(x = as.double(long.x),
                 y = as.double(lati.x),
                 xend = as.double(long.y),
                 yend = as.double(lati.y),
                 linewidth = mean_per_senson
             ),
             color = "grey50",
             arrow = arrow(angle = 40, length = unit(0.1, "inches"),
                           ends = "last", type = "open"),
             alpha = 0.5,
             curvature = 0.35)+
  # guides(fill = F, linewidth = F)+
  # guides(fill = F)+
  guides(linewidth = guide_legend(nrow = 1, title.position = "top"),
         fill = guide_legend(nrow = 1, title.position = "top"))+
  theme_void()+
  geom_point(data = data.frame(tmp), aes(x = long, y = lati),shape = 21, size = 3,fill = "grey50")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+
  # scale_x_continuous(limits=c(-170,170))+
  theme(legend.direction = 'horizontal')+
  scale_linewidth_continuous("Viral movement between each pair of sub-locations per season", range = c(0.5,3.5),
                             # breaks = seq(1,5,1),
                             # labels = c("(0.1, 0.5]","(0.5, 1.0]","(1.0, 2.0]","(2.0, 5.0]","(5.0, 10.0]"),
                             limits = c(0.5,3.05)
  )+
  scale_fill_manual("Sub-regions within Southeastern Asia", values = colors3[c(4,5,6)])+
  labs(subtitle = "c. During COVID-19 pandemic season")-> p3

#==plot (panel f)==
phylo_simi <- read.csv("../data/others/phylo_similarity.csv")

ggplot(data = phylo_simi) +
  annotate("rect", xmin = 2.5, xmax = 3.5,
           ymin = -0.01, ymax = 0.5, alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = 13.5, xmax = 14.5,
           ymin = -0.01, ymax = 0.5, alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  geom_ribbon(aes(season_type, ymin = hpd_low, ymax = hpd_upp, fill = compare_location, group = compare_location), alpha = 0.2)+
  geom_line(aes(x = season_type, y = mean, color = compare_location, group = compare_location))+
  geom_point(aes(x = season_type, y = mean, color = compare_location, group = compare_location))+
  scale_x_discrete(expand = c(0.01,0))+
  scale_y_continuous(limits = c(-0.01,0.5), expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(x = "Influenza seasons", y = "Pairwise PhyloSor similarity", tag = "f")+
  scale_fill_manual("Pairs of locations within Southeastern Asia",values = colors3[c(2,1,5)])+
  scale_color_manual("Pairs of locations within Southeastern Asia",values = colors3[c(2,1,5)])+
  guides(fill = guide_legend(title.position = "top"))-> p2_tmp

#==plot (panel e)==
h3n2_jump <- rbind(read.delim("../phylogeographic_analysis/jump_SEA/h3n2_p1_sea_Jumps_0_2.5M.txt", sep = ","),
                   read.delim("../phylogeographic_analysis/jump_SEA/h3n2_p1_sea_Jumps_2.5_5M.txt", sep = ","))
h3n2_jump$date <- decimal2Date(decimal_date(as.Date("2023-11-27")) - h3n2_jump$time)
h3n2_jump$period[h3n2_jump$date < as.Date("2010-07-01") & h3n2_jump$date >= as.Date("2009-07-01")] <- "2009/2010"
h3n2_jump$period[h3n2_jump$date < as.Date("2011-07-01") & h3n2_jump$date >= as.Date("2010-07-01")] <- "2010/2011"
h3n2_jump$period[h3n2_jump$date < as.Date("2012-07-01") & h3n2_jump$date >= as.Date("2011-07-01")] <- "2011/2012"
h3n2_jump$period[h3n2_jump$date < as.Date("2013-07-01") & h3n2_jump$date >= as.Date("2012-07-01")] <- "2012/2013"
h3n2_jump$period[h3n2_jump$date < as.Date("2014-07-01") & h3n2_jump$date >= as.Date("2013-07-01")] <- "2013/2014"
h3n2_jump$period[h3n2_jump$date < as.Date("2015-07-01") & h3n2_jump$date >= as.Date("2014-07-01")] <- "2014/2015"
h3n2_jump$period[h3n2_jump$date < as.Date("2016-07-01") & h3n2_jump$date >= as.Date("2015-07-01")] <- "2015/2016"
h3n2_jump$period[h3n2_jump$date < as.Date("2017-07-01") & h3n2_jump$date >= as.Date("2016-07-01")] <- "2016/2017"
h3n2_jump$period[h3n2_jump$date < as.Date("2018-07-01") & h3n2_jump$date >= as.Date("2017-07-01")] <- "2017/2018"
h3n2_jump$period[h3n2_jump$date < as.Date("2019-07-01") & h3n2_jump$date >= as.Date("2018-07-01")] <- "2018/2019"
h3n2_jump$period[h3n2_jump$date < as.Date("2020-07-01") & h3n2_jump$date >= as.Date("2019-07-01")] <- "2019/2020"
h3n2_jump$period[h3n2_jump$date < as.Date("2021-07-01") & h3n2_jump$date >= as.Date("2020-07-01")] <- "2020/2021"

h3n2_jump1 <- h3n2_jump %>%
  filter(!is.na(period)) %>%
  group_by(treeId, period, startLocation, endLocation) %>%
  summarise(n = n())

length(unique(h3n2_jump1$treeId))
length(unique(h3n2_jump1$startLocation))
length(unique(h3n2_jump1$endLocation))

h3n2_tmp <- rbind(data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2009/2010"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2010/2011"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2011/2012"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2012/2013"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2013/2014"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2014/2015"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2015/2016"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2016/2017"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2017/2018"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2018/2019"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2019/2020"),
                  data.frame(treeId = rep(unique(h3n2_jump1$treeId), each = 22*22),
                             startLocation = rep(unique(h3n2_jump1$startLocation), 22*91),
                             endLocation = rep(rep(unique(h3n2_jump1$endLocation),91), each = 22)) %>% mutate(period = "2020/2021"))

h3n2_jump2 <- left_join(h3n2_tmp, h3n2_jump1)
nrow(h3n2_jump2[!is.na(h3n2_jump2$n),])
h3n2_jump2$n[is.na(h3n2_jump2$n)] <- 0

h3n2_jump2 <- h3n2_jump2 %>%
  group_by(period, startLocation, endLocation) %>%
  summarise(mean = mean(n)) %>%
  filter(startLocation != endLocation) 

#==3 regions (directional)==
h3n2_jump4 <- h3n2_jump2
h3n2_jump4$startLocation2[h3n2_jump4$startLocation %in% c("TaiwanChina", "Japan", "NorthChina", "CentralChina","SouthChina", "HongkongChina", "SouthKorea","MacaoChina")] <- "East Asia"
h3n2_jump4$startLocation2[h3n2_jump4$startLocation %in% c("Bangladesh", "Bhutan","Nepal","SriLanka","India")] <- "South Asia"
h3n2_jump4$startLocation2[h3n2_jump4$startLocation %in% c("MalaysiaBrunei","Laos","Cambodia","Vietnam","Singapore", "Myanmar","Philippines","Thailand","IndonesiaEastTimor")] <- "Southeast Asia"
h3n2_jump4$endLocation2[h3n2_jump4$endLocation %in% c("TaiwanChina", "Japan", "NorthChina", "CentralChina","SouthChina", "HongkongChina", "SouthKorea","MacaoChina")] <- "East Asia"
h3n2_jump4$endLocation2[h3n2_jump4$endLocation %in% c("Bangladesh", "Bhutan","Nepal","SriLanka","India")] <- "South Asia"
h3n2_jump4$endLocation2[h3n2_jump4$endLocation %in% c("MalaysiaBrunei","Laos","Cambodia","Vietnam","Singapore", "Myanmar","Philippines","Thailand","IndonesiaEastTimor")] <- "Southeast Asia"

h3n2_jump5 <- h3n2_jump4 %>% 
  mutate(move_dir = paste0(startLocation2,"_",endLocation2)) %>%
  # filter(startLocation2 != endLocation2) %>%
  group_by(period, move_dir) %>%
  summarise(mean = sum(mean)) 

h3n2_jump5 <- as.data.frame(spread(h3n2_jump5, move_dir, mean))
row.names(h3n2_jump5) <- h3n2_jump5[,1]
h3n2_jump5 <- h3n2_jump5[,-1]
row.names(h3n2_jump5)[c(1,12)] <- c("During\nA/H1N1\npandemic\nseason", "During\nCOVID-19\npandemic\nseason")

d <- dist(h3n2_jump5)
fit <- cmdscale(d, eig = TRUE, k = 2)

# Extract (x, y) coordinates
x <- fit$points[, 1]
y <- fit$points[, 2]

# Create a scatter plot
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Multidimensional Scaling Results", type = "n")
text(x, y, labels = row.names(h3n2_jump5))

mds_data <- as.data.frame(fit$points)
mds_data$label <- row.names(mds_data)
ggplot(mds_data) +
  geom_point(aes(x = V1, y = V2),shape = 21, size = 3, fill = c(colors[1],rep(colors[2],10),colors[1]), color = "grey90") +
  ggrepel::geom_text_repel(aes(x = V1, y = V2, label = label), 
                           vjust = 0.8, size = 2.5)+
  theme_bw()+
  scale_x_continuous(limits = c(-65,35),breaks = seq(-65,35,20))+
  scale_y_continuous(limits = c(-65,35),breaks = seq(-65,35,20))+
  theme(
    axis.text = element_blank(),
    plot.tag = element_text(size = 12))+
  labs(tag = "e", x = "Dimension 1", y = "Dimension 2")-> p1_4

#==plot (panel d)==
h3n2_even <- data.frame(rbind(t(read.delim("../phylogeographic_analysis/log_file/h3n2_p1_sea_log.tsv"))))
h3n2_even$type <- row.names(h3n2_even)
h3n2_even <- h3n2_even[,c(1,8,12)]
colnames(h3n2_even) <- h3n2_even[1,]
h3n2_even$type <- "H3N2"

even <- rbind(h3n2_even)
rownames(even) <- NULL
even2 <- even %>% 
  filter(mean != "mean") %>%
  filter(str_detect(Summary.Statistic,"clock.epoch")) %>%
  filter(!str_detect(Summary.Statistic,"1|5")) %>%
  mutate(mean = as.numeric(mean)) %>%
  mutate(`95% HPD interval` = str_remove_all(`95% HPD interval`,"\\[|\\]")) %>%
  mutate(HPD_low = as.numeric(str_trim(sapply(str_split(`95% HPD interval`,","), function(x) x[1])))) %>%
  mutate(HPD_upp = as.numeric(str_trim(sapply(str_split(`95% HPD interval`,","), function(x) x[2]))))

ggplot(data = even2) +
  geom_errorbar(aes(x = Summary.Statistic, ymin = HPD_low, ymax = HPD_upp), width = 0.1, position = position_dodge(width = 0.7), size = 0.3)+
  geom_point(aes(x = Summary.Statistic, y = mean), shape = 21, size = 3, alpha = 1, position = position_dodge(width = 0.7), fill = colors[c(2)] )+
  scale_x_discrete(labels = c("During A/H1N1\npandemic season", 
                              "During inter-\npandemic period","During COVID-19\npandemic season"))+
  scale_y_continuous(limits = c(-0.1,4),expand = c(0,0))+
  theme_bw()+
  theme( 
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.background = element_rect(fill = "transparent",color = "transparent"),
    panel.border = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    plot.background =  element_rect(fill = "transparent",color = "transparent"),
    legend.background = element_rect(fill = "transparent",color = "transparent"))+
  labs(x = "", y = "Dispersal rate\n(events per lineage per year)",tag = "d")-> p4

#==output==
pdf("Fig5.pdf", height = 7, width = 11)
((p1 + p2 + p3) + plot_layout(nrow = 1,ncol = 3,heights = c(1,1,1), guides = "collect")&theme(legend.position = "bottom")) -> fig1
viewport(x = 0, y = 0.5, height = 0.5,  width = 0.75,just = c("left", "bottom")) -> vp1
viewport(x = 0.35, y = 0, height = 0.51, width = 0.65,just = c("left", "bottom")) -> vp2
viewport(x = 0.74, y = 0.495, height = 0.5, width = 0.26,just = c("left", "bottom")) -> vp4
viewport(x = 0, y = 0, height = 0.51, width = 0.35,just = c("left", "bottom")) -> vp5
print(fig1, vp = vp1)
print(p4, vp = vp4)
print(p1_4, vp = vp5)
print(p2_tmp, vp = vp2)
dev.off()