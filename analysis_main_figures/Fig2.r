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
library(coda)
library(zoo)
library(ggsci)

#==define function==
summarise_hpd_lower <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[1])
}

summarise_hpd_upper <- function(x) {
  if(length(x) <= 1) {
    return(x[1]);
  }
  return(HPDinterval(as.mcmc(x),prob = 0.95)[2])
}

#==define week==
date_match <- data.frame(date_plot = seq(as.Date("1995-01-05"),as.Date("2024-01-18"),7))
date_match$ISO_YEAR <- isoyear(date_match$date_plot)
date_match$ISO_WEEK <- isoweek(date_match$date_plot)

date_match1 <- data.frame(date_plot = seq(as.Date("1995-01-05"),as.Date("2024-01-18"),7))
date_match1$ISO_year <- isoyear(date_match1$date_plot)
date_match1$ISO_week <- isoweek(date_match1$date_plot)

#==read h3n2 data==
mj_h3n2_p1 <- rbind(read.delim("../phylogeographic_analysis/jump_two_state/h3n2_p1_Jumps_0_1M.txt", sep = ","),
                    read.delim("../phylogeographic_analysis/jump_two_state/h3n2_p1_Jumps_1_2M.txt", sep = ","))
mj_h3n2_p1$date <- decimal2Date(decimal_date(as.Date("2023-12-13")) - mj_h3n2_p1$time)

#==data cleaning==
jump1 <- mj_h3n2_p1 %>%
  mutate(tran_type = paste0(startLocation,"_", endLocation)) %>%
  mutate(ISO_YEAR = isoyear(date), ISO_WEEK = isoweek(date)) %>%
  left_join(date_match) %>%
  group_by(date_plot, treeId, tran_type) %>%
  summarise(jump_no = n())

length(unique(jump1$treeId))
length(unique(jump1$tran_type))
length(unique(jump1$date_plot))
length(seq(as.Date("2006-12-28"),as.Date("2023-12-07"),7))
jump1_tmp <- rbind(data.frame(treeId = rep(unique(jump1$treeId), each = 885),
                              date_plot = rep(seq(as.Date("2006-12-28"),as.Date("2023-12-07"),7), 901)) %>% mutate(tran_type = "Southeast_Asia_Temperate_region"),
                   data.frame(treeId = rep(unique(jump1$treeId), each = 885),
                              date_plot = rep(seq(as.Date("2006-12-28"),as.Date("2023-12-07"),7), 901)) %>% mutate(tran_type = "Temperate_region_Southeast_Asia"))
jump1_tmp <- left_join(jump1_tmp, jump1)
nrow(jump1_tmp[!is.na(jump1_tmp$jump_no),])
jump1_tmp$jump_no[is.na(jump1_tmp$jump_no)] <- 0

#1.Movement per week
jump2 <- jump1_tmp %>%
  group_by(date_plot, tran_type) %>%
  summarise(mean = mean(jump_no),
            sd = sd(jump_no),
            upp = summarise_hpd_upper(jump_no),
            low = summarise_hpd_lower(jump_no)) 
range(jump2$date_plot)

jump2_tmp <- data.frame(date_plot = rep(seq(as.Date("2006-12-28"),as.Date("2023-12-07"),7), 2),
                        tran_type = rep(c("Southeast_Asia_Temperate_region", "Temperate_region_Southeast_Asia"), each = 885)) %>%
  left_join(jump2)
jump2_tmp$mean[is.na(jump2_tmp$mean)] <- 0

#2.Overall movement
jump2_overall <- jump1_tmp %>%
  filter(date_plot <= as.Date("2021-06-30") & date_plot >= as.Date("2009-07-01")) %>%
  group_by(treeId, tran_type) %>%
  summarise(overall_no = sum(jump_no)) %>%
  ungroup() %>%
  group_by(tran_type) %>%
  summarise(mean = mean(overall_no),
            upp = summarise_hpd_upper(overall_no),
            low = summarise_hpd_lower(overall_no)) 

#==plot==
jump2_tmp <- jump2_tmp[order(jump2_tmp$tran_type, jump2_tmp$date_plot),]
jump3 <- jump2_tmp[,1:6] %>% group_by(tran_type) %>%
  mutate(mean_roll = rollmean(mean, k=5, fill=NA, align='center'))

#==Supplementary figure
ggplot(data = jump3[jump3$date_plot <= as.Date("2023-06-30") & jump3$date_plot >= as.Date("2007-07-01"),]) +
  annotate("rect", xmin = c(as.Date("2007-07-01"), as.Date("2010-07-01"), as.Date("2021-07-01")), 
           xmax = c(as.Date("2009-06-30"),as.Date("2020-07-01"), as.Date("2023-07-01")),
           ymin = 0, ymax = 4, alpha = 0.2, fill = c("white","grey50","white"))+
  annotate("rect", xmin = as.Date("2009-07-01"), xmax = as.Date("2010-06-30"),
           ymin = 0, ymax = 4, alpha = 0.3, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-06-30"),
           ymin = 0, ymax = 4, alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  geom_line(aes(x = date_plot, y = mean_roll, group = tran_type),color = "black")+
  geom_ribbon(aes(x = date_plot, ymax = mean_roll, ymin = 0, fill = tran_type, group = tran_type), alpha = 0.7)+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y", expand = c(0.01,0))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 4), expand = c(0,0), labels = seq(0,1,0.25))+
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    # axis.title.x = element_blank(),
    legend.position = "none")+
  scale_fill_manual("Markov Jumps direction", values = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(1,8)],
                    labels = c("From Southeastern Asia to temperate regions",
                               "From temperate regions to Southeastern Asia"))+
  annotate("text", x = c(as.Date("2008-07-01"), as.Date("2015-07-01"),as.Date("2022-07-01")), y = 2.2,size =3,
           label = c("Before\nA/H1N1\npandemic\nseason (P1)","Interpandemic\nseasons\n(P3)","After\nCOVID-19\npandemic\nseason (P5)"))+
  annotate("text", x = c( as.Date("2010-01-01"),as.Date("2021-01-01")), y = 3.2,size =3,
           label = c("A/H1N1\npandemic\nseason (P2)","COVID-19\npandemic\nseason (P4)"))+
  guides(fill = guide_legend(title.position = "top",ncol = 1))+
  labs(x = "Date", y = "Relative viral movement intensity", subtitle = "A/H3N2", tag = "a")-> p1

#==panel a
ggplot(jump2_overall) +
  geom_bar(aes(x = tran_type, y = mean), stat = "identity", 
           fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(1,8)], width = 0.5, color = "black")+
  geom_errorbar(aes(x = tran_type, ymin = low, ymax = upp), width = 0.2)+
  theme_bw()+
  scale_x_discrete("Markov Jumps direction")+
  scale_y_continuous("Relative viral movement intensity", limits = c(0,600), breaks = seq(0,600,120), labels = seq(0,1,0.2))+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank())+
  labs(subtitle = "A/H3N2", tag = "a")-> p1_m

#==panel c
jump4 <- left_join(jump3, date_match) %>% filter(tran_type == "Southeast_Asia_Temperate_region")

tmp1 <- jump4 %>% filter((ISO_YEAR == 2011 & ISO_WEEK >= 27) | (ISO_YEAR == 2012 & ISO_WEEK <= 26))
peak1 <- tmp1$ISO_WEEK[which.max(tmp1$mean_roll)]
tmp2 <- jump4 %>% filter((ISO_YEAR == 2012 & ISO_WEEK >= 27) | (ISO_YEAR == 2013 & ISO_WEEK <= 26))
peak2 <- tmp2$ISO_WEEK[which.max(tmp2$mean_roll)]
tmp3 <- jump4 %>% filter((ISO_YEAR == 2013 & ISO_WEEK >= 27) | (ISO_YEAR == 2014 & ISO_WEEK <= 26))
peak3 <- tmp3$ISO_WEEK[which.max(tmp3$mean_roll)]
tmp4 <- jump4 %>% filter((ISO_YEAR == 2014 & ISO_WEEK >= 27) | (ISO_YEAR == 2015 & ISO_WEEK <= 26))
peak4 <- tmp4$ISO_WEEK[which.max(tmp4$mean_roll)]
tmp5 <- jump4 %>% filter((ISO_YEAR == 2015 & ISO_WEEK >= 27) | (ISO_YEAR == 2016 & ISO_WEEK <= 26))
peak5 <- tmp5$ISO_WEEK[which.max(tmp5$mean_roll)]
tmp6 <- jump4 %>% filter((ISO_YEAR == 2016 & ISO_WEEK >= 27) | (ISO_YEAR == 2017 & ISO_WEEK <= 26))
peak6 <- tmp6$ISO_WEEK[which.max(tmp6$mean_roll)]
tmp7 <- jump4 %>% filter((ISO_YEAR == 2017 & ISO_WEEK >= 27) | (ISO_YEAR == 2018 & ISO_WEEK <= 26))
peak7 <- tmp7$ISO_WEEK[which.max(tmp7$mean_roll)]
tmp8 <- jump4 %>% filter((ISO_YEAR == 2018 & ISO_WEEK >= 27) | (ISO_YEAR == 2019 & ISO_WEEK <= 26))
peak8 <- tmp8$ISO_WEEK[which.max(tmp8$mean_roll)]
tmp9 <- jump4 %>% filter((ISO_YEAR == 2019 & ISO_WEEK >= 27) | (ISO_YEAR == 2020 & ISO_WEEK <= 26))
peak9 <- tmp9$ISO_WEEK[which.max(tmp9$mean_roll)]
sort(c(peak1, peak2, peak4, peak6, peak7, peak8, peak9)) 
#peak: 49 week (from week 27 of one year to week 26 of next year) (outlier weeks: week 3 & 5)
peak_set <- 49

#align with the median week
tmp1 <- tmp1 %>% 
  mutate(date_plot1 = date_plot + (peak_set - peak1)*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2011 & ISO_week >= 27) | (ISO_year == 2012 & ISO_week <= 26)) %>%
  mutate(season = "2011/2012")

tmp2 <- tmp2 %>% 
  mutate(date_plot1 = date_plot + (peak_set - peak2)*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2012 & ISO_week >= 27) | (ISO_year == 2013 & ISO_week <= 26))%>%
  mutate(season = "2012/2013")

tmp4 <- tmp4 %>% 
  mutate(date_plot1 = date_plot + (peak_set - peak4)*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2014 & ISO_week >= 27) | (ISO_year == 2015 & ISO_week <= 26)) %>%
  mutate(season = "2014/2015")

tmp6 <- tmp6 %>% 
  mutate(date_plot1 = date_plot + (peak_set - peak6)*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2016 & ISO_week >= 27) | (ISO_year == 2017 & ISO_week <= 26)) %>%
  mutate(season = "2016/2017")

tmp7 <- tmp7 %>% 
  mutate(date_plot1 = date_plot + (peak_set - peak7)*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2017 & ISO_week >= 27) | (ISO_year == 2018 & ISO_week <= 26)) %>%
  mutate(season = "2017/2018")

tmp8 <- tmp8 %>% 
  mutate(date_plot1 = date_plot + (peak_set - (peak8 + nrow(tmp8) ))*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2018 & ISO_week >= 27) | (ISO_year == 2019 & ISO_week <= 26)) %>%
  mutate(season = "2018/2019")

tmp9 <- tmp9 %>% 
  mutate(date_plot1 = date_plot + (peak_set - peak9)*7) %>% 
  left_join(date_match1, by = c("date_plot1" = "date_plot")) %>%
  filter((ISO_year == 2019 & ISO_week >= 27) | (ISO_year == 2020 & ISO_week <= 26)) %>%
  mutate(season = "2019/2020")

#combine data
jump5 <- rbind(tmp1, tmp2,  tmp4, tmp6, tmp7, tmp8,tmp9)
factor(jump5$ISO_week, levels = c(27:52,1:26)) -> jump5$ISO_week
jump5$ISO_week1 <- as.numeric(jump5$ISO_week)
jump6 <- jump5 %>% 
  group_by(ISO_week) %>% 
  summarise(average_count = mean(mean_roll)) %>% 
  mutate(type = "Average Markov Jumps after aligning at the\nmedian peak for the interpandemic seasons")
jump6 <- jump6 %>% 
  mutate(average_count_roll = rollmean(average_count, k=5, fill=NA, align='center'))
tmp_pandemic <- as.data.frame(jump4) %>% 
  filter(((ISO_YEAR == 2009 & ISO_WEEK >= 27) | (ISO_YEAR == 2010 & ISO_WEEK <= 26))|
           ((ISO_YEAR == 2020 & ISO_WEEK >= 27) | (ISO_YEAR == 2021 & ISO_WEEK <= 26))) %>%
  mutate(type = ifelse(ISO_YEAR <= 2011, "During the H1N1 pandemic", "During the COVID-19 pandemic")) %>%
  select(c("ISO_WEEK","mean_roll", "type")) %>%
  rename(c("ISO_week" = "ISO_WEEK", "average_count" = "mean_roll", "type"= "type")) %>%
  filter(!ISO_week %in% c(53))
# filter(!ISO_week %in% c(25:28))
jump7 <- rbind(jump6[,c(1,2,3)],tmp_pandemic) %>%
  mutate(x_text = as.character(ISO_week)) %>%
  mutate(x_text = ifelse(x_text %in% c("10","20","30","40","50"),x_text, "")) %>%
  mutate(ISO_week1 = as.numeric(ISO_week))

#plot
factor(jump7$type, levels = unique(jump7$type)[c(2,1,3)]) -> jump7$type
ggplot(jump7) +
  geom_line(data = jump5, aes(x = ISO_week, y =  mean_roll, group = season), color = "grey")+
  geom_line(aes(x = ISO_week, y = average_count,color = type, group = type), linewidth = 1.5)+
  theme_bw()+
  geom_vline(xintercept = 23, linetype = 2)+
  geom_segment(aes(x= 15 , y= 2.5 , xend= 23 , yend= 2.5 ),
               arrow = arrow(length=unit(0.2, 'cm')),lwd= 0.2)+
  annotate("text", x = 14, y = 2.5,size =4, label = "Median\nweek", hjust = 1)+
  scale_x_discrete(breaks = unique(jump7$ISO_week[jump7$x_text != ""]),
                   labels = unique(jump7$x_text[jump7$x_text != ""]))+
  theme(legend.position = c(0.75, 0.85),
        panel.grid = element_blank(),
        # panel.grid.minor.x = element_blank()
  )+
  scale_y_continuous(limits = c(-0.2,3), expand = c(0,0), breaks = seq(0,3,3/5), labels = seq(0,1,0.2))+
  scale_color_manual("Periods",values = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6,4,7)],
                     labels = c("During the A/H1N1 pandemic season", 
                                "During the interpandemic seasons (average)", 
                                "During the COVID-19 pandemic season"))+
  guides(color = guide_legend(title.position = "top",ncol = 1, nrow = 3))+
  labs(x = "ISO Week\n", y = "Relative viral movement\n(Southeastern Asia to temperate regions)", tag = "c", subtitle = "A/H3N2")-> p1_r
