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
library(treeio)
library(zoo)
library(patchwork)
library(ggsci)
library(scales)

#==define week==
date_match <- data.frame(date = seq(as.Date("1995-01-02"),as.Date("2024-01-15"),7))
date_match$ISO_YEAR <- isoyear(date_match$date)
date_match$ISO_WEEK <- isoweek(date_match$date)

date_match1 <- data.frame(date1 = seq(as.Date("1995-01-02"),as.Date("2024-01-15"),7))
date_match1$ISO_year <- isoyear(date_match1$date1)
date_match1$ISO_week <- isoweek(date_match1$date1)

#==read epi data==
epi_final <- read.csv("../data/epi_data/epi_data_0430.csv") %>%
  mutate(date = as.Date(date)) 

#==data cleaning & plot (panel a)==
epi_final1 <- epi_final[,1:3]
epi_final1 <- epi_final1[order(epi_final1$region_final1, epi_final1$date),]
epi_final1$diff <- c(NA,diff(epi_final1$date))
epi_final2 <- epi_final1[,1:4] %>%
  group_by(region_final1) %>%
  mutate(specimen_roll = rollmean(specimen, k=5, fill=NA, align='center'))

# pal_npg("nrc", alpha =1)(10)
# show_col(pal_lancet(palette = c("lanonc"), alpha = 0.9)(9))
# show_col(pal_lancet(palette = c("lanonc"), alpha = 0.2)(9))

ggplot(data = epi_final2[epi_final2$date >= as.Date("2007-01-01") & epi_final2$date <= as.Date("2023-12-25"),]) +
  annotate("rect", xmin = as.Date("2009-07-01"), xmax = as.Date("2010-06-30"),
           ymin = 90, ymax = 400000, alpha = 0.3, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6)])+
  annotate("rect", xmin = as.Date("2020-07-01"), xmax = as.Date("2021-06-30"),
           ymin = 90, ymax = 400000, alpha = 0.2, fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7)])+
  annotate("rect", xmin = c(as.Date("2010-07-01"),as.Date("2012-07-01"),as.Date("2014-07-01"),as.Date("2016-07-01"),as.Date("2018-07-01")),
           xmax = c(as.Date("2011-06-30"),as.Date("2013-06-30"),as.Date("2015-06-30"),as.Date("2017-06-30"),as.Date("2019-06-30")),
           ymin = 90, ymax = 400000, alpha = 0.2, fill = "grey50")+
  annotate("rect", xmin = c(as.Date("2011-07-01"),as.Date("2013-07-01"),as.Date("2015-07-01"),as.Date("2017-07-01"),as.Date("2019-07-01")),
           xmax = c(as.Date("2012-06-30"),as.Date("2014-06-30"),as.Date("2016-06-30"),as.Date("2018-06-30"),as.Date("2020-06-30")),
           ymin = 90, ymax = 400000, alpha = 0.2, fill = "grey75")+
  geom_vline(xintercept = c(as.Date("2010-07-01"),as.Date("2020-07-01")), linetype =2, color = "grey30") +
  geom_segment(aes(x= as.Date("2010-07-01") , y= 500 , xend= as.Date("2013-06-30") , yend= 500 ), 
               arrow = arrow(length=unit(0.2, 'cm')),lwd= 0.2)+
  geom_segment(aes(x= as.Date("2020-07-01") , y= 500 , xend= as.Date("2017-06-30") , yend= 500 ), 
               arrow = arrow(length=unit(0.2, 'cm')),lwd= 0.2)+
  annotate("text", x = as.Date("2015-07-01"), y = 500,size =3, label = "Interpandemic seasons")+
  annotate("text", x = as.Date("2009-10-18"), y = 290,size =2.5, label = "A/H1N1\npandemic\nseason")+
  annotate("text", x = as.Date("2021-03-15"), y = 290,size =2.5, label = "COVID-19\npandemic\nseason")+
  geom_line(aes(x = date+3, y = specimen_roll, color = region_final1))+
  theme_bw()+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y", expand = c(0.01, 0.01))+
  theme(legend.title = element_blank(),
        plot.margin = margin(0,0.2,0,0, "cm"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  scale_y_log10(limits = c(90, 400000), expand = c(0,0),breaks = c(100, 1000,10000,100000), labels =c(expression("10"^2),expression("10"^3),expression("10"^4),expression("10"^5)) )+
  scale_color_manual(values = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(2,1,5)])+
  labs(x = "Date", y = "Number of specimens\nfor influenza testing")-> p1

#==data cleaning by regions==
epi_final3 <- epi_final[order(epi_final$region_final1, epi_final$date),]
epi_final3$diff <- c(NA,diff(epi_final3$date))
epi_final4 <- epi_final3[,1:10] %>%
  group_by(region_final1) %>%
  mutate(specimen_roll = rollmean(specimen, k=5, fill=NA, align='center'),
         h3n2_roll = rollmean(h3n2, k=5, fill=NA, align='center'),
         h3n2_LL_roll = rollmean(h3n2_LL, k=5, fill=NA, align='center'),
         h3n2_UL_roll = rollmean(h3n2_UL, k=5, fill=NA, align='center'),
         bv_roll = rollmean(bv, k=5, fill=NA, align='center'),
         bv_LL_roll = rollmean(bv_LL, k=5, fill=NA, align='center'),
         bv_UL_roll = rollmean(bv_UL, k=5, fill=NA, align='center'))
epi_final4 <- left_join(epi_final4, date_match)

#====NH (H3N2)====[take panel c as an example]
factor(epi_final4$ISO_WEEK, levels = c(27:53,1:26)) -> epi_final4$ISO_WEEK1

tmp0 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2010 & ISO_WEEK >= 27) | (ISO_YEAR == 2011 & ISO_WEEK <= 26))
peak0 <- tmp0$ISO_WEEK[which.max(tmp0$h3n2_roll*100/tmp0$specimen_roll)]
tmp1 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2011 & ISO_WEEK >= 27) | (ISO_YEAR == 2012 & ISO_WEEK <= 26))
peak1 <- tmp1$ISO_WEEK[which.max(tmp1$h3n2_roll*100/tmp1$specimen_roll)]
tmp2 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2012 & ISO_WEEK >= 27) | (ISO_YEAR == 2013 & ISO_WEEK <= 26))
peak2 <- tmp2$ISO_WEEK[which.max(tmp2$h3n2_roll*100/tmp2$specimen_roll)]
tmp3 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2013 & ISO_WEEK >= 27) | (ISO_YEAR == 2014 & ISO_WEEK <= 26))
peak3 <- tmp3$ISO_WEEK[which.max(tmp3$h3n2_roll*100/tmp3$specimen_roll)]
tmp4 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2014 & ISO_WEEK >= 27) | (ISO_YEAR == 2015 & ISO_WEEK <= 26))
peak4 <- tmp4$ISO_WEEK[which.max(tmp4$h3n2_roll*100/tmp4$specimen_roll)]
tmp5 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2015 & ISO_WEEK >= 27) | (ISO_YEAR == 2016 & ISO_WEEK <= 26))
peak5 <- tmp5$ISO_WEEK[which.max(tmp5$h3n2_roll*100/tmp5$specimen_roll)]
tmp6 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2016 & ISO_WEEK >= 27) | (ISO_YEAR == 2017 & ISO_WEEK <= 26))
peak6 <- tmp6$ISO_WEEK[which.max(tmp6$h3n2_roll*100/tmp6$specimen_roll)]
tmp7 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2017 & ISO_WEEK >= 27) | (ISO_YEAR == 2018 & ISO_WEEK <= 26))
peak7 <- tmp7$ISO_WEEK[which.max(tmp7$h3n2_roll*100/tmp7$specimen_roll)]
tmp8 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2018 & ISO_WEEK >= 27) | (ISO_YEAR == 2019 & ISO_WEEK <= 26))
peak8 <- tmp8$ISO_WEEK[which.max(tmp8$h3n2_roll*100/tmp8$specimen_roll)]
tmp9 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2019 & ISO_WEEK >= 27) | (ISO_YEAR == 2020 & ISO_WEEK <= 26))
peak9 <- tmp9$ISO_WEEK[which.max(tmp9$h3n2_roll*100/tmp9$specimen_roll)]

pd1 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2009 & ISO_WEEK >= 27) | (ISO_YEAR == 2010 & ISO_WEEK <= 26))
pd2 <- epi_final4 %>% filter(region_final1 == "Temperate regions (NH)") %>%
  filter((ISO_YEAR == 2020 & ISO_WEEK >= 27) | (ISO_YEAR == 2021 & ISO_WEEK <= 26))

#choose median of peak week (from 27 week of one year to 26 week of next year) 
sort(c(peak0, peak1, peak2, peak3, peak4,peak5, peak6, peak7, peak8, peak9)) #peak: week 1 (or week 2)
peak_set <- 1

#align with the median week
tmp0 <- tmp0 %>% 
  mutate(date1 = date + (nrow(tmp0) + peak_set - peak0)*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2010 & ISO_week >= 27) | (ISO_year == 2011 & ISO_week <= 26)) %>%
  mutate(season = "2010/2011")

tmp1 <- tmp1 %>% 
  mutate(date1 = date + (peak_set - (peak1))*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2011 & ISO_week >= 27) | (ISO_year == 2012 & ISO_week <= 26)) %>%
  mutate(season = "2011/2012")

tmp2 <- tmp2 %>% 
  mutate(date1 = date + (nrow(tmp2) +peak_set - peak2)*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2012 & ISO_week >= 27) | (ISO_year == 2013 & ISO_week <= 26))%>%
  mutate(season = "2012/2013")

tmp3 <- tmp3 %>%
  mutate(date1 = date + (peak_set - (peak3 ))*7) %>%
  left_join(date_match1) %>%
  filter((ISO_year == 2013 & ISO_week >= 27) | (ISO_year == 2014 & ISO_week <= 26))%>%
  mutate(season = "2013/2014")

tmp4 <- tmp4 %>% 
  mutate(date1 = date + (nrow(tmp4) +peak_set - peak4)*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2014 & ISO_week >= 27) | (ISO_year == 2015 & ISO_week <= 26)) %>%
  mutate(season = "2014/2015")

tmp5 <- tmp5 %>% 
  mutate(date1 = date + (peak_set - (peak5))*7)  %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2015 & ISO_week >= 27) | (ISO_year == 2016 & ISO_week <= 26)) %>%
  mutate(season = "2015/2016")
tmp5 <- tmp5 %>% filter(ISO_week != 53)

tmp6 <- tmp6 %>% 
  mutate(date1 = date + (peak_set - peak6)*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2016 & ISO_week >= 27) | (ISO_year == 2017 & ISO_week <= 26)) %>%
  mutate(season = "2016/2017")

tmp7 <- tmp7 %>% 
  mutate(date1 = date + (peak_set - peak7)*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2017 & ISO_week >= 27) | (ISO_year == 2018 & ISO_week <= 26)) %>%
  mutate(season = "2017/2018")

tmp8 <- tmp8 %>% 
  mutate(date1 = date + (peak_set - (peak8))*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2018 & ISO_week >= 27) | (ISO_year == 2019 & ISO_week <= 26)) %>%
  mutate(season = "2018/2019")

tmp9 <- tmp9 %>% 
  mutate(date1 = date + (nrow(tmp9) +peak_set - peak9)*7) %>% 
  left_join(date_match1) %>%
  filter((ISO_year == 2019 & ISO_week >= 27) | (ISO_year == 2020 & ISO_week <= 26)) %>%
  mutate(season = "2019/2020")

#combine data
epi_average <- rbind(tmp0, tmp1, tmp2,  tmp3, tmp4, tmp5, tmp6, tmp7, tmp8,tmp9) %>% mutate(rate = h3n2_roll*100/specimen_roll)
factor(epi_average$ISO_week, levels = c(27:52,1:26)) -> epi_average$ISO_week
epi_average$ISO_week1 <- as.numeric(epi_average$ISO_week)

epi_average1 <- epi_average %>% group_by(ISO_week) %>% 
  summarise(average_rate = mean(rate), sd_rate = sd(rate)) %>% 
  mutate(type = "Mean of the positivity rate after aligning at the median\npeak for the interpandemic seasons")

tmp_pandemic <- as.data.frame(epi_final4) %>%
  filter(region_final1 == "Temperate regions (NH)") %>%
  filter(((ISO_YEAR == 2009 & ISO_WEEK >= 27) | (ISO_YEAR == 2010 & ISO_WEEK <= 26))|
           ((ISO_YEAR == 2020 & ISO_WEEK >= 27) | (ISO_YEAR == 2021 & ISO_WEEK <= 26))) %>%
  mutate(type = ifelse(ISO_YEAR <= 2011, "H1N1 pandemic", "COVID-19 pandemic")) %>%
  mutate(rate = h3n2_roll*100/specimen_roll)%>%
  select(c("ISO_WEEK","rate", "type")) %>%
  rename(c("ISO_week" = "ISO_WEEK","average_rate" = "rate", "type"= "type")) %>%
  mutate(sd_rate = NA) %>%
  filter(!ISO_week %in% c(53))
# filter(!ISO_week %in% c(25:28))

epi_average3 <- rbind(epi_average1, tmp_pandemic) %>%
  mutate(x_text = as.character(ISO_week)) %>%
  mutate(x_text = ifelse(x_text %in% c("10","20","30","40","50"), x_text, "")) %>%
  mutate(ISO_week1 = as.numeric(ISO_week))

ggplot(epi_average3) +
  geom_line(data = epi_average, aes(x = ISO_week, y = rate, group = season), color = "grey")+
  geom_line(aes(x = ISO_week, y = average_rate, color = type, group = type), linewidth = 1)+
  theme_bw()+
  geom_vline(xintercept = 27, linetype = 2)+
  # scale_x_continuous(breaks = c(7,17,27,37,47),
  #                    labels = c(-20,-10,0,10,20),
  #                    sec.axis = sec_axis(trans=~.* 1,
  #                                        name= "Week (pandemic seasons)",
  #                                        breaks = unique(epi_average3$ISO_week1[epi_average3$x_text != ""]),
  #                                        labels = unique(epi_average3$x_text[epi_average3$x_text != ""]) ))+
  scale_x_discrete(breaks = unique(epi_average3$ISO_week[epi_average3$x_text != ""]),
                   labels = unique(epi_average3$x_text[epi_average3$x_text != ""]))+
  geom_segment(aes(x= "45" , y= 22 , xend= "1" , yend= 22 ), 
               arrow = arrow(length=unit(0.2, 'cm')),lwd= 0.2)+
  annotate("text", x = "44", y = 22,size =4, label = "Median\nweek", hjust = 1)+
  theme(legend.position = "none",
        panel.grid = element_blank())+
  scale_y_continuous(limits = c(-1,26), expand = c(0,0))+
  scale_color_manual("Periods",values = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(7,6,4)])+
  guides(color = guide_legend(title.position = "top",ncol = 1, nrow = 3))+
  labs(x = "ISO Week", y = "Positivity rate of A/H3N2 (%)",  subtitle = "Temperate regions (NH)")-> p2_1
