#==load packages==
library(ggsci)
library(ggplot2)
library(stringr)
library(lubridate)
library(data.table)
library(tidyverse)
library(grid)

#==define color==
color <- pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)
color1 <- pal_npg("nrc", alpha =1)(10)

#==read data==
strains2 <- readRDS("../data/antigenic_distance/h3n2_p1_antigenic_dis.rds")
strains2$strain_type1 <- strains2$strain_type
strains2$strain_type1[strains2$strain_type1 %in% c("Southeast_Asia_nonpersistent", "Southeast_Asia_persistent")] <- "Southeast_Asia"

#==plot (panel a)==
ggplot(data = strains2,aes(x = distance, y = Date1))+
  annotate("rect", xmin = 0.5, xmax = 6.05,ymin = seq(2010.5, 2019.5, 2),
           ymax = seq(2011.5, 2020.5, 2), alpha = 0.2, fill = "grey50")+
  annotate("rect", xmin = 0.5, xmax = 6.05,ymin = seq(2011.5, 2019.5, 2),
           ymax = seq(2012.5, 2020.5, 2), alpha = 0.2, fill = "grey75")+
  annotate("rect", xmin = 0.5, xmax = 6.05,ymin = c(2009.5, 2020.5),
           ymax = c(2010.5, 2021.5), alpha = c(0.3,0.2), fill = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6:7)])+
  geom_point(aes(color = strain_type1), alpha = 0.6)+
  geom_smooth(method = "lm", color = "grey", se = T)+
  theme_bw()+
  scale_color_manual("Geographic locations", values = c(color[c(2,1)]),
                     labels = c("Southeastern Asia","Temperate regions"))+
  scale_y_continuous(breaks = seq(2007,2024,1), limits = c(2006.5, 2024.5), expand = c(0,0))+
  scale_x_continuous(limits = c(0.5, 6.05), expand = c(0,0), breaks = seq(1,6,1))+
  theme(legend.position = c(0.8,0.08),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.grid.minor.y = element_blank())+
  labs(x = "Sequence-based antigenic distance\nfrom A/Wisconsin/67/2005", y = "Time",  tag = "a")-> p1

#==data cleaning==
fit <- lm(strains2$Date1 ~ strains2$distance)
fit$coefficients
strains2$predict_y <- predict(fit, data.frame(x = strains2$distance))
strains2$predict_x <- (strains2$Date1 - 2004.102294)/3.507084

strains2$distance_lag <- strains2$distance - strains2$predict_x
strains2$time_lag <- strains2$predict_y - strains2$Date1

strains3 <- strains2 %>% filter(Date >= as.Date("2007-07-01")) %>% filter(Date <= as.Date("2023-06-30"))       
cut(strains3$Date, breaks = c(as.Date("2006-12-31"), as.Date("2009-06-30"),
                              as.Date("2010-06-30"), as.Date("2020-06-30"),
                              as.Date("2021-06-30"), as.Date("2023-07-01")),
    labels = c("Before A/H1N1 pandemic", "A/H1N1 pandemic season",
               "Interpandemic period", "COVID-19 pandemic season", "After COVID-19 pandemic"), right = T) -> strains3$season
cut(strains3$Date, breaks = c(as.Date("2007-06-30"), as.Date("2008-06-30"),
                              as.Date("2009-06-30"), as.Date("2010-06-30"),
                              as.Date("2011-06-30"), as.Date("2012-06-30"),
                              as.Date("2013-06-30"), as.Date("2014-06-30"),
                              as.Date("2015-06-30"), as.Date("2016-06-30"),
                              as.Date("2017-06-30"), as.Date("2018-06-30"),
                              as.Date("2019-06-30"), as.Date("2020-06-30"),
                              as.Date("2021-06-30"), as.Date("2022-06-30"),as.Date("2023-07-01")),
    labels = c("2007/2008", "2008/2009","2009/2010", "2010/2011",
               "2011/2012", "2012/2013","2013/2014", "2014/2015",
               "2015/2016", "2016/2017","2017/2018", "2018/2019",
               "2019/2020", "2020/2021","2021/2022", "2022/2023"), right = T) -> strains3$season1

strains3$strain_type1 <- as.character(strains3$strain_type)
strains3$strain_type2 <- as.character(strains3$Location)
strains3$strain_type2[str_detect(strains3$strain_type1,"Southeast_Asia")] <- 
  strains3$strain_type1[str_detect(strains3$strain_type1,"Southeast_Asia")]
strains3$strain_type2[strains3$strain_type2 == "Southeast_Asia_nonpersistent"] <- "Southeastern Asia (NPL)"
strains3$strain_type2[strains3$strain_type2 == "Southeast_Asia_persistent"] <- "Southeastern Asia (PL)"
strains3$strain_type2[strains3$strain_type2 == "NorthAmerica"] <- "North America"
strains3$strain_type2[strains3$strain_type2 == "SouthAmerica"] <- "South America"
strains3$strain_type2[strains3$strain_type2 == "Africa"] <- "Southern Africa"
strains3$strain_type2 <- paste0("      ", strains3$strain_type2)

strains3$strain_type1[str_detect(strains3$strain_type1,"Southeast_Asia")] <- "Southeastern Asia"
strains3$strain_type1[strains3$strain_type1 == "Temperate_region"] <- "Temperate regions"

#==plot (panel b)==
stat1 <- strains3 %>% filter(Date >= as.Date("2009-07-01") & Date <= as.Date("2021-06-30")) %>% 
  group_by(strain_type1) %>% summarise(n = n(), mean_dis = mean(distance_lag), 
                                       sd_dis = sd(distance_lag),
                                       mean_time = mean(time_lag), 
                                       sd_time = sd(time_lag)) %>%
  rename("strain_type" = "strain_type1")
stat2 <- strains3 %>% filter(Date >= as.Date("2009-07-01") & Date <= as.Date("2021-06-30")) %>% 
  group_by(strain_type2) %>% summarise(n = n(), mean_dis = mean(distance_lag), 
                                       sd_dis = sd(distance_lag),
                                       mean_time = mean(time_lag), 
                                       sd_time = sd(time_lag)) %>%
  rename("strain_type" = "strain_type2")
stat3 <- strains3 %>% filter(Date >= as.Date("2009-07-01") & Date <= as.Date("2021-06-30")) %>% 
  group_by(strain_type1, season) %>% summarise(n = n(), mean_dis = mean(distance_lag),
                                               sd_dis = sd(distance_lag),
                                               mean_time = mean(time_lag),
                                               sd_time = sd(time_lag)) %>%
  rename("strain_type" = "strain_type1")

stat4 <- strains3 %>% filter(Date >= as.Date("2009-07-01") & Date <= as.Date("2021-06-30")) %>% 
  group_by(strain_type2, season) %>% summarise(n = n(), mean_dis = mean(distance_lag), 
                                               sd_dis = sd(distance_lag),
                                               mean_time = mean(time_lag),
                                               sd_time = sd(time_lag)) %>%
  rename("strain_type" = "strain_type2")

stat5 <- strains3 %>% group_by(strain_type1, season1) %>% summarise(mean_dis = mean(distance_lag), 
                                                                    sd_dis = sd(distance_lag),
                                                                    mean_time = mean(time_lag),
                                                                    sd_time = sd(time_lag)) %>%
  rename("strain_type" = "strain_type1")

stat6 <- strains3 %>% group_by(strain_type2, season1) %>% summarise(mean_dis = mean(distance_lag), 
                                                                    sd_dis = sd(distance_lag),
                                                                    mean_time = mean(time_lag),
                                                                    sd_time = sd(time_lag)) %>%
  rename("strain_type" = "strain_type2")

stat_plot1 <- rbind(stat1, stat2)
stat_plot2 <- rbind(stat3, stat4)
stat_plot3 <- rbind(stat5, stat6)

factor(stat_plot1$strain_type, levels = rev(c("Southeastern Asia", "      Southeastern Asia (PL)","      Southeastern Asia (NPL)",
                                              "Temperate regions",  "      Oceania",  "      North America","      Europe", "      South America","      Southern Africa"))) -> stat_plot1$strain_type
factor(stat_plot2$strain_type, levels = rev(c("Southeastern Asia", "      Southeastern Asia (PL)","      Southeastern Asia (NPL)",
                                              "Temperate regions",  "      Oceania",  "      North America","      Europe", "      South America","      Southern Africa"))) -> stat_plot2$strain_type
factor(stat_plot3$strain_type, levels = rev(c("Southeastern Asia", "      Southeastern Asia (PL)","      Southeastern Asia (NPL)",
                                              "Temperate regions",  "      Oceania",  "      North America","      Europe", "      South America","      Southern Africa"))) -> stat_plot3$strain_type
ggplot() +
  geom_vline(xintercept = 0, linetype = 2, color = "red")+
  geom_point(data = stat_plot2, aes(x = mean_dis, y = strain_type, color = season, size = n))+
  geom_point(data = stat_plot1, aes(x = mean_dis, y = strain_type,fill = "grey"), color = "black",  shape =23, alpha = 0.5, size = 5)+
  scale_x_continuous(limits = c(-0.4, 0.4), expand = c(0,0), breaks = c(-0.4,-0.3,-0.2,-0.1,0, 0.1,0.2,0.3,0.4),
                     labels = c(-0.4,-0.3,-0.2,-0.1,0, 0.1,0.2,0.3,0.4),sec.axis = sec_axis(trans=~.* 1,
                                                                                            name= "Months",
                                                                                            breaks = c(-1/3.507084, -0.5/3.507084,0, 0.5/3.507084, 1/3.507084),
                                                                                            labels = c(-12, -6, 0, 6, 12)))+
  theme_bw()+
  theme(axis.text.y = element_text(hjust = 0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.spacing.y = unit(-0.01, 'cm')
  )+
  guides(color = guide_legend(nrow = 3, title.position = "top", order = 1, direction = "vertical"),
         size = guide_legend(ncol = 2, title.position = "top"),
         fill = guide_legend(order = 2, direction = "vertical", title = element_blank() ))+
  scale_color_manual("Periods", values = pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)[c(6,4,7)])+
  # scale_fill_manual("", values = "grey", labels = "A")+
  scale_fill_manual("", values = "grey", labels = "Average distance\nmonths across\nall periods")+
  scale_size_continuous("Number of sequences", limit = c(1,1500), range = c(1.5,5), breaks = c(10,50,100,500,1000,1500))+
  labs(x = "Mean sequence-based antigenic distance\nto the fitted line", y = "",tag = "b")-> p2

#==output==
pdf("Fig4.pdf", height = 6, width = 12)
viewport(x = 0, y = 0, height = 1,  width = 0.5,just = c("left", "bottom")) -> vp1
viewport(x = 0.495, y = 0, height = 1, width = 0.5,just = c("left", "bottom")) -> vp2
print(p2, vp = vp2)
print(p1, vp = vp1)
dev.off()