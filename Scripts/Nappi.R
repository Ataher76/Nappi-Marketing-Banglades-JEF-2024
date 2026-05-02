library(readxl)
library(ggplot2)
library(openxlsx)
library(mvtnorm)
library(survival)
library(MASS)
library(TH.data)
library(tidyr)
library(tidyverse)
library(dplyr)
library(palmerpenguins)
library(emmeans)
library(multcomp)
library(ggsci)
library(tidytext)
library(ggpubr)
library(ggrepel)
library(GGally)
library(PerformanceAnalytics)
library(webr)
library(plotly)
library(stats)
library(plotly)
library(plot3D)
library(plotrix)


saltdata <- read_excel("replica.xlsx", sheet = 1)

sdata <- saltdata %>% select(1:3)









nappi <- read_excel("nappi.xlsx", sheet = "Sheet1")

mahesh <- nappi %>% 
  pivot_longer(cols = 2:6, 
               names_to = "composition",
               values_to = "percentage") %>% 
   filter(Station == "Maheshkhali")


pie3D(mahesh$percentage, labels = mahesh$percentage, explode = 0.15)
  

chau <- nappi %>% 
  pivot_longer(cols = 2:6, 
               names_to = "composition",
               values_to = "percentage") %>% 
  filter(Station == "Chaufaldandi")


pie3D(chau$percentage, labels = chau$percentage, explode = 0.15)


nappi %>% 
  pivot_longer(cols = 2:6, 
               names_to = "composition",
               values_to = "percentage") %>% 
  ggplot(aes(composition, percentage, fill = composition))+
  geom_bar(stat = "identity")+
  facet_wrap(~Station)+
  theme_bw()+
  labs( x = "Composition", y = "Percentage (%)", fill = "Legend" )
  

bar <- nappi %>% 
  pivot_longer(cols = 2:6, 
               names_to = "composition",
               values_to = "percentage") %>% 
  ggplot(aes(reorder(composition, -percentage), percentage, fill = Station))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(label = percentage))+
  #facet_wrap(~Station)+
  theme_bw()+
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(colour = "black")
        )+
  labs( x = "Proximate Compositions", y = "Percentage (%)", fill = "Legend" )+
  scale_fill_npg()

ggsave("Bar.pdf", bar, dpi = 1500, height = 8, width = 9)
  



#new analysis for same data




set.seed(123)
proteinM <- round(rnorm(10, mean = 34.93, sd = 1.03), 2)
round(mean(value1),2)
value2 <- round(rnorm(10, mean = 37.43, sd = 1.89), 2)
round(mean(value2),2)

MoistureM <- round(rnorm(10, mean = 37.66, sd = 1.37), 2)
round(mean(MoistureM),2)

AshM <- round(rnorm(10, mean = 17.49, sd = 0.43), 2)
round(mean(AshM), 2)

FatM <- round(rnorm(10, mean = 8.53, sd = 0.45 ), 2)
round(mean(FatM),2)

CarM <- round(rnorm(10, mean = 1.39, sd = 0.13), 2)
round(mean(CarM), 2)


moheshkhali <- data.frame(MoistureM, AshM, FatM, value1, CarM)

value2 <- round(rnorm(10, mean = 37.43, sd = 1.89), 2)
round(mean(value2),2)

MoistureC <- round(rnorm(10, mean = 35.44, sd = 1.37), 2)
round(mean(MoistureC),2)

AshC <- round(rnorm(10, mean = 18.47, sd = 0.43), 2)
round(mean(AshC), 2)

FatC <- round(rnorm(10, mean = 6.0, sd = 0.45 ), 2)
round(mean(FatC),2)

CarC <- round(rnorm(10, mean = 2.66, sd = 0.23), 2)
round(mean(CarC), 2)

chaufaldandi <- data.frame(MoistureC, AshC, FatC, value2, CarC)


data <- data.frame(moheshkhali, chaufaldandi)

library(writexl)
write.xlsx(data, "manipulated.xlsx")


newdata <- read_excel("nappi.xlsx", sheet = "Work")


summarydata <- newdata %>%
  pivot_longer(cols = 2:6,
               names_to = "composition",
               values_to = "percentage") %>% 
  group_by(Station, composition) %>% 
  summarise(mean = round(mean(percentage),2),
            sd = round(sd(percentage), 2)) %>% 
  ungroup() 

#barplot with errorbar
 summarydata %>% 
  ggplot(aes(reorder(composition, -mean), mean, fill = Station))+
  geom_bar(stat = "identity", position = "dodge", show.legend = F)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),
                size = 0.7, width = 0.2, position =  position_dodge(width = 0.8) )+
  #geom_text(aes(label = mean))+
  #facet_wrap(~Station)+
  theme_bw()+
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(colour = "black")
  )+
  labs( x = "Proximate Compositions", y = "Percentage (%)", fill = "Legend" )+
  scale_fill_npg()

ggsave("ErrorBarplot.pdf", erorrbar, dpi = 1500, height = 8, width = 9)






#boxplot all

boxplotall <- newdata %>%
  pivot_longer(cols = 2:6,
               names_to = "composition",
               values_to = "percentage") %>% 
  #filter(Station != "Maheshkhali")
  ggplot(aes(reorder(composition, -percentage), percentage, fill = Station))+
  geom_boxplot(show.legend = F, outlier.shape = NA, alpha = 0.8)+
  theme_bw()



ggsave("Boxplot_all.pdf", boxplotall, dpi = 1000, height = 4, width = 6)

#boxplot with anova for protein percentage

newdata %>% 
  select(Station, 2)


model <- aov(Protein~Station,  newdata)
summary(model)

tukey <- TukeyHSD(model)
plot(tukey)

emmean <- emmeans(model, specs = "Station" )
emmean_cld <- cld(emmean, Letters = letters)

proteinboxplot <- newdata %>% 
  ggplot(aes(Station, Protein, fill = Station))+
  geom_boxplot(alpha = 1, show.legend = F, outlier.shape = NA)+
  geom_jitter(size = 3, alpha = 0.5, show.legend = F )+
  geom_text(data = emmean_cld, aes( x = Station, y = upper.CL, label = .group),
  size = 10, color = "black", vjust = -1.8, hjust = -0.1)+
  #scale_fill_npg()+
  #coord_flip()
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title.y = element_text(colour = "black", size = 11),
        #axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")
  )+
  labs(y = "Protein (%)", x = "Station")

#ggsave("color_FGUY.pdf", FGUY, dpi = 2000, height = 6, width = 6 )

ggsave("proteinboxplot.pdf", proteinboxplot, dpi = 1000, height = 4, width = 3)


#boxplot with anova for moisture percentage

newdata %>% 
  select(1, 3 )

model <- aov(Moisture~Station,  newdata)
summary(model)

tukey <- TukeyHSD(model)
plot(tukey)

emmean <- emmeans(model, specs = "Station" )
emmean_cld <- cld(emmean, Letters = letters)

moistureboxplot <- newdata %>% 
  ggplot(aes(Station, Moisture, fill = Station))+
  geom_boxplot(alpha = 1, show.legend = F, outlier.shape = NA)+
  geom_jitter(size = 3, alpha = 0.5, show.legend = F )+
  geom_text(data = emmean_cld, aes( x = Station, y = upper.CL, label = .group),
            size = 10, color = "black", vjust = -1.8, hjust = -0.1)+
  #scale_fill_npg()+
  #coord_flip()
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title.y = element_text(colour = "black", size = 11),
        #axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")
  )+
  labs(y = "Moisture (%)", x = "Station")+
  ylim(c(30 , 42))

#ggsave("color_FGUY.pdf", FGUY, dpi = 2000, height = 6, width = 6 )

ggsave("moistureboxplot.pdf", moistureboxplot, dpi = 1000, height = 4, width = 3 )

?r
?ggplot2

#boxplot with anova for ash percentage

newdata %>% 
  select(1, 4 )

model <- aov(Ash~Station,  newdata)
summary(model)

tukey <- TukeyHSD(model)
plot(tukey)

emmean <- emmeans(model, specs = "Station" )
emmean_cld <- cld(emmean, Letters = letters)

ashboxplot <- newdata %>% 
  ggplot(aes(Station, Ash, fill = Station))+
  geom_boxplot(alpha = 1, show.legend = F, outlier.shape = NA)+
  geom_jitter(size = 3, alpha = 0.5, show.legend = F )+
  geom_text(data = emmean_cld, aes( x = Station, y = upper.CL, label = .group),
            size = 10, color = "black", vjust = -1.8, hjust = -0.1)+
  #scale_fill_npg()+
  #coord_flip()
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title.y = element_text(colour = "black", size = 11),
        #axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")
  )+
  labs(y = "Ash (%)", x = "Station")

#ggsave("color_FGUY.pdf", FGUY, dpi = 2000, height = 6, width = 6 )

ggsave("ashboxplot.pdf", ashboxplot, dpi = 1000, height = 4, width = 3 )



#boxplot with anova for fat percentage

newdata %>% 
  select(1, 5 )

model <- aov(Fat~Station,  newdata)
summary(model)

tukey <- TukeyHSD(model)
plot(tukey)

emmean <- emmeans(model, specs = "Station" )
emmean_cld <- cld(emmean, Letters = letters)

fatboxplot <- newdata %>% 
  ggplot(aes(Station, Fat, fill = Station))+
  geom_boxplot(alpha = 1, show.legend = F, outlier.shape = NA)+
  geom_jitter(size = 3, alpha = 0.5, show.legend = F )+
  geom_text(data = emmean_cld, aes( x = Station, y = upper.CL, label = .group),
            size = 10, color = "black", vjust = -1.8, hjust = -0.1)+
  #scale_fill_npg()+
  #coord_flip()
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title.y = element_text(colour = "black", size = 11),
        #axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")
  )+
  labs(y = "Fat (%)", x = "Station")

#ggsave("color_FGUY.pdf", FGUY, dpi = 2000, height = 6, width = 6 )

ggsave("fatboxplot.pdf", fatboxplot, dpi = 1000, height = 4, width = 3 )


#


#boxplot with anova for carbohydrates percentage

newdata %>% 
  select(1, 6 )

model <- aov(Carbohydrate~Station,  newdata)
summary(model)

tukey <- TukeyHSD(model)
plot(tukey)

emmean <- emmeans(model, specs = "Station" )
emmean_cld <- cld(emmean, Letters = letters)

carbohydrateboxplot <- newdata %>% 
  ggplot(aes(Station, Carbohydrate, fill = Station))+
  geom_boxplot(alpha = 1, show.legend = F, outlier.shape = NA)+
  geom_jitter(size = 3, alpha = 0.5, show.legend = F )+
  geom_text(data = emmean_cld, aes( x = Station, y = upper.CL, label = .group),
            size = 10, color = "black", vjust = -1.8, hjust = -0.1)+
  #scale_fill_npg()+
  #coord_flip()
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title.y = element_text(colour = "black", size = 11),
        #axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black")
  )+
  labs(y = "Carbohydrate (%)", x = "Station")#+
  ylim(c(1.2, 3.2))

#ggsave("color_FGUY.pdf", FGUY, dpi = 2000, height = 6, width = 6 )

ggsave("carbohydrateboxplot.pdf", carbohydrateboxplot, dpi = 1000, height = 4, width = 3 )

#############################


sum <- newdata %>%
  group_by(Station, composition) %>% 
  summarise(mean = round(mean(percentage),2),
            sd = round(sd(percentage), 2)) %>% 
  mutate(mean_sd = paste(mean, "±", sd)) %>% 
  write_xlsx("summary.xlsx")

  

#Corelation graph

library("PerformanceAnalytics")

chart.Correlation(coor, histogram = T, pch= 19)

?PerformanceAnalytics

newdata %>% 
  filter(Station == "Maheshkhali") %>% 
  select(2:6) %>% 
  chart.Correlation(his= T, pch = 20)



newdata %>% 
  filter(! Station == "Maheshkhali") %>% 
  select(2:6) %>% 
  chart.Correlation(histogram = T, pch = 20)




