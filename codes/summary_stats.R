
library("cowplot")
library("ggplot2")  

#read the data file

stat_sum <- read.csv("data/run_summary.csv", 
                     sep = ",", 
                     header = TRUE)

head(stat_sum)


Bases<-ggplot(data=stat_sum)+ 
  (aes(x =Run, y=Bases))+
  geom_boxplot(fill = "#C0C0C0", outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.2, size=0.3)+
  theme_bw()+
  labs( x= "Ion A530 chip",
        y= "Number of bases")+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "bottom")+ 
  theme(legend.title = element_blank())

Q20.<-ggplot(data=stat_sum)+ 
  (aes(x =Run, y=Q20.))+
  geom_boxplot(fill = "#C0C0C0", outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.2, size=0.3)+
  theme_bw()+
  labs( x= "Ion A530 chip",
        y= "Number of bases")+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "bottom")+ 
  theme(legend.title = element_blank())

Reads<-ggplot(data=stat_sum)+ 
  (aes(x =Run, y=Reads))+
  geom_boxplot(fill = "#C0C0C0", outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.2, size=0.3)+
  theme_bw()+
  labs( x= "Ion A530 chip",
        y= "Number of reads")+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "bottom")+ 
  theme(legend.title = element_blank())

Read_lenght<-ggplot(data=stat_sum)+ 
  (aes(x =Run, y=Read_lenght))+
  geom_boxplot(fill = "#C0C0C0", outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.2, size=0.3)+
  theme_bw()+
  labs( x= "Ion A530 chip",
        y= "Read length (bp)")+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "bottom")+ 
  theme(legend.title = element_blank())

Summ<-plot_grid(Bases+theme(legend.position="none"),
                 Q20.+theme(legend.position="none"), 
                 Reads+theme(legend.position="none"),
                 Read_lenght+theme(legend.position="none"),
                 align="vh",
                 labels = c("A", "B", "C", "D"),
                 hjust = -1,
                 vjust= 2,
                 nrow = 2)

Summ_final<-plot_grid(Summ, ncol = 1, rel_heights = c(0.8, .05))

Summ_final

summary(stat_sum)

ggsave(file="figures/figS1_summary_stat.pdf", 
       width=8.3, height=4.2, units="in", dpi=300)
