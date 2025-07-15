library("readxl")
library(dplyr)
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
setwd(path)

temperature<-read.csv("data/field_data/temp_maize.csv")

mintemp<-temperature[temperature$Type=="Min_Daily_Temp",]
mintemp<-mintemp[,-3]
colnames(mintemp)[3]<-"minTemperature"

maxtemp<-temperature[temperature$Type=="Max_Daily_Temp",]
maxtemp<-maxtemp[,-3]
colnames(maxtemp)[3]<-"maxTemperature"

temperature<-merge(mintemp[,-2],maxtemp[,-2],by="Date")
temperature$meanTemperature<-rowMeans(temperature[,2:3])


temperature$Date <- as.Date(temperature$Date)


get_avg_temp <- function(harvest_date) {
  temperature %>%
    filter(Date >= (harvest_date - 7) & Date <= harvest_date) %>%
    summarise(avg_temp = mean(meanTemperature, na.rm = TRUE)) %>%
    pull(avg_temp)
}




AQcurve2<-read.csv("data/processed_data/Testing_AQcurves_years22&23_plot.csv")
AQcurve2$Plot<-sub("_.*", "", AQcurve2$Plot_rep)
AQcurve2<-AQcurve2[,-3]


AQcurve1a<-read.csv("data/processed_data/Testing_Asat21_plot.csv")
AQcurve1a$Plot<-sub("_.*", "", AQcurve1a$Plot_rep)
AQcurve1a<-AQcurve1a[,-3]

AQcurve1b<-read.csv("data/processed_data/Training_Asat21_plot.csv")
AQcurve1b$Plot<-sub("_.*", "", AQcurve1b$Plot_rep)
AQcurve1b<-AQcurve1b[,-3]
AQcurve1<-rbind(AQcurve1a,AQcurve1b)


years<-c(2021,2022,2023)
for (y in 1:3){
  harv_dates<-read_excel(paste0("data/field_data/",years[y],"_sample_dates.xlsx"))
  # harv_dates$`Sampling date`<-gsub("/","-",harv_dates$`Sampling date`)
  colnames(harv_dates)<-c("Plot","Date")
  
  harv_dates$Date <- as.Date(harv_dates$Date)
  df<-merge(harv_dates,temperature,by="Date")
  
  df_week <- harv_dates %>%
    rowwise() %>%
    mutate(weekly_mean = get_avg_temp(Date)) %>%
    ungroup()
  
  
  if (years[y]==2021){
    AQcurve<-AQcurve1
  }else{
    AQcurve<-AQcurve2
  }
  
  temp<-merge(df,AQcurve[AQcurve$Year==years[y],],by="Plot")
  temp<-temp[order(temp$Accession),]
  colnames(temp)[colnames(temp)=="A_sat"]<-"PAR_1800"
  
  
  temp.w<-merge(df_week,AQcurve[AQcurve$Year==years[y],],by="Plot")
  temp.w<-temp.w[order(temp.w$Accession),]
  colnames(temp.w)[colnames(temp.w)=="A_sat"]<-"PAR_1800"
  
  
  temp_ave<-aggregate(temp, by=list(temp$Accession), "mean")
  colnames(temp_ave)[1]="Accession"
  
  temp_ave.w<-aggregate(temp.w, by=list(temp.w$Accession), "mean")
  colnames(temp_ave.w)[1]="Accession"
  
  
  temp_ave<-temp_ave[,c(1,4,5,6,7,9)]
  temp_ave.w<-temp_ave.w[,c(1,4,5,7)]

  
  # temp_ave<-aggregate(temp, by=list(temp$Plot,temp$Accession), "mean")
  # colnames(temp_ave)[2]="Accession"
  # temp_ave<-aggregate(temp_ave, by=list(temp_ave$Accession), "mean")
  # temp_ave<-temp_ave[,c(1,6,7,8,10,11)]

  if (y==1){
    Asat_field<-temp
    Asat_field_ave<-temp_ave
    
    Asat_field.w<-temp.w
    Asat_field_ave.w<-temp_ave.w
  }else{
    Asat_field<-rbind(Asat_field,temp[,1:8])
    Asat_field_ave<-rbind(Asat_field_ave,temp_ave)
    
    Asat_field.w<-rbind(Asat_field.w,temp.w[,1:6])
    Asat_field_ave.w<-rbind(Asat_field_ave.w,temp_ave.w)
  }
}



write.csv(Asat_field,file="data/processed_data/Testing_Asat3years_fieldcond_plot.csv",row.names = F)
write.csv(Asat_field_ave,file="data/processed_data/Testing_Asat3years_fieldcond_accession.csv",row.names = F)


write.csv(Asat_field.w,file="data/processed_data/Testing_Asat3years_weeklyTEMP_plot.csv",row.names = F)
write.csv(Asat_field_ave.w,file="data/processed_data/Testing_Asat3years_weeklyTEMP_accession.csv",row.names = F)


############# ACI 21 ######################
curve1<-read.csv("data/processed_data/Training_ACI21_plot.csv")
curve1$Plot<-sub("_.*", "", curve1$Plot_rep)
curve1<-curve1[,-3]

curve2<-read.csv("data/processed_data/Testing_ACI21_plot.csv")
curve2$Plot<-sub("_.*", "", curve2$Plot_rep)
curve2<-curve2[,-3]

curve3<-rbind(curve1,curve2)

harv_dates<-read_excel(paste0("data/field_data/2021_sample_dates.xlsx"))
# harv_dates$`Sampling date`<-gsub("/","-",harv_dates$`Sampling date`)
colnames(harv_dates)<-c("Plot","Date")

df<-merge(harv_dates,temperature,by="Date")


temp<-merge(df,curve3[curve3$Year==2021,],by="Plot")
temp<-temp[order(temp$Accession),]

temp_ave<-aggregate(temp, by=list(temp$Accession), "mean")
colnames(temp_ave)[1]="Accession"

temp_ave<-temp_ave[,c(1,4,5,6,7,seq(9,19))]

write.csv(temp,file="data/processed_data/ACI21_all_fieldcond_plot.csv",row.names = F)
write.csv(temp_ave,file="data/processed_data/ACI21_all_fieldcond_accession.csv",row.names = F)


####################################################

fieldcond<-read.csv("./data/processed_data/ACI21_all_fieldcond_plot.csv")
fieldcond$Year<-as.factor(fieldcond$Year)


g<-ggplot(fieldcond, aes(y=CO2_1250,x=meanTemperature, color=Year)) +
  geom_point(size=1)  + theme_bw()+
  geom_smooth(method=lm)+
  theme(legend.position = "bottom",legend.text = element_text(face="bold"),legend.title=element_text(face="bold"),
        axis.text=element_text(face="bold",color = "black"), axis.title=element_text(face="bold")) +
  scale_color_brewer(palette = "Dark2") +
  stat_cor(method="pearson",cor.coef.name = "r", label.y.npc='center')+
  labs(x="Temperature", y = "A")
