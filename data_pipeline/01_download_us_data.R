
# Build the raw U.S. monthly database used by the Matlab replication code.
#
# This script is sourced by build_us_database.R, which defines start.date,
# end.date, and indic.download, and loads the required packages.

fred_key <- Sys.getenv("FRED_API_KEY")
if (nzchar(fred_key)) {
  fredr_set_key(fred_key)
} else {
  message("FRED_API_KEY is not set. fredr will use any key already configured in the R session.")
}

start_date <- as.Date(start.date)
end_date   <- as.Date(end.date)
f <- function(ticker,freq){
  fredr(series_id = ticker,
        observation_start = start_date,observation_end = end_date,
        frequency = freq,aggregation_method = "avg")
}

list.variables <- c("FEDFUNDS","DTB4WK","DTB3","THREEFY1","THREEFY2","THREEFY5","THREEFY10",
                    "THREEFYTP5","THREEFYTP10",
                    #"THREEFY3","THREEFY5","THREEFY7","THREEFY10",
                    "CPIAUCSL","BBKMGDP","DFII5","DFII10","DGS30")

for(i in 1:length(list.variables)){
  data.var <- f(list.variables[i],"m")
  eval(parse(text = gsub(" ","",paste("data.var.frame = data.frame(date=data.var$date,",
                                      list.variables[i],"=data.var$value)",
                                      sep=""))))
  if(i==1){
    DATA = data.var.frame
  }else{
    DATA = merge(DATA,data.var.frame,by="date",all=TRUE)
  }
}

list.q.variables <- c("GDPPOT","GDPC1")
for(i in 1:length(list.q.variables)){
  data.var <- f(list.q.variables[i],"q")
  eval(parse(text = gsub(" ","",paste("data.var.frame = data.frame(date=data.var$date,",
                                      list.q.variables[i],"=data.var$value)",
                                      sep=""))))
  data.var.frame$date <- as.Date(paste(format(data.var$date,"%Y"),"-",
                                       as.numeric(format(data.var$date,"%m"))+2,"-01",sep=""))
  DATA = merge(DATA,data.var.frame,by="date",all=TRUE)
}


DATA$z <- log(DATA$GDPC1/DATA$GDPPOT)
#DATA$z <- log(DATA$GDPC1/DATA$GDPPOT) - mean(DATA$z,na.rm = TRUE)

lag <- 1
DATA$pi <- NaN
DATA$pi[(lag+1):dim(DATA)[1]] <- log(DATA$CPIAUCSL[(lag+1):dim(DATA)[1]]/
                                       DATA$CPIAUCSL[1:(dim(DATA)[1]-lag)])

DATA$dy <- log(1 + DATA$BBKMGDP/12/100)

plot(DATA$date,DATA$z,pch=19,col="black");grid()


# Download Survey of Professional Forecasters data: ============================

# 10 year GDP growth------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_rgdp10_level.xlsx",
                "data/Data_US/mean_rgdp10_level.xlsx")
}
SPF <- readxl::read_xlsx(path="data/Data_US/mean_rgdp10_level.xlsx",na="#N/A")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
indic1st <- which(SPF$date==start.date)
SPF      <- SPF[indic1st:dim(SPF)[1],]
SPF.GDP10<- data.frame(date=SPF$date,RGDP10=SPF$RGDP10)
# Shift dates:
SPF.GDP10$date <- as.Date(paste(format(SPF.GDP10$date,"%Y"),"-",
                                as.numeric(format(SPF.GDP10$date,"%m"))+1,"-01",sep=""))


# 1 year GDP growth-------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_rgdp_level.xlsx",
                "data/Data_US/mean_rgdp_level.xlsx", mode = "wb")
}
SPF <- readxl::read_xlsx(path="data/Data_US/mean_rgdp_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,RGDP1=100*log(as.numeric(SPF$RGDP6)/as.numeric(SPF$RGDP2)))
indic1st <- which(SPF$date==start.date)
SPF <- SPF[indic1st:dim(SPF)[1],]
SPF.GDP<- data.frame(date=SPF$date,RGDP1=SPF$RGDP1)
SPF.GDP$date <- as.Date(paste(format(SPF.GDP$date,"%Y"),"-",
                              as.numeric(format(SPF.GDP$date,"%m"))+1,"-01",sep=""))

# 1 year CPI -------------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_cpi_level.xlsx",
                "data/Data_US/mean_cpi_level.xlsx")
}
SPF <- readxl::read_xlsx(path="data/Data_US/mean_cpi_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,CPI1=as.numeric(SPF$CPI6))
indic1st <- which(SPF$date==start.date)
SPF <- SPF[indic1st:dim(SPF)[1],]
SPF.CPI <- data.frame(date=SPF$date,CPI1=SPF$CPI1)
SPF.CPI$date <- as.Date(paste(format(SPF.CPI$date,"%Y"),"-",
                              as.numeric(format(SPF.CPI$date,"%m"))+1,"-01",sep=""))

# 10 year CPI-------------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_cpi10_level.xlsx",
                "data/Data_US/mean_cpi10_level.xlsx")
}
SPF <- readxl::read_xlsx(path="data/Data_US/mean_cpi10_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,CPI10=as.numeric(SPF$CPI10))
indic1st <- which(SPF$date==start.date)
SPF <- SPF[indic1st:dim(SPF)[1],]
SPF.CPI10<- data.frame(date=SPF$date,CPI10=SPF$CPI10)
SPF.CPI10$date <- as.Date(paste(format(SPF.CPI10$date,"%Y"),"-",
                                as.numeric(format(SPF.CPI10$date,"%m"))+1,"-01",sep=""))


# 1 year TBILL -----------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_tbill_level.xlsx",
                "data/Data_US/mean_tbill_level.xlsx")
}
SPF <- readxl::read_xlsx(path="data/Data_US/mean_tbill_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,BILL1=as.numeric(SPF$TBILL6))
indic1st <- which(SPF$date==start.date)
SPF <- SPF[indic1st:dim(SPF)[1],]
SPF.BILL <- data.frame(date=SPF$date,BILL1=SPF$BILL1)
SPF.BILL$date <- as.Date(paste(format(SPF.BILL$date,"%Y"),"-",
                               as.numeric(format(SPF.BILL$date,"%m"))+1,"-01",sep=""))


# 10 year TBILL-----------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_bill10_level.xlsx",
                "data/Data_US/mean_bill10_level.xlsx")
}
SPF <- readxl::read_xlsx(path="data/Data_US/mean_bill10_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,BILL10=as.numeric(SPF$BILL10))
indic1st <- which(SPF$date==start.date)
SPF <- SPF[indic1st:dim(SPF)[1],]
SPF.BILL10<- data.frame(date=SPF$date,BILL10=SPF$BILL10)
SPF.BILL10$date <- as.Date(paste(format(SPF.BILL10$date,"%Y"),"-",
                                 as.numeric(format(SPF.BILL10$date,"%m"))+1,"-01",sep=""))


# LW estimates -----------------------------------------------------------------
if(indic.download==1){
  download.file("https://www.newyorkfed.org/medialibrary/media/research/economists/williams/data/Laubach_Williams_current_estimates.xlsx",
                "data/Data_US/Laubach_Williams_current_estimates.xlsx", mode = "wb")
  download.file("https://www.newyorkfed.org/medialibrary/media/research/economists/williams/data/Holston_Laubach_Williams_current_estimates.xlsx",
                "data/Data_US/Holston_Laubach_Williams_current_estimates.xlsx", mode = "wb")
}
LW <- readxl::read_xlsx(path="data/Data_US/Laubach_Williams_current_estimates.xlsx",
                        sheet="data",skip=5)
LW$date <- as.Date(substr(LW$Date,1,10))
LW <- data.frame(date=LW$date,rstarLW=as.numeric(LW$rstar...8))
indic1st <- which(LW$date==start.date)
LW <- LW[indic1st:dim(LW)[1],]
LW$date <- as.Date(paste(format(LW$date,"%Y"),"-",
                         as.numeric(format(LW$date,"%m"))+2,"-01",sep=""))

HLW <- readxl::read_xlsx(path="data/Data_US/Holston_Laubach_Williams_current_estimates.xlsx",
                        sheet="HLW Estimates",skip=5)
HLW$date <- as.Date(substr(HLW$Date,1,10))
HLW <- data.frame(date=HLW$date,rstarHLW=as.numeric(HLW$US...11))
indic1st <- which(HLW$date==start.date)
HLW <- HLW[indic1st:dim(HLW)[1],]
HLW$date <- as.Date(paste(format(HLW$date,"%Y"),"-",
                         as.numeric(format(HLW$date,"%m"))+2,"-01",sep=""))


# Nominal yield data from GSW 2006:
if(indic.download==1){
  download.file("https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv",
                "data/Data_US/feds200628.csv")
}
DAT <- csv.get(file = "data/Data_US/feds200628.csv", skip = 9)
library(tidyverse)
library(zoo)
DAT$Date <- as.Date(DAT$Date)
DAT <- arrange(DAT,Date)
DAT$mdate <- paste(format(DAT$Date,"%Y"),"-",
                   format(DAT$Date,"%m"),sep="")
DAT_mthly <- DAT %>%
  group_by(mdate) %>% #filter(Date == max(Date,na.rm = TRUE))
  summarise_all(function(x){mean(x,na.rm=TRUE)})
DAT_mthly$Date <- as.Date(
  paste(format(DAT_mthly$Date,"%Y"),"-",
        format(DAT_mthly$Date,"%m"),"-01",sep=""))
DAT_GSW_nom <- data.frame(date=DAT_mthly$Date,
                          SVENY01=DAT_mthly$SVENY01,
                          SVENY02=DAT_mthly$SVENY02,
                          SVENY03=DAT_mthly$SVENY03,
                          SVENY04=DAT_mthly$SVENY04,
                          SVENY05=DAT_mthly$SVENY05,
                          SVENY06=DAT_mthly$SVENY06,
                          SVENY07=DAT_mthly$SVENY07,
                          SVENY08=DAT_mthly$SVENY08,
                          SVENY09=DAT_mthly$SVENY09,
                          SVENY10=DAT_mthly$SVENY10,
                          SVENY11=DAT_mthly$SVENY11,
                          SVENY12=DAT_mthly$SVENY12,
                          SVENY13=DAT_mthly$SVENY13,
                          SVENY14=DAT_mthly$SVENY14,
                          SVENY15=DAT_mthly$SVENY15,
                          SVENY16=DAT_mthly$SVENY16,
                          SVENY17=DAT_mthly$SVENY17,
                          SVENY18=DAT_mthly$SVENY18,
                          SVENY19=DAT_mthly$SVENY19,
                          SVENY20=DAT_mthly$SVENY20,
                          SVENY21=DAT_mthly$SVENY21,
                          SVENY22=DAT_mthly$SVENY22,
                          SVENY23=DAT_mthly$SVENY23,
                          SVENY24=DAT_mthly$SVENY24,
                          SVENY25=DAT_mthly$SVENY25,
                          SVENY26=DAT_mthly$SVENY26,
                          SVENY27=DAT_mthly$SVENY27,
                          SVENY28=DAT_mthly$SVENY28,
                          SVENY29=DAT_mthly$SVENY29,
                          SVENY30=DAT_mthly$SVENY30)

# Real yield data from GSW 2008:
if(indic.download==1){
  download.file("https://www.federalreserve.gov/data/yield-curve-tables/feds200805.csv",
                "data/Data_US/feds200805.csv")
}
DAT <- csv.get(file = "data/Data_US/feds200805.csv", skip = 18)
DAT$Date <- as.Date(DAT$Date)
DAT <- arrange(DAT,Date)
DAT$mdate <- paste(format(DAT$Date,"%Y"),"-",
                   format(DAT$Date,"%m"),sep="")
DAT_mthly <- DAT %>%
  group_by(mdate) %>% #filter(Date == max(Date,na.rm = TRUE))
  summarise_all(function(x){mean(x,na.rm=TRUE)})
DAT_mthly$Date <- as.Date(
  paste(format(DAT_mthly$Date,"%Y"),"-",
        format(DAT_mthly$Date,"%m"),"-01",sep=""))
DAT_GSW_real <- data.frame(date=DAT_mthly$Date,
                           TIPSY02=DAT_mthly$TIPSY02,
                           TIPSY03=DAT_mthly$TIPSY03,
                           TIPSY04=DAT_mthly$TIPSY04,
                           TIPSY05=DAT_mthly$TIPSY05,
                           TIPSY06=DAT_mthly$TIPSY06,
                           TIPSY07=DAT_mthly$TIPSY07,
                           TIPSY08=DAT_mthly$TIPSY08,
                           TIPSY09=DAT_mthly$TIPSY09,
                           TIPSY10=DAT_mthly$TIPSY10,
                           TIPSY11=DAT_mthly$TIPSY11,
                           TIPSY12=DAT_mthly$TIPSY12,
                           TIPSY13=DAT_mthly$TIPSY13,
                           TIPSY14=DAT_mthly$TIPSY14,
                           TIPSY15=DAT_mthly$TIPSY15,
                           TIPSY16=DAT_mthly$TIPSY16,
                           TIPSY17=DAT_mthly$TIPSY17,
                           TIPSY18=DAT_mthly$TIPSY18,
                           TIPSY19=DAT_mthly$TIPSY19,
                           TIPSY20=DAT_mthly$TIPSY20)


# # PTR from Bauer and Rudebusch: 
# pistar_PTR <- read.csv("Data/pistar_PTR.csv")
# pistar_PTR$date <- as.Date(
#   paste(pistar_PTR$Year,"-",3*pistar_PTR$Quarter,"-01",sep=""))
# DAT_pistar <- data.frame(date=pistar_PTR$date,
#                          PTR=pistar_PTR$pistar_PTR)

# Source: FRB/US https://www.federalreserve.gov/econres/files/data_only_package.zip
if(indic.download==1){
  ptr_zip <- file.path("data", "Data_US", "data_only_package.zip")
  download.file("https://www.federalreserve.gov/econres/files/data_only_package.zip",
                ptr_zip, mode = "wb")
  unzip(ptr_zip, files = "data_only_package/HISTDATA.TXT",
        exdir = file.path("data", "Data_US", "Data_Fed_PTR"))
  file.rename(file.path("data", "Data_US", "Data_Fed_PTR", "data_only_package", "HISTDATA.TXT"),
              file.path("data", "Data_US", "Data_Fed_PTR", "HISTDATA.TXT"))
  unlink(file.path("data", "Data_US", "Data_Fed_PTR", "data_only_package"), recursive = TRUE)
}
HISTDATA <- read.csv("data/Data_US/Data_Fed_PTR/HISTDATA.TXT")
HISTDATA$quarter <- as.numeric(substr(HISTDATA$OBS,6,6))
HISTDATA$date <- as.Date(paste(substr(HISTDATA$OBS,1,4),"-",
                               HISTDATA$quarter*3,"-01",sep=""))
DAT_pistar <- data.frame(date=HISTDATA$date,
                         PTR=HISTDATA$PTR)


# Synthetic breakeven inflation rates:
if(indic.download==1){
  download.file("https://newyorkfed.org/medialibrary/media/research/blog/groen_tips/Synthetic_TIPS_Breakeven_Rates.xlsx",
                "data/Data_US/Synthetic_TIPS_Breakeven_Rates.xlsx", mode = "wb")
}
Synthetic_BEIR <- read_excel("data/Data_US/Synthetic_TIPS_Breakeven_Rates.xlsx",
                             skip = 13)
Synthetic_BEIR$Date <- as.Date(Synthetic_BEIR$...1)
Synthetic_BEIR$Date <- as.Date(
  paste(as.numeric(format(Synthetic_BEIR$Date,"%Y")),"-",
        as.numeric(format(Synthetic_BEIR$Date,"%m")),"-01",sep=""))
Synthetic_BEIR <- data.frame(date=Synthetic_BEIR$Date,
                             RR10Y=Synthetic_BEIR$`10-yr real rate, backcast`,
                             RR20Y=Synthetic_BEIR$`20-yr real rate, backcast`)


# Merge all data frames:

DATA <- merge(DATA,SPF.GDP10,by="date",all=TRUE)
DATA <- merge(DATA,SPF.GDP,by="date",all=TRUE)
DATA <- merge(DATA,SPF.CPI10,by="date",all=TRUE)
DATA <- merge(DATA,SPF.CPI,by="date",all=TRUE)
DATA <- merge(DATA,SPF.BILL10,by="date",all=TRUE)
DATA <- merge(DATA,SPF.BILL,by="date",all=TRUE)
DATA <- merge(DATA,DAT_GSW_nom,by="date",all=TRUE)
DATA <- merge(DATA,DAT_GSW_real,by="date",all=TRUE)
DATA <- merge(DATA,DAT_pistar,by="date",all=TRUE)
DATA <- merge(DATA,LW,by="date",all=TRUE)
DATA <- merge(DATA,Synthetic_BEIR,by="date",all=TRUE)

DATA <- DATA[complete.cases(DATA$DTB3),]


par(mfrow=c(1,2))
plot(DATA$date,DATA$z,pch=19,ylim=c(-.1,.04))
points(DATA$date,DATA$RGDP1)

# plot(DATA$date,DATA$z,type="l")
# lines(DATA$date,DATA$DFII10/100,col="red")
# lines(DATA$date,DATA$DFII5/100,col="blue")

# plot(DATA$date,DATA$RGDP1,type="l",ylim=c(-.02,.07))
# lines(DATA$date,DATA$DFII10/100,col="red")
# lines(DATA$date,DATA$DFII5/100,col="blue")

plot(DATA$date,DATA$pi,type="l")
points(DATA$date,DATA$CPI1/12/100,col="blue")
points(DATA$date,DATA$CPI10/12/100,col="red")

plot(DATA$date,DATA$dy,type="l")
points(DATA$date,DATA$RGDP1/12/100,col="blue")
points(DATA$date,DATA$RGDP10/12/100,col="red")

DATA$TRP10 <- DATA$THREEFY10 - DATA$BILL10

DATA$IRP10 <- (DATA$THREEFY10 - DATA$TIPSY10) - DATA$CPI10
plot(DATA$date,DATA$THREEFYTP10,type="l")
points(DATA$date,DATA$IRP10,col="red")
points(DATA$date,DATA$TRP10)

# Lower bound for nominal rates:
DATA$i_bar <- 0

save(DATA,LW,HLW,file="data/Data_US/data.Rda")
