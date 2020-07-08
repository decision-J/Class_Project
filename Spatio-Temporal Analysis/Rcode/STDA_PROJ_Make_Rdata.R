library(sp)
library(gstat)
library(fields)
library(classInt)
library(maps)
library(dplyr)

library(ggmap)
library(ggplot2)
library(rgdal)

dir = "C:/Users/JYW/OneDrive - 연세대학교 (Yonsei University)/정해용/2020 1학기/시공간분석/Final project/dataset"
# Should change the directory


###########################################################
data = read.csv(paste0(dir , "/Case.csv"), as.is = T)
data = data[data$province == "Daegu"| data$province == "Gyeongsangbuk-do"|data$province ==  "Gyeongsangnam-do",]
data = subset(data, select=c("province","confirmed","latitude","longitude","city","infection_case"))
data = data[data[,3]!="-",]
data$confirmed = as.numeric(data$confirmed)
data$latitude = as.numeric(data$latitude)
data$longitude = as.numeric(data$longitude)

proj = mapproject(data$longitude, data$latitude, projection="albers", parameters=c(29.5,45.5), orientation=NULL)
data$longitude = proj$x
data$latitude = proj$y

patient = read.csv(paste0(dir , "/PatientInfo.csv"), as.is = T)
patient = patient[patient$province == "Daegu"| patient$province == "Gyeongsangbuk-do"|patient$province ==  "Gyeongsangnam-do",]
patient = subset(patient, select=c("patient_id","province","city", "confirmed_date"))
patient = patient[patient[,3]!="",]
patient = patient[patient[,4]!="",]

###########################################################
total <- 
  rgdal::readOGR(
    dsn = dir, 
    layer = 'total',
    encoding = 'CP949')

totalDf <- fortify(model = total)
total@data$id <- rownames(x = total@data)

###########################################################
totalDf = totalDf %>% left_join(total@data, by = "id")
totalDf <- totalDf[totalDf$CTP_ENG_NM == "Daegu"|totalDf$CTP_ENG_NM == "Gyeongsangbuk-do"|totalDf$CTP_ENG_NM == "Gyeongsangnam-do", ]
totalDf <- totalDf[totalDf$SIG_ENG_NM != "Ulleung-gun", ]
totalDf <- totalDf[order(totalDf$id, totalDf$order), ]
table(totalDf$CTP_ENG_NM) %>% sort()

colnames(totalDf) <- c("long", "lat", "order", "hole" ,     
                       "piece", "id", "group", "CTPRVN_CD", 
                       "province", "SIG_CD", "city")
totalDf$city = as.character(totalDf$city)
totalDf$province = as.character(totalDf$province)

proj = mapproject(totalDf$long, totalDf$lat, projection="albers", parameters=c(29.5,45.5), orientation=NULL)
totalDf$long = proj$x
totalDf$lat = proj$y

###########################################################
# Data join

# Add 1 
patient1 = aggregate(patient$patient_id, by=list(patient$city,patient$province), FUN=length)
colnames(patient1) <- c("city","province","confirmed")
patient1 = patient1[!duplicated(patient1[,c('city','province')]),] # 중복 제거

totalDf = totalDf %>% left_join(patient1, by=c("city","province"))

# Add 2
patient$date = paste0(substr(patient$confirmed_date, 6,7), substr(patient$confirmed_date, 9,10))
patient_date <- patient %>% group_by(confirmed_date, city) 
patient_date <- patient_date[!duplicated(patient_date$city), ]
patient_date = subset(patient_date, select=c("province","city", "confirmed_date"))

totalDf = totalDf %>% left_join(patient_date, by=c("city","province"))

# Add 3
data1 = aggregate(data$confirmed, by=list(data$city,data$province), FUN=sum)
colnames(data1) <- c("city","province","case")
totalDf = totalDf %>% left_join(data1, by=c("province","city"))

# Add 4
data = data %>% left_join(patient_date, by=c("province","city"))

###########################################################

save(data,patient,totalDf, file=paste0(dir, "/STDA_proj.Rdata"))

