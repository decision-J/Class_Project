rm(list=ls())

library(sp)
library(gstat)
library(fields)
library(classInt)
library(maps)
library(dplyr)

library(ggmap)
library(ggplot2)
library(rgdal)

dir = "C:/Users/JYW/OneDrive - 연세대학교 (Yonsei University)/정해용/2020 1학기/시공간분석/Final project/"
# Should change the directory

load(paste0(dir, "dataset/STDA_proj.Rdata"))

#### Date에 따른 전파 현황 plot <- TARGET
MAP_invasion_date <- ggplot(totalDf,
                      aes(x = long,
                          y = lat,
                          group = group)) +
  geom_polygon(color = "black", aes(fill = as.factor(confirmed_date))) +
  scale_fill_manual(values = rev(heat.colors(n = length(unique(totalDf$confirmed_date)))),
                    name = "First date") +
  coord_fixed(1.3)
x11()
MAP_invasion_date

#### 집단발병 case point 찍기 <- DATA
MAP_point_case <-
  ggplot(data = totalDf,add=TRUE,
         mapping = aes(x = long,
                       y = lat,
                       group = group,
                       color=data$confirmed)) +
  geom_polygon(fill = 'lightgrey',
               color = 'black') +
  geom_point(data = data,
             mapping = aes(x = longitude,
                           y = latitude,
                           group = confirmed)) +
  scale_color_gradientn(colours = rainbow(5), name = "confirmed")
x11()
MAP_point_case


#### 
source(paste0(dir, 'Rcode/STDA_PROJ_function.R'))
load(paste0(dir, "dataset/STDA_proj.Rdata"))
data = na.omit(data)

data$date <- as.numeric(as.Date(data$confirmed_date))

coord = cbind(data$longitude,data$latitude)
dates = data$date
confirmed = data$confirmed


# adaptive mcmc algorithm
# specify batches, batch length, target acceptance rate
amcmc=list(n.batch=2000,batch.length=20000,accept.rate=0.15)
out.grad = localgrad(verbose=T,dates=dates,covar = confirmed,coord.longlat=coord,n.samp=200000,thin=5, Albers=c(29.5,45.5), amcmc=amcmc)
summary.localgrad(out.grad, show=T)

# posterior mean of parameters
par.samp = apply(out.grad$post.samp,2,mean)
par.samp
as.Date(par.samp[1], origin="1970-01-01")

par.samp.CI <- t(apply(out.grad$post.samp, 2, function(x) {quantile(x, c(0.025, 0.975)) }))
par.samp.CI
as.Date(par.samp.CI[1,], origin="1970-01-01")

# test for long range jump
out.jump = jump.scan(dates=dates, coord.longlat=coord, params = par.mean, Albers=c(29.5,45.5) )

# gradient plot
plotgrad(out.grad, merge.data=data, database=totalDf)
plotgrad_plus(out.grad, map = MAP_invasion_date, merge.data=data, database=totalDf)




#### Simulation to Korea
dir = "C:/Users/JYW/OneDrive - 연세대학교 (Yonsei University)/정해용/2020 1학기/시공간분석/Final project/dataset"

###########################################################
data.kor = read.csv(paste0(dir , "/Case.csv"), as.is = T)
data.kor = subset(data.kor, select=c("province","confirmed","latitude","longitude","city","infection_case"))
data.kor = data.kor[data.kor[,3]!="-",]
data.kor$confirmed = as.numeric(data.kor$confirmed)
data.kor$latitude = as.numeric(data.kor$latitude)
data.kor$longitude = as.numeric(data.kor$longitude)

proj = mapproject(data.kor[,4], data.kor[,3], projection="albers", parameters=c(29.5,45.5), orientation=NULL)
data.kor$longitude = proj$x[1:nrow(data.kor)]
data.kor$latitude = proj$y[1:nrow(data.kor)]

patient.kor = read.csv(paste0(dir , "/patientInfo.csv"), as.is = T)
patient.kor = subset(patient.kor, select=c("patient_id","province","city", "confirmed_date"))
patient.kor = patient.kor[patient.kor[,3]!="",]
patient.kor = patient.kor[patient.kor[,4]!="",]

###########################################################
total <- 
  rgdal::readOGR(
    dsn = dir, 
    layer = 'total',
    encoding = 'CP949')

totalDf.kor <- fortify(model = total)
total@data$id <- rownames(x = total@data)

###########################################################
totalDf.kor = totalDf.kor %>% left_join(total@data, by = "id")
totalDf.kor <- totalDf.kor[order(totalDf.kor$id, totalDf.kor$order), ]
table(totalDf.kor$CTP_ENG_NM) %>% sort()

colnames(totalDf.kor) <- c("long", "lat", "order", "hole" ,     
                           "piece", "id", "group", "CTPRVN_CD", 
                           "province", "SIG_CD", "city")
totalDf.kor$city = as.character(totalDf.kor$city)
totalDf.kor$province = as.character(totalDf.kor$province)

proj = mapproject(totalDf.kor$long, totalDf.kor$lat, projection="albers", parameters=c(29.5,45.5), orientation=NULL)
totalDf.kor$long = proj$x
totalDf.kor$lat = proj$y

###########################################################
# data.kor join

# Add 1 
patient.kor1 = aggregate(patient.kor$patient_id, by=list(patient.kor$city,patient.kor$province), FUN=length)
colnames(patient.kor1) <- c("city","province","confirmed")
patient.kor1 = patient.kor1[!duplicated(patient.kor1[,c('city','province')]),] # 중복 제거

totalDf.kor = totalDf.kor %>% left_join(patient.kor1, by=c("city","province"))

# Add 2
patient.kor$date = paste0(substr(patient.kor$confirmed_date, 6,7), substr(patient.kor$confirmed_date, 9,10))
patient.kor_date <- patient.kor %>% group_by(confirmed_date, city) 
patient.kor_date <- patient.kor_date[!duplicated(patient.kor_date$city), ]
patient.kor_date = subset(patient.kor_date, select=c("province","city", "confirmed_date"))

totalDf.kor = totalDf.kor %>% left_join(patient.kor_date, by=c("city","province"))

# Add 3
data.kor1 = aggregate(data.kor$confirmed, by=list(data.kor$city,data.kor$province), FUN=sum)
colnames(data.kor1) <- c("city","province","case")
totalDf.kor = totalDf.kor %>% left_join(data.kor1, by=c("province","city"))
totalDf.kor = totalDf.kor[totalDf.kor$province != "Daegu"& totalDf.kor$province != "Gyeongsangbuk-do"&totalDf.kor$province !=  "Gyeongsangnam-do",]

# Add 4
data.kor = data.kor %>% left_join(patient.kor_date, by=c("province","city"))
############################

data.kor = data.kor[data.kor$province != "Daegu"& data.kor$province != "Gyeongsangbuk-do"&data.kor$province !=  "Gyeongsangnam-do",]
data.kor = na.omit(data.kor)

data.kor$date <- as.numeric(as.Date(data.kor$confirmed_date))

coord.new = cbind(data.kor$longitude,data.kor$latitude)
dates.new = data.kor$date
confirmed.new = data.kor$confirmed


# adaptive mcmc algorithm
# specify batches, batch length, target acceptance rate
amcmc=list(n.batch=2000,batch.length=20000,accept.rate=0.15)
out.grad2 = localgrad(dates=dates,covar=confirmed.new, coord.longlat=coord,newcoord.longlat=coord.new,n.samp=200000,thin=5, amcmc=amcmc)
summary.localgrad(out.grad2, show=T)


# posterior mean of parameters
par.mean2 = apply(out.grad2$post.samp,2,mean)

korea <- map_data(map = 'world', region = c('South Korea'))
korea <- korea[-grep("Do", korea$subregion),]

proj = mapproject(korea$long, korea$lat, projection="albers", parameters=c(29.5,45.5), orientation=NULL)
korea$long = proj$x
korea$lat = proj$y

# gradient plot
plotgrad2(out.grad2, merge.data=data.kor, database=korea)

MAP_invasion_date.kor <- ggplot(totalDf.kor,
                            aes(x = long,
                                y = lat,
                                group = group)) +
  geom_polygon(color = "black", aes(fill = as.factor(confirmed_date))) +
  scale_fill_manual(values = rev(heat.colors(n = length(unique(totalDf.kor$confirmed_date)))),
                    name = "First date") +
  coord_fixed(1.3)
x11()
MAP_invasion_date.kor