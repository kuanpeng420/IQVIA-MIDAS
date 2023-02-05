#### IQVIA biosimilar paper analysis section
#### Kuan Peng 2021/9/13
pacman::p_load(data.table,tidyverse,tableone,writexl,readxl,comorbidity,broom,survival,ggplot2,ggpubr,
       plyr,foreign,tsModel,lmtest,Epi,mice,splines,vcd,AER,sandwich)
MIDAS <- read_excel("MIDAS.new.xlsx",sheet = 3,guess_max = 10000,progress = TRUE)
MIDAS <- data.table(MIDAS)
#filter essential drugs
Es_drug <- str_to_upper(c("infliximab"))
# Es_drug <- str_to_upper(c("ADALIMUMAB","INFLIXIMAB" ,"ETANERCEPT", "CERTOLIZUMAB PEGOL","GOLIMUMAB"))
ES <- MIDAS[`Molecule List`%in%Es_drug]
#filter out country
ES <- ES[!str_detect(Panel,regex("Russian",ignore_case = TRUE))]
ES <- ES[is.na(Biosimilar),Biosimilar:="N"]

Demo <- colnames(MIDAS[,Panel:`Relative protection expiry (year)`])

sale <- colnames(MIDAS)[str_detect(colnames(MIDAS),"MNF")]
ES_s_o <- ES[,.SD,by=c("Molecule List","Panel","NFC123","Corporation","ATC3","International Product" ),.SDcol=c(Demo,sale)]
ES_s_o <- ES_s_o[,.(quarter_sale=sapply(.SD,sum,na.rm = TRUE),year=substr(colnames(.SD),nchar(colnames(.SD))-3,nchar(colnames(.SD))),quater=substr(colnames(.SD),nchar(colnames(.SD))-5,nchar(colnames(.SD))-4)),by=c("Molecule List","Panel","NFC123","Corporation","ATC3","Biosimilar","International Product" ),.SDcol=sale]
countries <- c()
test <- (strsplit(ES_s_o$Panel, " "))
for (i in 1:length(test)) {
        tt <- unlist(test[[i]])[1]
        countries <- c(countries,tt)
}
ES_s_o <- ES_s_o[,nation:=countries]

Consumption <- colnames(MIDAS)[str_detect(colnames(MIDAS),"Standard")]
ES_c_o <- ES[,.SD,by=c("Molecule List","Panel","NFC123","Corporation","ATC3","International Product" ),.SDcol=c(Demo,Consumption)]
ES_c_o <- ES_c_o[,.(quarter_consum=sapply(.SD,sum,na.rm = TRUE),quater=substr(colnames(.SD),nchar(colnames(.SD))-5,nchar(colnames(.SD))-4),year=substr(colnames(.SD),nchar(colnames(.SD))-3,nchar(colnames(.SD)))),by=c("Molecule List","Panel","NFC123","Corporation","ATC3","Biosimilar","International Product" ),.SDcol=Consumption]
countries <- c()
test <- (strsplit(ES_c_o$Panel, " "))
for (i in 1:length(test)) {
        tt <- unlist(test[[i]])[1]
        countries <- c(countries,tt)
}
ES_c_o <- ES_c_o[,nation:=countries]
ES_p <- left_join(ES_c_o,ES_s_o,by=c("Biosimilar","Molecule List","nation","year","quater","Panel","ATC3","Corporation","NFC123","International Product" ))
ES_p <- ES_p[!quarter_sale==0][!quarter_consum==0]

ES_p[,`conversion factor`:=1]#standardized to 2020 https://www.usinflationcalculator.com/
ES_p[year==2010,`conversion factor`:=1.19]
ES_p[year==2011,`conversion factor`:=1.15]
ES_p[year==2012,`conversion factor`:=1.13]
ES_p[year==2013,`conversion factor`:=1.11]
ES_p[year==2014,`conversion factor`:=1.09]
ES_p[year==2015,`conversion factor`:=1.09]
ES_p[year==2016,`conversion factor`:=1.08]
ES_p[year==2017,`conversion factor`:=1.06]
ES_p[year==2018,`conversion factor`:=1.03]
ES_p[year==2019,`conversion factor`:=1.01]
ES_p[year==2020,`conversion factor`:=1]
ES_p[,quarter_sale:=round(quarter_sale*`conversion factor`)]
ES_p[,price:=round(quarter_sale/quarter_consum)]
es_remove <- ES_p[price<=10&str_detect(nation,"AUSTRALIA")]
a1 <- sum(es_remove$quarter_consum)
ES_p <- ES_p[price>10]# 1 dollar
# a2 <- sum(ES_p[str_detect(nation,"AUSTRALIA")]$quarter_consum)
a2 <- sum(ES_p$quarter_consum)

a1/a2*100#0.01% data removed
ES_p <- dcast(ES_p,`Molecule List`+nation+year+quater~Biosimilar,value.var = c("quarter_sale","quarter_consum"),fun.aggregate = sum)

eligible <- distinct(ES_p[!quarter_sale_Y==0][year<2019][!quarter_sale_N==0],`Molecule List`,nation)
#only want inflixmab data
eligible <- eligible[`Molecule List`=="INFLIXIMAB"]
ES_eli <- data.table()
for (i in 1:length(eligible$`Molecule List`)) {
        t <- ES_p[`Molecule List`=="INFLIXIMAB"&nation==eligible[[2]][i],]
        ES_eli <- rbind(ES_eli,t)
}
ES_eli[,qrt:=as.yearqtr(paste0(ES_eli$year,"-",ES_eli$quater))]

pop <- setDT(read_excel("Population.xlsx"))
pop <- melt(pop,measure.vars = c("2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"),variable.name = "year",value.name = "pop")
ES_eli <- merge(ES_eli,pop,by.x =c("nation","year"),by.y=c("Country","year") ,all.x = T)
ES_eli <- ES_eli[,`:=`(price_N=round(quarter_sale_N/quarter_consum_N),price_Y=round(quarter_sale_Y/quarter_consum_Y))]
ES_eli <- ES_eli[,`:=`(quarter_consum_N_std=quarter_consum_N/pop,quarter_consum_Y_std=quarter_consum_Y/pop)]
ES_eli <- ES_eli[,overall:=quarter_consum_N_std+quarter_consum_Y_std]
ES_eli <- ES_eli[,share:=round(quarter_consum_Y/(quarter_consum_Y+quarter_consum_N)*100,digits = 1)]
ES_eli <- ES_eli[,reduction:=round((price_N-price_Y)/price_N*100,digits = 1)]
ES_eli <- ES_eli[is.nan(reduction),reduction:=0]

##additional analysis
for (i in 1:8) {
    print(n[i])
    z <- ES_eli[,sum(overall),by=c("year","nation")][nation==n[i]]$V1[11]/ES_eli[,sum(overall),by=c("year","nation")][nation==n[i]]$V1[1]
    print(z)
}

##function to plot price trend----
res_plot <- function(model,y=800,data,datanew,ES_ITS,region){
        ggplot(data,aes(x=qrt,y=price_N))+geom_point(color="deepskyblue3",alpha=0.6)+geom_point(aes(x=qrt,y=price_Y),
        color="darkorange3",alpha=0.6)+geom_rect(mapping=aes(xmin=ES_ITS[time_after_inter==1,qrt],
        xmax=as.yearqtr("2020 Q4"),ymin=0,ymax=y),alpha=0.002)+
        geom_line(aes(y=predict(model,type="response",data)),color="springgreen4",size=0.7)+
        geom_line(aes(y=predict(model,type="response",datanew)),color="springgreen4",lty=2,size=0.7)+
        xlab("Quarter")+ylab("INF price (2020 USD)")+ggtitle(region)+scale_x_yearqtr(format = "%YQ%q", 
                    limits=c(as.yearqtr("2010 Q1"), as.yearqtr("2020 Q4")))+geom_label(
    label="Originator", 
    x=as.yearqtr(2010+1),
    y=y,
    label.padding = unit(0.2, "lines"), 
    label.size = 0.2,
    color = "black",fill="deepskyblue3")+geom_label(
    label="Biosimilar", 
    x=as.yearqtr(2010+4),
    y=y,
    label.padding = unit(0.2, "lines"),
    label.size = 0.2,
    color = "black",fill="darkorange3")+theme_classic()
}

res_plot2 <- function(model,y=800,data,datanew,ES_ITS,region){
        ggplot(data,aes(x=qrt,y=overall))+geom_point(color="firebrick3",alpha=0.8)+geom_rect(mapping=aes(xmin=ES_ITS[time_after_inter==1,qrt],
        xmax=as.yearqtr("2020 Q4"),ymin=0,ymax=y),alpha=0.002)+
        geom_line(aes(y=predict(model,type="response",data)),color="springgreen4",size=0.7)+
        geom_line(aes(y=predict(model,type="response",datanew)),color="springgreen4",lty=2,size=0.6)+
        xlab("Quarter")+ylab("Overall INF consumption/1000 inhabitants")+ggtitle(region)+scale_x_yearqtr(format = "%YQ%q", 
                    limits=c(as.yearqtr("2010 Q1"), as.yearqtr("2020 Q4")))+theme_classic()
    }

res_plot3 <- function(model,y=800,data,datanew,ES_ITS,region){
  ggplot(data,aes(x=qrt,y=price_N))+geom_point(color="deepskyblue3",alpha=0.6)+geom_rect(mapping=aes(xmin=ES_ITS[time_after_inter==1,qrt],xmax=as.yearqtr("2020 Q4"),ymin=0,ymax=y),alpha=0.002)+
    geom_line(aes(y=predict(model,type="response",data)),color="springgreen4",size=0.7)+
    geom_line(aes(y=predict(model,type="response",datanew)),color="springgreen4",lty=2,size=0.7)+
    xlab("Quarter")+ylab("INF price (2020 USD)")+ggtitle(region)+scale_x_yearqtr(format = "%YQ%q", 
                                                                                 limits=c(as.yearqtr("2010 Q1"), as.yearqtr("2020 Q4")))+theme_classic()
}

##newey atuocorrelation adjustment

newey_adjust <- function(x=model1){
        bw <- bwNeweyWest(x) # if set lag = NULL, floor(bwNeweyWest(x)) is used as lag
        
        est <- c(coef(x)["time"], coef(x)["Biosimilarintro"], coef(x)["time_after_inter"])
        se <- c(sqrt(diag(NeweyWest(x, prewhite = F, lag = ifelse(floor(bw)>3,3,floor(bw)))))["time"],
                sqrt(diag(NeweyWest(x, prewhite = F, lag = ifelse(floor(bw)>3,3,floor(bw)))))["Biosimilarintro"],
                sqrt(diag(NeweyWest(x, prewhite = F, lag = ifelse(floor(bw)>3,3,floor(bw)))))["time_after_inter"])
        
        lb <- est -1.96 * se
        ub <- est +1.96 * se
        
        table <- cbind(round(est, digits = 6), round(lb, digits = 6), round(ub, digits = 6), round(se, digits = 6))
        p_val <- coeftest(x, vcov. = NeweyWest(x, prewhite = F, lag = ifelse(floor(bw)>3,3,floor(bw)))) # Test SE get P-value  
        # list( BW = bw, Table = table, Summary = p_val)
        p_val
}

t.dif <- function(a,b,c){
   y1=predict(a,type="response",b)
y2=predict(a,type="response",c)
t.test(y1[!y1==y2],y2[!y1==y2]) 
}

# Price ITS----
ES_ITS1 <- ES_eli[nation=="HONG"]
ES_ITS1 <-ES_ITS1[,time:=1:nrow(ES_ITS1)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS1$Biosimilarintro==1)
ES_ITS1 <- ES_ITS1[,time_after_inter:=ifelse(time<=b,0,time-b)]
data1 <- ES_ITS1
datanew1 <- data1
datanew1$Biosimilarintro=0
model1 <- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data1)
newey_adjust(model1)  


ES_ITS2 <- ES_eli[nation=="INDIA"]
ES_ITS2 <-ES_ITS2[,time:=1:nrow(ES_ITS2)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS2$Biosimilarintro==1)
ES_ITS2 <- ES_ITS2[,time_after_inter:=ifelse(time<=b,0,time-b)]
data2 <- ES_ITS2
datanew2 <- data2
datanew2$Biosimilarintro=0
model2<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data2)
newey_adjust(model2)  


ES_ITS3 <- ES_eli[nation=="JAPAN"]
ES_ITS3 <-ES_ITS3[,time:=1:nrow(ES_ITS3)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS3$Biosimilarintro==1)
ES_ITS3 <- ES_ITS3[,time_after_inter:=ifelse(time<=b,0,time-b)]
data3 <- ES_ITS3
datanew3 <- data3
datanew3$Biosimilarintro=0
model3<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data3)
newey_adjust(model3)  


ES_ITS4 <- ES_eli[nation=="KOREA"]
ES_ITS4 <-ES_ITS4[,time:=1:nrow(ES_ITS4)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS4$Biosimilarintro==1)
ES_ITS4 <- ES_ITS4[,time_after_inter:=ifelse(time<=b,0,time-b)]
data4 <- ES_ITS4
datanew4 <- data4
datanew4$Biosimilarintro=0
model4<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data4)
newey_adjust(model4)  


ES_ITS5 <- ES_eli[nation=="MALAYSIA"]
ES_ITS5 <-ES_ITS5[,time:=1:nrow(ES_ITS5)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
ES_ITS5[year==2019&quarter_sale_Y==0,Biosimilarintro:=1]
b=44-sum(ES_ITS5$Biosimilarintro==1)
ES_ITS5 <- ES_ITS5[,time_after_inter:=ifelse(time<=b,0,time-b)]
data5 <- ES_ITS5
datanew5 <- data5
datanew5$Biosimilarintro=0
model5<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data5)
newey_adjust(model5)  


ES_ITS6 <- ES_eli[nation=="THAILAND"]
ES_ITS6 <-ES_ITS6[,time:=1:nrow(ES_ITS6)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS6$Biosimilarintro==1)
ES_ITS6 <- ES_ITS6[,time_after_inter:=ifelse(time<=b,0,time-b)]
data6 <- ES_ITS6
datanew6 <- data6
datanew6$Biosimilarintro=0
model6<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data6)
newey_adjust(model6)  


ES_ITS7 <- ES_eli[nation=="VIETNAM"]
ES_ITS7 <-ES_ITS7[,time:=1:nrow(ES_ITS7)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=34-sum(ES_ITS7$Biosimilarintro==1)
ES_ITS7 <- ES_ITS7[,time_after_inter:=ifelse(time<=b,0,time-b)]
data7 <- ES_ITS7
datanew7 <- data7
datanew7$Biosimilarintro=0
model7<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data7)
newey_adjust(model7)  



ES_ITS8 <- ES_eli[nation=="US"]
ES_ITS8 <-ES_ITS8[,time:=1:nrow(ES_ITS8)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS8$Biosimilarintro==1)
ES_ITS8 <- ES_ITS8[,time_after_inter:=ifelse(time<=b,0,time-b)]
data8 <- ES_ITS8
datanew8 <- data8
datanew8$Biosimilarintro=0
model8<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data8)
newey_adjust(model8)  

medicaid <- data.table(read_xlsx("MedicaidUSA_revised.xlsx"))
medicaid <- medicaid[year>=2010&year<=2020]
medicaid$Qtr <- as.character(medicaid$Qtr)
medicaid$year <- as.character(medicaid$year)
medicaid[,`conversion factor`:=1]#standardized to 2020 https://www.usinflationcalculator.com/
medicaid[year==2010,`conversion factor`:=1.19]
medicaid[year==2011,`conversion factor`:=1.15]
medicaid[year==2012,`conversion factor`:=1.13]
medicaid[year==2013,`conversion factor`:=1.11]
medicaid[year==2014,`conversion factor`:=1.09]
medicaid[year==2015,`conversion factor`:=1.09]
medicaid[year==2016,`conversion factor`:=1.08]
medicaid[year==2017,`conversion factor`:=1.06]
medicaid[year==2018,`conversion factor`:=1.03]
medicaid[year==2019,`conversion factor`:=1.01]
medicaid[year==2020,`conversion factor`:=1]
setnames(medicaid,"Qtr","quater")
medicaid$quater <- paste0(medicaid$quater," ")
medicaid <- merge(medicaid,ES_ITS8[,.(year,quater,qrt,time,Biosimilarintro,time_after_inter,pop)],by = c("year","quater"),all = T)
setnames(medicaid,"PricePerUnit","price_N")
medicaid$price_Y=0
medicaidnew <- medicaid
medicaidnew$Biosimilarintro=0
medicaid[,price_N:=price_N*`conversion factor`]
model_us <- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, medicaid)
newey_adjust(model_us)  
tidy(newey_adjust(model_us),conf.int = T,conf.level = 0.95)
# setnames(medicaid,"TotalUnits","overall")
# medicaid[,overall:=round(overall/pop,4)]
# model_us_con<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, medicaid)
# newey_adjust(model_us_con)  





ad1 <- res_plot3(model_us,1300,medicaid,medicaidnew,medicaid,"US - Medicaid")
ad2 <- res_plot3(model8,1300,data8,datanew8,ES_ITS8,"US - MIDAS")
figuread <- ggarrange(ad2,ad1)
tiff("Medicaid_price.tiff",width = 3600,height = 2000,res = 300,compression = "jpeg")
figuread 
dev.off()




ES_ITS9 <- ES_eli[nation=="UK"]
ES_ITS9 <-ES_ITS9[,time:=1:nrow(ES_ITS9)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS9$Biosimilarintro==1)
ES_ITS9 <- ES_ITS9[,time_after_inter:=ifelse(time<=b,0,time-b)]
data9 <- ES_ITS9
datanew9 <- data9
datanew9$Biosimilarintro=0
model9<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data9)
newey_adjust(model9)  

ES_ITS10 <- ES_eli[nation=="CANADA"]
ES_ITS10 <-ES_ITS10[,time:=1:nrow(ES_ITS10)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS10$Biosimilarintro==1)
ES_ITS10 <- ES_ITS10[,time_after_inter:=ifelse(time<=b,0,time-b)]
data10 <- ES_ITS10
datanew10 <- data10
datanew10$Biosimilarintro=0
model10<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data10)
newey_adjust(model10)  

ES_ITS11 <- ES_eli[nation=="AUSTRALIA"]
ES_ITS11 <-ES_ITS11[,time:=1:nrow(ES_ITS11)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS11$Biosimilarintro==1)
ES_ITS11 <- ES_ITS11[,time_after_inter:=ifelse(time<=b,0,time-b)]
data11 <- ES_ITS11
datanew11 <- data11
datanew11$Biosimilarintro=0
model11<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data11)
newey_adjust(model11)  

#Professor leung's comments

# 
# rsp <- data.table()
# for (i in 1:11) {
# r <- get(paste0("ES_ITS",i))[c(1,(get(paste0("ES_ITS",i))[,which(time_after_inter==1)]-1),get(paste0("ES_ITS",i))[,which(time_after_inter==1)],nrow(get(paste0("ES_ITS",i)))),.(Sales=quarter_sale_N+quarter_sale_Y,Consumption=overall,Area=nation,b=get(paste0("ES_ITS",i))[,which(time_after_inter==1)],a=44-get(paste0("ES_ITS",i))[,which(time_after_inter==1)])]
# rsp <- rbind(rsp,r)
# }
# rsp$ab <- rep(c("B1","B2","A1","A2"),11)
# rsp <- dcast(rsp,Area+a+b~ab,value.var = c("Sales","Consumption"))
# 
# rsp[,CQGR_sales:=((Sales_A2/Sales_A1)^(1/a)-1)*100]
# rsp[,CQGR_consumption:=((Consumption_A2/Consumption_A1)^(1/a)-1)*100]
# 
# write_xlsx(rsp,"Total_spending.xlsx")

# budget <- data.table()
# for (i in 1:11) {
#    t <- get(paste0("ES_ITS",i))[str_detect(time_after_inter,"^1$|^2$|^3$|^4$"),.(nation,time,Actual_sales=quarter_sale_N+quarter_sale_Y,quarter_consum_N,quarter_consum_Y)]
# t <- t[,Actual_consum:=quarter_consum_N+quarter_consum_Y]
# tt <- predict(get(paste0("model",i)),type="response",get(paste0("datanew",i)))
# t[,pre_price:=tt[t$time]]
# t[,wo_sales:=round(Actual_consum*pre_price)]
# t <- t[,lapply(.SD,sum),nation]
# t <- t[,.(nation,Actual_sales,Actual_consum,wo_sales)]
# t[,reduction:=round((wo_sales-Actual_sales)/wo_sales*100)] 
# t[,diff:=wo_sales-Actual_sales]
# budget <- rbind(budget,t)
# }
# write_xlsx(budget,"budget_saving.xlsx")



#### Tables and figures 1----
p1 <- res_plot(model1,1300,data1,datanew1,ES_ITS1,"HONG KONG")
p2 <- res_plot(model2,1300,data2,datanew2,ES_ITS2,"INDIA")
p3 <- res_plot(model3,1300,data3,datanew3,ES_ITS3,"JAPAN")
p4 <- res_plot(model4,1300,data4,datanew4,ES_ITS4,"KOREA")
p5 <- res_plot(model5,1300,data5,datanew5,ES_ITS5,"MALAYSIA")
p6 <- res_plot(model6,1300,data6,datanew6,ES_ITS6,"THAILAND")
p7 <- res_plot(model7,1300,data7,datanew7,ES_ITS7,"VIETNAM")
p8 <- res_plot(model8,1300,data8,datanew8,ES_ITS8,"US")
p9 <- res_plot(model9,1300,data9,datanew9,ES_ITS9,"UK")
p10 <- res_plot(model10,1300,data10,datanew10,ES_ITS10,"CANADA")
p11 <- res_plot(model11,1300,data11,datanew11,ES_ITS11,"AUSTRALIA")

figure <- ggarrange(p11+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    p10+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    p1+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    p2+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    p3+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    p4+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    p9+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
    p8+theme(axis.title.x = element_blank(),axis.title.y = element_blank()) ,
                    ncol = 4, nrow = 2)

figure <- annotate_figure(figure,
                
                bottom = text_grob("Quarter",
                                     size = 20),
                left = text_grob("INF price (2020 USD)", rot = 90,size = 20),
                )
    
    
tiff("Inf_price.tiff",width = 4800,height = 2000,res = 300,compression = "jpeg")
figure
dev.off()

r1 <- tidy(newey_adjust(model1),conf.int = T,conf.level = 0.95)
r2 <- tidy(newey_adjust(model2),conf.int = T,conf.level = 0.95)
r3 <- tidy(newey_adjust(model3),conf.int = T,conf.level = 0.95)
r4 <- tidy(newey_adjust(model4),conf.int = T,conf.level = 0.95)
r5 <- tidy(newey_adjust(model5),conf.int = T,conf.level = 0.95)
r6 <- tidy(newey_adjust(model6),conf.int = T,conf.level = 0.95)
r7 <- tidy(newey_adjust(model7),conf.int = T,conf.level = 0.95)
r8 <- tidy(newey_adjust(model8),conf.int = T,conf.level = 0.95)
r9 <- tidy(newey_adjust(model9),conf.int = T,conf.level = 0.95)
r10 <- tidy(newey_adjust(model10),conf.int = T,conf.level = 0.95)
r11 <- tidy(newey_adjust(model11),conf.int = T,conf.level = 0.95)

r <- list("Hong Kong"=r1,"India"=r2,"Japan"=r3,"Korea"=r4,
        "Malaysia"=r5,"Thailand"=r6,"Vietnam"=r7,"Us"=r8,"Uk"=r9,"Canada"=r10,"Australia"=r11)
write_xlsx(r,"inf_price.xlsx")

t1 <- tidy(t.dif(model1,data1,datanew1))
t2 <- tidy(t.dif(model2,data2,datanew2))
t3 <- tidy(t.dif(model3,data3,datanew3))
t4 <- tidy(t.dif(model4,data4,datanew4))
t5 <- tidy(t.dif(model5,data5,datanew5))
t6 <- tidy(t.dif(model6,data6,datanew6))
t7 <- tidy(t.dif(model7,data7,datanew7))
t8 <- tidy(t.dif(model8,data8,datanew8))
t9 <- tidy(t.dif(model9,data9,datanew9))
t10 <- tidy(t.dif(model10,data10,datanew10))
t11 <- tidy(t.dif(model11,data11,datanew11))



#Take 2020 as standard year, calculate the price reduction compared biosimilar versus bio-originator

ES_reduction1 <- ES_eli[year==2020][,.(ori=sum(price_N,na.rm = T)/4,sim=sum(price_Y,na.rm = T)/4),by=c("nation","year")][,reduction:=round((ori-sim)/ori*100)][!nation%in%c("TAIWAN","NEW","SINGAPORE")]
ES_reduction2 <- ES_eli[!quarter_consum_Y==0][,index:=min(qrt),by=c("nation")][index==qrt][!nation%in%c("TAIWAN","NEW","SINGAPORE")][,reduction_intro:=round((price_N-price_Y)/price_N*100)]
ES_reduction <- merge(ES_reduction1,ES_reduction2[,.(nation,reduction_intro,year,price_N,price_Y)],by = c("nation"),all.x = T)
ES_reduction[,`:=`(ori_re=round((price_N-ori)/price_N*100),sim_re=round((price_Y-sim)/price_Y*100))]
mean(ES_reduction[nation%in%c("INDIA","CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$reduction)
sd(ES_reduction[nation%in%c("INDIA","CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$reduction)
mean(ES_reduction[nation%in%c("INDIA","CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$reduction_intro)
sd(ES_reduction[nation%in%c("INDIA","CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$reduction_intro)

mean(ES_reduction[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$ori)
mean(ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$ori)
sd(ES_reduction[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$ori)
sd(ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$ori)
t.test(ES_reduction[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$ori,ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$ori)


mean(ES_reduction[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$sim)
mean(ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$sim)
sd(ES_reduction[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$sim)
sd(ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$sim)
t.test(ES_reduction[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$sim,ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$sim)

mean(ES_reduction$reduction)
sd(ES_reduction$reduction)
mean(ES_reduction$reduction_intro)
sd(ES_reduction$reduction_intro)

write_xlsx(ES_reduction,"Price reduction.xlsx")

# Consumption ITS----
ES_ITS1 <- ES_eli[nation=="HONG"]
ES_ITS1 <-ES_ITS1[,time:=1:nrow(ES_ITS1)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS1$Biosimilarintro==1)
ES_ITS1 <- ES_ITS1[,time_after_inter:=ifelse(time<=b,0,time-b)]
data1 <- ES_ITS1
datanew1 <- data1
datanew1$Biosimilarintro=0
model1 <- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data1)
newey_adjust(model1)  


ES_ITS2 <- ES_eli[nation=="INDIA"]
ES_ITS2 <-ES_ITS2[,time:=1:nrow(ES_ITS2)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS2$Biosimilarintro==1)
ES_ITS2 <- ES_ITS2[,time_after_inter:=ifelse(time<=b,0,time-b)]
data2 <- ES_ITS2
datanew2 <- data2
datanew2$Biosimilarintro=0
model2<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data2)
newey_adjust(model2)  


ES_ITS3 <- ES_eli[nation=="JAPAN"]
ES_ITS3 <-ES_ITS3[,time:=1:nrow(ES_ITS3)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS3$Biosimilarintro==1)
ES_ITS3 <- ES_ITS3[,time_after_inter:=ifelse(time<=b,0,time-b)]
data3 <- ES_ITS3
datanew3 <- data3
datanew3$Biosimilarintro=0
model3<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data3)
newey_adjust(model3)  


ES_ITS4 <- ES_eli[nation=="KOREA"]
ES_ITS4 <-ES_ITS4[,time:=1:nrow(ES_ITS4)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS4$Biosimilarintro==1)
ES_ITS4 <- ES_ITS4[,time_after_inter:=ifelse(time<=b,0,time-b)]
data4 <- ES_ITS4
datanew4 <- data4
datanew4$Biosimilarintro=0
model4<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data4)
newey_adjust(model4)  


ES_ITS5 <- ES_eli[nation=="MALAYSIA"]
ES_ITS5 <-ES_ITS5[,time:=1:nrow(ES_ITS5)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
ES_ITS5[year==2019&quarter_sale_Y==0,Biosimilarintro:=1]
b=44-sum(ES_ITS5$Biosimilarintro==1)
ES_ITS5 <- ES_ITS5[,time_after_inter:=ifelse(time<=b,0,time-b)]
data5 <- ES_ITS5
datanew5 <- data5
datanew5$Biosimilarintro=0
model5<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data5)
newey_adjust(model5)  


ES_ITS6 <- ES_eli[nation=="THAILAND"]
ES_ITS6 <-ES_ITS6[,time:=1:nrow(ES_ITS6)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS6$Biosimilarintro==1)
ES_ITS6 <- ES_ITS6[,time_after_inter:=ifelse(time<=b,0,time-b)]
data6 <- ES_ITS6
datanew6 <- data6
datanew6$Biosimilarintro=0
model6<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data6)
newey_adjust(model6)  


ES_ITS7 <- ES_eli[nation=="VIETNAM"]
ES_ITS7 <-ES_ITS7[,time:=1:nrow(ES_ITS7)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=34-sum(ES_ITS7$Biosimilarintro==1)
ES_ITS7 <- ES_ITS7[,time_after_inter:=ifelse(time<=b,0,time-b)]
data7 <- ES_ITS7
datanew7 <- data7
datanew7$Biosimilarintro=0
model7<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data7)
newey_adjust(model7)  


ES_ITS8 <- ES_eli[nation=="US"]
ES_ITS8 <-ES_ITS8[,time:=1:nrow(ES_ITS8)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS8$Biosimilarintro==1)
ES_ITS8 <- ES_ITS8[,time_after_inter:=ifelse(time<=b,0,time-b)]
data8 <- ES_ITS8
datanew8 <- data8
datanew8$Biosimilarintro=0
model8<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data8)
newey_adjust(model8)  

ES_ITS9 <- ES_eli[nation=="UK"]
ES_ITS9 <-ES_ITS9[,time:=1:nrow(ES_ITS9)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS9$Biosimilarintro==1)
ES_ITS9 <- ES_ITS9[,time_after_inter:=ifelse(time<=b,0,time-b)]
data9 <- ES_ITS9
datanew9 <- data9
datanew9$Biosimilarintro=0
model9<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data9)
newey_adjust(model9)  

ES_ITS10 <- ES_eli[nation=="CANADA"]
ES_ITS10 <-ES_ITS10[,time:=1:nrow(ES_ITS10)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS10$Biosimilarintro==1)
ES_ITS10 <- ES_ITS10[,time_after_inter:=ifelse(time<=b,0,time-b)]
data10 <- ES_ITS10
datanew10 <- data10
datanew10$Biosimilarintro=0
model10<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data10)
newey_adjust(model10)  

ES_ITS11 <- ES_eli[nation=="AUSTRALIA"]
ES_ITS11 <-ES_ITS11[,time:=1:nrow(ES_ITS11)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS11$Biosimilarintro==1)
ES_ITS11 <- ES_ITS11[,time_after_inter:=ifelse(time<=b,0,time-b)]
data11 <- ES_ITS11
datanew11 <- data11
datanew11$Biosimilarintro=0
model11<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data11)
newey_adjust(model11)  

#### Tables and figures 2----

p1 <- res_plot2(model1,2.5,data1,datanew1,ES_ITS1,"HONG KONG")
p2 <- res_plot2(model2,0.003,data2,datanew2,ES_ITS2,"INDIA")
p3 <- res_plot2(model3,2.5,data3,datanew3,ES_ITS3,"JAPAN")
p4 <- res_plot2(model4,2.5,data4,datanew4,ES_ITS4,"KOREA")
p5 <- res_plot2(model5,0.04,data5,datanew5,ES_ITS5,"MALAYSIA")
p6 <- res_plot2(model6,0.03,data6,datanew6,ES_ITS6,"THAILAND")
p7 <- res_plot2(model7,0.03,data7,datanew7,ES_ITS7,"VIETNAM")
p8 <- res_plot2(model8,10,data8,datanew8,ES_ITS8,"US")
p9 <- res_plot2(model9,10,data9,datanew9,ES_ITS9,"UK")
p10 <- res_plot2(model10,10,data10,datanew10,ES_ITS10,"CANADA")
p11 <- res_plot2(model11,10,data11,datanew11,ES_ITS11,"AUSTRALIA")

figure1 <- ggarrange(p11+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    p10+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    p9+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    p8+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
                    nrow = 1,ncol = 4)
figure2 <- ggarrange(p3+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    p4+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
    p1+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
                     nrow = 1,ncol = 4)
figure3 <- ggarrange(p2+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
                    nrow = 1,ncol = 4)

figure <- ggarrange(figure1,figure2,figure3,
    nrow = 3)

figure <- annotate_figure(figure,
                
                bottom = text_grob("Quarter",
                                     size = 20),
                left = text_grob("Overall INF consumption/1000 inhabitants", rot = 90,size = 20),
                )

tiff("Inf_sale.tiff",width = 4800,height = 3600,res = 300,compression = "jpeg")
figure
dev.off()

r1 <- tidy(newey_adjust(model1),conf.int = T,conf.level = 0.95)
r2 <- tidy(newey_adjust(model2),conf.int = T,conf.level = 0.95)
r3 <- tidy(newey_adjust(model3),conf.int = T,conf.level = 0.95)
r4 <- tidy(newey_adjust(model4),conf.int = T,conf.level = 0.95)
r5 <- tidy(newey_adjust(model5),conf.int = T,conf.level = 0.95)
r6 <- tidy(newey_adjust(model6),conf.int = T,conf.level = 0.95)
r7 <- tidy(newey_adjust(model7),conf.int = T,conf.level = 0.95)
r8 <- tidy(newey_adjust(model8),conf.int = T,conf.level = 0.95)
r9 <- tidy(newey_adjust(model9),conf.int = T,conf.level = 0.95)
r10 <- tidy(newey_adjust(model10),conf.int = T,conf.level = 0.95)
r11 <- tidy(newey_adjust(model11),conf.int = T,conf.level = 0.95)

r <- list("Hong Kong"=r1,"India"=r2,"Japan"=r3,"Korea"=r4,
        "Malaysia"=r5,"Thailand"=r6,"Vietnam"=r7,"Us"=r8,"Uk"=r9,"Canada"=r10,"Australia"=r11)
write_xlsx(r,"inf_sale.xlsx")

#consumption difference

ES_consum <- ES_eli[year==2020][!nation%in%c("TAIWAN","NEW","SINGAPORE")][,.(sum(overall),sum(quarter_consum_N_std),sum(quarter_consum_Y_std)),by=c("nation","year")]
# write_xlsx(ES_consum,"consumption.xlsx")
mean(ES_consum[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$V1)
sd(ES_consum[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$V1)
mean(ES_consum[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$V1)
sd(ES_consum[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$V1)
t.test(ES_consum[nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$V1,ES_reduction[!nation%in%c("CANADA","HONG","JAPAN","KOREA","UK","US","AUSTRALIA")]$V1)

# Sensitivity analysis of US----
ES_ITS11 <- ES_eli[nation=="US"]
ES_ITS11 <-ES_ITS11[,time:=1:nrow(ES_ITS11)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS11$Biosimilarintro==1)+4
ES_ITS11 <- ES_ITS11[,time_after_inter:=ifelse(time<=b,0,time-b)][time_after_inter==0,Biosimilarintro:=0]
data11 <- ES_ITS11
datanew11 <- data11
datanew11$Biosimilarintro=0
model11<- glm(price_N ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data11)
newey_adjust(model11) 
p11 <- res_plot(model11,1000,data11,datanew11,ES_ITS11,"US Infliximab originator price trend (One-year gap period)")


ES_ITS12 <- ES_eli[nation=="US"]
ES_ITS12 <-ES_ITS12[,time:=1:nrow(ES_ITS12)][,Biosimilarintro:=ifelse(quarter_sale_Y==0,0,1)]
b=44-sum(ES_ITS12$Biosimilarintro==1)+4
ES_ITS12 <- ES_ITS12[,time_after_inter:=ifelse(time<=b,0,time-b)][time_after_inter==0,Biosimilarintro:=0]
data12 <- ES_ITS12
datanew12 <- data12
datanew12$Biosimilarintro=0
model12<- glm(overall ~ Biosimilarintro + time+time_after_inter+factor(quater), family = gaussian, data12)
newey_adjust(model12) 
p12 <- res_plot2(model12,10,data12,datanew12,ES_ITS12,"US infliximab consumption trend (One-year gap period)")

figure.gap <- ggarrange(p11,p12,
                    ncol = 1, nrow = 2)
tiff("Inf_gap.tiff",width = 1600,height = 2000,res = 300,compression = "jpeg")
figure.gap
dev.off()


r11 <- tidy(newey_adjust(model11),conf.int = T,conf.level = 0.95)
r12 <- tidy(newey_adjust(model12),conf.int = T,conf.level = 0.95)

r.gap <- list("Us-gap-price"=r11,"Us-gap-sale"=r12)
write_xlsx(r.gap,"inf_gap.xlsx")



# Consumption distribution----
es_con_graph <- melt(ES_eli,id.vars = c("nation","year"),measure.vars = c("quarter_consum_N_std","quarter_consum_Y_std"),variable.name = "Agent")
es_con_graph <- es_con_graph[,.(value=sum(value)),by=c("nation","year","Agent")]
es_con_graph <- es_con_graph[Agent=="quarter_consum_N_std",Agent:="Originator"]
es_con_graph <- es_con_graph[Agent=="quarter_consum_Y_std",Agent:="Biosimilar"]
es_con_graph <- es_con_graph[nation=="HONG",nation:="HONG KONG"]
test <- ES_eli[,round(sum(quarter_consum_Y)/(sum(quarter_consum_Y)+sum(quarter_consum_N))*100,digits = 1),by=c("nation","year")]
test <- test[nation=="HONG",nation:="HONG KONG"]
es_con_graph <- merge(es_con_graph,test,all.x = T)
es_con_graph <- es_con_graph[V1==0,V1:=NA]
es_con_graph <- es_con_graph[Agent=="Biosimilar",V1:=NA]
es_con_graph[,Height:=sum(value),by=c("nation","year")]
es_con_graph$V1 <- as.character(es_con_graph$V1)
es_con_graph[!is.na(V1),V1:=paste0(V1,"%")]
es_con_graph <- es_con_graph[nation%in%c("AUSTRALIA","CANADA","HONG KONG","INDIA","JAPAN","KOREA","UK","US")]

getplot <- function(w,t,d){
    ggplot(data=es_con_graph[nation==w],aes(y=value,x=year))+geom_bar(aes(fill=Agent),stat="identity")+
    geom_text(aes(label=V1),size=3,y=es_con_graph[nation==w]$Height,hjust=-0.2)+ggtitle(w)+
    theme(axis.text.x = element_text(size=rel(1.2)))+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
     theme(axis.text.x=element_text(angle=90,hjust=0.1))+scale_fill_discrete(name="Agent") +
        coord_flip()+scale_y_continuous(expand = expansion(mult =c(0, 0.3)),limits = c(t,d))+
        theme(legend.title = element_text( size=10, face="bold"),legend.text = element_text( size=10, face="bold"))
}
z1 <- getplot("AUSTRALIA",t=0,d=40)
z2 <- getplot("CANADA",t=0,d=40)
z3 <- getplot("UK",t=0,d=40)
z4 <- getplot("US",t=0,d=40)
z5 <- getplot("JAPAN",t=0,d=10)
z6 <- getplot("KOREA",t=0,d=10)
z7 <- getplot("HONG KONG",t=0,d=10)
z8 <- getplot("INDIA",t=0,d=0.012)

figurez1 <- ggarrange(z1+ theme(legend.position="none"),z2+ theme(legend.position="none"),z3+ theme(legend.position="none"),z4+ theme(legend.position="none"),
                    nrow = 4,ncol = 1)
figurez2 <- ggarrange(z5+ theme(legend.position="none"),z6+ theme(legend.position="none"),z7+ theme(legend.position="none"),
                     nrow = 4,ncol = 1)
figurez3 <- ggarrange(z8+ theme(legend.position="none"),
                    nrow = 4,ncol = 1)

library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(z1)

figurez <- ggarrange(figurez1,figurez2,figurez3,ncol = 3,legend.grob = legend,legend = "right")

figurezz <- annotate_figure(figurez,
                
                bottom = text_grob("Overall INF consumption/1000 inhabitants",
                                     size = 20,face = "bold"),
                left = text_grob("Year", rot = 90,size = 20),
                )

tiff("Inf_uptake.tiff",width = 4800,height = 3600,res = 300,compression = "jpeg")
figurezz
dev.off()
