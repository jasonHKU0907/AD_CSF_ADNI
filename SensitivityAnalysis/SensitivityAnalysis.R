#run differential analysis
library(tidyr)
library(dplyr)
m.data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")

m.data$PTGENDER=as.factor(m.data$PTGENDER)
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DXBLYES","+AGE+PTGENDER+PTEDUCAT+APOEcarrier"))
    Proname=x
    glm1=lm(FML,data=df)
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("Dementia","MCI")) {
  df=filter(m.data,m.data$DXBLYES%in%c(j,"CN"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                       
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BICN.csv")
  write.csv(result,filename)
  
}



for (j in c("Dementia")) {
  df=filter(m.data,m.data$DXBLYES%in%c(j,"MCI"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                       
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BIMCI.csv")
  write.csv(result,filename)
  
}



## 
library(tidyr)
library(dplyr)
m.data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")

m.data$PTGENDER=as.factor(m.data$PTGENDER)
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
m.data$DXBLYES[which(m.data$DX_AT_status=="CN A+T-")]="DCN A+T-"   
m.data$DXBLYES[which(m.data$DX_AT_status=="CN A+T+")]="ECN A+T+"
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DX_AT_status","+AGE+PTGENDER+PTEDUCAT+APOEcarrier"))
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("CN A+T-","CN A+T+","MCI A+T+","ZDementia A+T+")) {
  df=filter(m.data,m.data$DX_AT_status%in%c(j,"CN A-T-"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                       
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BICN A-T-.csv")
  write.csv(result,filename)
  
}




for (j in c("CN A+T+")) {
  df=filter(m.data,m.data$DX_AT_status%in%c(j,"CN A+T-"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                       
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BICN A+T-.csv")
  write.csv(result,filename)
  
}


for (j in c("MCI A+T+")) {
  df=filter(m.data,m.data$DX_AT_status%in%c(j,"CN A+T+"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                       
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BICN A+T+.csv")
  write.csv(result,filename)
  
}




for (j in c("ZDementia A+T+")) {
  df=filter(m.data,m.data$DX_AT_status%in%c(j,"MCI A+T+"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                       
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BIMCI A+T+.csv")
  write.csv(result,filename)
  
}



## subgroup
library(tidyr)
library(dplyr)
data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")
m.data<-data.frame(dplyr::filter(data,APOEcarrier=="NO"))  #APOEcarrier=="YES"
m.data$PTGENDER=as.factor(m.data$PTGENDER)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DXBLYES","+AGE+PTGENDER+PTEDUCAT")) 
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Analytes","Beta","P","t")
    return(Uni_glm) 
  } 

for (j in c("Dementia")) {
  df=filter(m.data,m.data$DXBLYES%in%c(j,"CN"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P, 
    method ="bonferroni"            
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("APOENO_",j,"BICN.csv")
  write.csv(result,filename)
  
}





## subgroup
library(tidyr)
library(dplyr)
data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")
m.data<-data.frame(dplyr::filter(data,PTGENDER=="Female"))  #PTGENDER=="Male"
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DXBLYES","+AGE+PTEDUCAT+APOEcarrier")) 
    Proname=x
    glm1=lm(FML,data=df)
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("Dementia")) {
  df=filter(m.data,m.data$DXBLYES%in%c(j,"CN"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P, 
    method ="bonferroni"              
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("Female_",j,"BICN.csv")
  write.csv(result,filename)
  
}



##subgroup
library(tidyr)
library(dplyr)
data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")
m.data<-data.frame(dplyr::filter(data,AGE>65))  #AGE<=65
m.data$PTGENDER=as.factor(m.data$PTGENDER)
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DXBLYES","+AGE+PTGENDER+PTEDUCAT+APOEcarrier")) 
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("Dementia")) {
  df=filter(m.data,m.data$DXBLYES%in%c(j,"CN"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"      
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("over65years_",j,"BICN.csv")
  write.csv(result,filename)
  
}





## subgroup
library(tidyr)
library(dplyr)
data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")
m.data<-data.frame(dplyr::filter(data,APOEcarrier=="NO"))  #APOEcarrier=="YES"
m.data$PTGENDER=as.factor(m.data$PTGENDER)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DX_AT_status1","+AGE+PTGENDER+PTEDUCAT")) 
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("ZA+T+")) {
  df=filter(m.data,m.data$DX_AT_status1%in%c(j,"CN A-T-"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P, 
    method ="bonferroni" 
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("APOENO_",j,"BICN A-T-.csv")
  write.csv(result,filename)
  
}



## subgroup
library(tidyr)
library(dplyr)
data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")
m.data<-data.frame(dplyr::filter(data,PTGENDER=="Female"))  #PTGENDER=="Male"
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DX_AT_status1","+AGE+PTEDUCAT+APOEcarrier"))
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("ZA+T+")) {
  df=filter(m.data,m.data$DX_AT_status1%in%c(j,"CN A-T-"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"           
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("Male_",j,"BICN A-T-.csv")
  write.csv(result,filename)
  
}



## subgroup
library(tidyr)
library(dplyr)
data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")
m.data<-data.frame(dplyr::filter(data,AGE>65))   #AGE<=65
m.data$PTGENDER=as.factor(m.data$PTGENDER)
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
variable.names=colnames(m.data)[9:7016] 
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "DX_AT_status1","+AGE+PTGENDER+PTEDUCAT+APOEcarrier")) 
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("ZA+T+")) {
  df=filter(m.data,m.data$DX_AT_status1%in%c(j,"CN A-T-"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"            
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("over65years_",j,"BICN A-T-.csv")
  write.csv(result,filename)
  
}




##autopsy
library(tidyr)
library(dplyr)
m.data=read.csv("C:/Users/yucuo/Desktop/CSF_SOMAscan_covariates_CSF707.csv")

m.data$PTGENDER=as.factor(m.data$PTGENDER)
m.data$APOEcarrier=as.factor(m.data$APOEcarrier)
variable.names=colnames(m.data)[9:7016]
m.data[9:7016]=scale(log(m.data[9:7016]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "autopsy_DX","+AGE+PTGENDER+PTEDUCAT+APOEcarrier"))
    Proname=x
    glm1=lm(FML,data=df)
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("ZAD")) {
  df=filter(m.data,m.data$autopsy_DX%in%c(j,"Non_AD"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("AGE","PTGENDERMale","PTEDUCAT","APOEcarrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  # bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P, 
    method ="bonferroni"
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BINon_AD.csv")
  write.csv(result,filename)
  
}



## PPMI
library(tidyr)
library(dplyr)
m.data=read.csv("C:/Users/yucuo/Desktop/PPMI_CSF_proteomics_HC174.csv")

m.data$SEX=as.factor(m.data$SEX)
m.data$APOE4carrier=as.factor(m.data$APOE4carrier)
variable.names=colnames(m.data)[154:4225] 
m.data[154:4225]=scale(log(m.data[154:4225]))

Uni_glm_model=
  function(x){
    FML=as.formula(paste0(x,"~", "diagnosis2","+age+SEX+education+APOE4carrier")) 
    Proname=x
    glm1=lm(FML,data=df) 
    LSUM=summary(glm1)
    beta=LSUM$coefficients[-1,1]
    p=LSUM$coefficients[-1,4]
    t=LSUM$coefficients[-1,3]
    Uni_glm=cbind(Proname,beta,p,t)
    Uni_glm=as.data.frame(Uni_glm)
    dimnames(Uni_glm)[[2]]=c("Protein","AD_Estimate","P","t")
    return(Uni_glm) 
  } 

for (j in c("Control_A+T+")) {
  df=filter(m.data,m.data$diagnosis2%in%c(j,"Control_A-T-"))
  Uni_glm_result=vector(mode="list",length=length(variable.names))
  Uni_glm_result=lapply(variable.names,Uni_glm_model) 
  result=data.frame()       
  Uni_glm2=vector(mode="list",length=length(variable.names))
  for (i in 1:length(variable.names)) {
    Uni_glm2[[i]]=Uni_glm_result[[i]][-which(row.names(Uni_glm_result[[i]])%in%c("age","SEXMale","education","APOE4carrierYES")),] 
    result=rbind(result,Uni_glm2[[i]]) 
  }
  #bonferroni correction
  result$P=as.numeric(result$P)
  result=result[order(result$P),]
  result$P_BH_adjust=p.adjust(
    result$P,  
    method ="bonferroni"                
  )
  setwd("C:/Users/yucuo/Desktop/")
  filename=paste0("LR_",j,"BIControl_A-T-.csv")
  write.csv(result,filename)
  
}








#run predictive performance
rm(list = ls())
setwd("E:/My works/Proteomics_AD/3_Predictive/")

#----- cox data -------
data <- read.csv('E:/My works/Proteomics_AD/DATA/CSF_SOMAscan_covariates_CSF707.csv')
Pt <- data[,c(1,3:7016)] # ADNI proteomics data
names(Pt[,c(1:30)])
data_progression <- read.csv('E:/My works/Proteomics_AD/3_Predictive/DATA/DATA_progression_status.csv')
names(data_progression)
data_cox <- merge(data_progression, Pt, by = 'id')
data_cox <- merge(data[,c('id','AGE','PTGENDER','PTEDUCAT','APOE4','APOEcarrier',
                          "CSF_ABETA40","CSF_ABETA42","CSF_TTAU","CSF_PTAU181",
                          "CSF_ABETA42_status1098","CSF_PTAU181_status26.64","AT_status","DXBLYES",
                          "DX_AT_status","DX_AT_status1")], data_cox, by = 'id')
names(data_cox[,c(1:40)])
summary(data_cox[,c(1:40)])
DXs <- colnames(data_cox)[grep('DX',colnames(data_cox))]
STA <- colnames(data_cox)[grep('status',colnames(data_cox))]
ls <- c(DXs,STA,'PTGENDER','APOE4','APOEcarrier')
for (i in ls) {
  data_cox[[i]]=as.factor(data_cox[[i]])
}
summary(data_cox[,c(1:40)])

# # Pt_panel
# Pt_panel  <- read.csv('E:/My works/Proteomics_AD/DATA/Protein_Panel.csv')
# data_cox <- merge(Pt_panel, data_cox, by = 'RID', all.x = T)
# Pt_ls <- names(Pt_panel[,c(2,3)])

# selected proteins
Pt_select <- read.csv('E:/My works/Proteomics_AD/12_Pt_select/bonferroni_significant_between_group_proteins_272.csv')
Pt_select <- Pt_select[!is.na(Pt_select$MARK),]
summary(Pt_select)
Pt_ls <- Pt_select$Analytes 

# median
for (i in Pt_ls) {
  Q2 <- quantile(scale(log(data_cox[[i]])), probs = seq(0, 1, 1/2), na.rm = T)
  data_cox[[paste0(i,'_C2')]] <- ifelse(scale(log(data_cox[[i]])) < Q2[2], 'Low','High')
  data_cox[[paste0(i,'_C2')]] <- factor(data_cox[[paste0(i,'_C2')]], levels = c('High','Low'),
                                        labels = c('High','Low'))
} 
Pt_ls_C2 <- c()
for (i in Pt_ls) {
  a <- paste0(i,'_C2')
  Pt_ls_C2 <- c(Pt_ls_C2,a)
} 


#----- Group comparison(COX) -------
library(survival)
library(plyr)
library(dplyr)

### C2 #####
# (1,2,4)DX_progression
Progression_ls <- c('DX_progression_MCI')
UniCox_model<-function(x){
  FML<-as.formula(paste0('basurv~',paste(paste0('relevel(',x,',','ref=','"Low"',')'),'AGE',"PTGENDER","APOEcarrier","PTEDUCAT",sep="+")))
  cox<-coxph(FML,data = DATA)
  sum<-summary(cox)
  beta<-sum$coefficients[1,1]
  HR<-sum$coefficients[1,2]
  CI<-paste0(round(sum$conf.int[1,3:4],2),collapse = "-")
  LCI<-sum$conf.int[1,3]
  UCI<-sum$conf.int[1,4]
  HR_95CI<-paste(round(HR,2)," (",CI,")", sep = "")
  PValue<-sum$coefficients[1,5]
  N=sum$n
  Nevent=sum$nevent
  Nevent_N=paste(Nevent,"/",N, sep = "")
  median_FL <- median(DATA$Years_D_progression)
  range <- paste0(round(range(DATA$Years_D_progression)[1],2),'-',round(range(DATA$Years_D_progression)[2],2))
  IQR <- paste0(round(quantile(DATA$Years_D_progression, 0.25),2),'-',round(quantile(DATA$Years_D_progression, 0.75),2))
  Unicox<-data.frame(UniCox<-data.frame('Analytes'=x,
                                        'Progression'=i,
                                        'Ref'='Low',
                                        'median_FL'=median_FL,
                                        'range_FL'=range,
                                        'IQR_FL'=IQR,
                                        'HR'=HR,
                                        'LCI'=LCI,
                                        'UCI'=UCI,
                                        'Nevent_N'=Nevent_N,
                                        'HR_95CI'=HR_95CI,
                                        'PValue'=PValue))
  return(Unicox)
}
Results <- data.frame()
for (i in Progression_ls) {
  for (n in Pt_ls_C2) {
    DATA <- data_cox
    DATA <- DATA[!is.na(DATA[[i]]),]
    DATA <- DATA[!is.na(DATA[[n]]),]
    basurv<-Surv(DATA$Years_D_progression, DATA[[i]]==1)
    DATA$basurv <- with(DATA, basurv)
    R <-UniCox_model(n)
    Results <- rbind(Results, R)
  }
}
Results01 <- Results


# (5)CN_CDR_progression
Progression_ls <- c('CN_CDR_progression')
UniCox_model<-function(x){
  FML<-as.formula(paste0('basurv~',paste(paste0('relevel(',x,',','ref=','"Low"',')'),'AGE',"PTGENDER","APOEcarrier","PTEDUCAT",sep="+")))
  cox<-coxph(FML,data = DATA)
  sum<-summary(cox)
  beta<-sum$coefficients[1,1]
  HR<-sum$coefficients[1,2]
  CI<-paste0(round(sum$conf.int[1,3:4],2),collapse = "-")
  LCI<-sum$conf.int[1,3]
  UCI<-sum$conf.int[1,4]
  HR_95CI<-paste(round(HR,2)," (",CI,")", sep = "")
  PValue<-sum$coefficients[1,5]
  N=sum$n
  Nevent=sum$nevent
  Nevent_N=paste(Nevent,"/",N, sep = "")
  median_FL <- median(DATA$Years_CN_CDR_progression)
  range <- paste0(round(range(DATA$Years_CN_CDR_progression)[1],2),'-',round(range(DATA$Years_CN_CDR_progression)[2],2))
  IQR <- paste0(round(quantile(DATA$Years_CN_CDR_progression, 0.25),2),'-',round(quantile(DATA$Years_CN_CDR_progression, 0.75),2))
  Unicox<-data.frame(UniCox<-data.frame('Analytes'=x,
                                        'Progression'=i,
                                        'Ref'='Low',
                                        'median_FL'=median_FL,
                                        'range_FL'=range,
                                        'IQR_FL'=IQR,
                                        'HR'=HR,
                                        'LCI'=LCI,
                                        'UCI'=UCI,
                                        'Nevent_N'=Nevent_N,
                                        'HR_95CI'=HR_95CI,
                                        'PValue'=PValue))
  return(Unicox)
}
Results <- data.frame()
for (i in Progression_ls) {
  for (n in Pt_ls_C2) {
    DATA <- data_cox
    DATA <- DATA[!is.na(DATA[[i]]),]
    DATA <- DATA[!is.na(DATA[[n]]),]
    basurv<-Surv(DATA$Years_CN_CDR_progression, DATA[[i]]==1)
    DATA$basurv <- with(DATA, basurv)
    R <-UniCox_model(n)
    Results <- rbind(Results, R)
  }
}
Results03 <- Results

Results <- rbind(Results01, Results02, Results03)
Results$Analytes <- sub('_C2','',Results$Analytes)
Results <- merge(Pt_select, Results, by = 'Analytes', all.y = T, no.dups = F)
write.csv(Results, 'E:/My works/Proteomics_AD/3_Predictive/Results/Results_Pt_progression_Cox_C2refLow.csv', row.names = F)


#----- KM plots --------
library("survival")
library("survminer")

### C2 #####
Progression_ls <- c('DX_progression_CN','DX_progression_MCI')
for (i in Progression_ls) {
  for (n in Pt_ls_C2) {
    DATA <- data_cox
    DATA <- DATA[!is.na(DATA[[i]]),]
    DATA <- DATA[!is.na(DATA[[n]]),]
    pdf(paste0('E:/My works/Proteomics_AD/3_Predictive/Results/KM/C2/','KM_',n,'_',i,'.pdf'))
    fit <- survfit(Surv(DATA$Years_D_progression, DATA[[i]]==1) ~ DATA[[n]], data=DATA)
    p1 <- ggsurvplot(fit,
                     title  = i,  
                     ylab = "Cumulative survival", 
                     xlab = "Follow-up time (years)", 
                     risk.table = "abs_pct",
                     # risk.table.col = "black",
                     # tables.y.text = F,
                     risk.table.y.text = T,  
                     risk.table.height = 0.25,  
                     fontsize = 4,  # 
                     risk.table.title = "Number at risk (%)",
                     xlim = c(-0.1, 11),
                     # ylim = c(0.25,1),
                     break.x.by = 2,
                     break.y.by = 0.25,
                     conf.int = T, 
                     conf.int.alpha = 0.2,
                     conf.int.style = 'ribbon',
                     pval = TRUE, 
                     pval.coord = c(2, 0.1),
                     pval.size = 5,
                     pval.method=TRUE,
                     pval.method.size=5,
                     pval.method.coord=c(0.2,0.1),
                     # surv.median.line = "hv",    # Add medians survival
                     legend = c(0.2, 0.35),
                     legend.title = " ",
                     legend.labs = c("High", "Low"),
                     font.legend = c(14, "plain"),
                     font.main = c(14, "plain"), 
                     font.x = c(14,  "plain"),   
                     font.y = c(14,  "plain"),   
                     font.tickslab = c(12, "plain"), 
                     size = 1,  
                     linetype = 1,
                     # ggtheme = theme_bw(),       # Change ggplot2 theme
                     palette = c("#BC3C29FF", "#0072B5FF"))    
    print(p1)
    dev.off()
  }
}











#run associations with AD endophenotypes
rm(list = ls())
setwd("E:/My works/Proteomics_AD/4_ADv_association/")

#----- (Cross-sectional) LR ------
#### data preparation #####
data0 <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA4_ADv_longitudinal.csv')
data0 <- data0[!is.na(data0$Pt_CompVis),]
data0 <- data0[data0$VISCODE =='bl',]
data0$CSF_ABETA42_40 <- data0$CSF_ABETA42/data0$CSF_ABETA40
data0$CSF_PTAU181_ABETA40 <- data0$CSF_PTAU181/data0$CSF_ABETA40
names(data0)
data0$DX_AT_status <- paste(data0$DX, data0$CSF_AT_status)

# data processing
data <- read.csv('E:/My works/Proteomics_AD/DATA/CSF_SOMAscan_covariates_CSF707.csv')
Pt <- data[,c(1,9:7016)] 
names(Pt[,c(1:30)])
data_progression <- read.csv('E:/My works/Proteomics_AD/3_Predictive/DATA/DATA_progression_status.csv')
names(data_progression)
data_all <- merge(data_progression, Pt, by = 'id')
data_all <- merge(data0, data_all, by = 'id')

DXs <- colnames(data_all)[grep('DX',colnames(data_all))]
STA <- colnames(data_all)[grep('status',colnames(data_all))]
ls <- c(DXs,STA,'PTGENDER','APOE4','APOEcarrier')
for (i in ls) {
  data_all[[i]]=as.factor(data_all[[i]])
}
summary(data_all[,c(1:40)])

# selected proteins
Pt_select <- read.csv('E:/My works/Proteomics_AD/12_Pt_select/bonferroni_significant_between_group_proteins_272.csv')
Pt_select <- Pt_select[!is.na(Pt_select$MARK),]
summary(Pt_select)
Pt_ls <- Pt_select$Analytes #protein list


#### Linear regression #####

library(car)
lmOutput <- function(fit){
  sum <-summary(fit)
  beta <-sum$coefficients[2,1]
  se <-sum$coefficients[2,2]
  tval <-sum$coefficients[2,3]
  PValue <-sum$coefficients[2,4]
  Biomarker.y <- i
  Analytes.x <- n
  N <- nrow(DATA)
  Output <- cbind(Analytes.x, Biomarker.y, beta, PValue, se, tval, N)
  Output
}
ADV_ls <- c("CSF_ABETA40","CSF_ABETA42","CSF_TTAU","CSF_PTAU181",'CSFNFL','CSF_PGRN_MSDCRT','CSF_sTREM2_MSDCRT',
            "PLASMAPTAU181","PLASMA_NFL","PLASMA_AB42",
            "PET_AV45_SUVR","PET_FDG_SUVR")
Cog_ls <- c("mPACCdigit","mPACCtrailsB","ADNI_MEM","ADNI_EF")
MRI_ls <- c("Hippocampus","Temporal_Meta")
{
  # (1) Total
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'Total'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_total <- lm_result
  write.csv(lm_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_total.csv', row.names = F)
  
  
  # (2) stratiried by CSF Ab42
  # A-
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$CSF_A_status == 'A-',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$CSF_A_status == 'A-',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$CSF_A_status == 'A-',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'CSF A-'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_An <- lm_result
  write.csv(lm_result_An, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_An.csv', row.names = F)
  
  # A+
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$CSF_A_status == 'A+',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$CSF_A_status == 'A+',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$CSF_A_status == 'A+',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'CSF A+'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_Ap <- lm_result
  write.csv(lm_result_Ap, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_Ap.csv', row.names = F)
  
  lm_result_A <- rbind(lm_result_Ap, lm_result_An)
  lm_result_A$P_Bon <- p.adjust(lm_result_A$PValue, method = 'bonferroni')
  lm_result_A$P_FDR <- p.adjust(lm_result_A$PValue, method = 'fdr')
  write.csv(lm_result_A, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_A.csv', row.names = F)
  
  
  # (3) CN/MCI/AD
  # CN
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'CN',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'CN',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'CN',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'CN'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_CN <- lm_result
  write.csv(lm_result_CN, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_CN.csv', row.names = F)
  
  # MCI
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'MCI',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'MCI',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'MCI',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'MCI'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_MCI <- lm_result
  write.csv(lm_result_MCI, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_MCI.csv', row.names = F)
  
  # ZDementia 
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'ZDementia',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'ZDementia',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DXBLYES == 'ZDementia',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'ZDementia'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_ZDementia <- lm_result
  write.csv(lm_result_ZDementia, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_ZDementia.csv', row.names = F)
  
  lm_result_S3 <- rbind(lm_result_CN, lm_result_MCI, lm_result_ZDementia)
  lm_result_S3$P_Bon <- p.adjust(lm_result_S3$PValue, method = 'bonferroni')
  lm_result_S3$P_FDR <- p.adjust(lm_result_S3$PValue, method = 'fdr')
  write.csv(lm_result_S3, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_S3.csv', row.names = F)
  
  
  # (4)AD continuum(CN A-T-,CN A+T-,CN A+T+,MCI A+T+, Dementia A+T+)
  lm_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DX_AT_status == 'CN A-T-' | data_all$DX_AT_status == 'CN A+T-' | data_all$DX_AT_status == 'CN A+T+' | data_all$DX_AT_status == 'MCI A+T+' | data_all$DX_AT_status == 'Dementia A+T+',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DX_AT_status == 'CN A-T-' | data_all$DX_AT_status == 'CN A+T-' | data_all$DX_AT_status == 'CN A+T+' | data_all$DX_AT_status == 'MCI A+T+' | data_all$DX_AT_status == 'Dementia A+T+',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls){
      DATA <- data_all[data_all$DX_AT_status == 'CN A-T-' | data_all$DX_AT_status == 'CN A+T-' | data_all$DX_AT_status == 'CN A+T+' | data_all$DX_AT_status == 'MCI A+T+' | data_all$DX_AT_status == 'Dementia A+T+',]
      DATA<-DATA[complete.cases(DATA[[n]]),]
      DATA<-DATA[complete.cases(DATA[[i]]),]
      fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
      Row <- lmOutput(fit)
      lm_result <- rbind(lm_result, Row)
    }
  }
  lm_result$Sample <- 'AD continuum'
  lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
  lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
  lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
  lm_result_ADC <- lm_result
  write.csv(lm_result_ADC, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional_ADC.csv', row.names = F)
  
}

# merge
lm_result <- rbind(lm_result_total, lm_result_A, lm_result_S3, lm_result_ADC)
write.csv(lm_result, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_crosssectional.csv', row.names = F)


#### MRI #####
data <- data_all
MRI  <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA_MRI_longitudinal.csv')
data <- merge(data, MRI[,c(1,65:151)], by = 'id', all.x = T)
CTX <- colnames(data)[grep('CV',colnames(data))]
CTX <- CTX[-c(1,2)]
SCT <- colnames(data)[grep('SV',colnames(data))]
SCT <- SCT[-c(1,2)]
# MRI_ls <- c(CTX,SCT)
imaging_match <- read.csv('E:/My works/Glymphatic/2_AD_variable_cross_sectional/Brain_region_formatting_MRI.csv')

library(car)
lmOutput <- function(fit){
  sum <-summary(fit)
  beta <-sum$coefficients[2,1]
  se <-sum$coefficients[2,2]
  tval <-sum$coefficients[2,3]
  PValue <-sum$coefficients[2,4]
  Biomarker.y <- i
  Analytes.x <- n
  N <- nrow(DATA)
  Output <- cbind(Analytes.x, Biomarker.y, beta, PValue, se, tval, N)
  Output
}
# (1) Total
lm_result <- data.frame()
for (i in CTX) {
  for (n in Pt_ls){
    DATA <- data
    DATA<-DATA[complete.cases(DATA[[n]]),]
    DATA<-DATA[complete.cases(DATA[[i]]),]
    fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
    Row <- lmOutput(fit)
    lm_result <- rbind(lm_result, Row)
  }
}
for (i in SCT) {
  for (n in Pt_ls){
    DATA <- data
    DATA<-DATA[complete.cases(DATA[[n]]),]
    DATA<-DATA[complete.cases(DATA[[i]]),]
    fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV, data = DATA)
    Row <- lmOutput(fit)
    lm_result <- rbind(lm_result, Row)
  }
}
lm_result$Sample <- 'Total'
lm_result <- merge(lm_result, imaging_match, by.x = 'Biomarker.y', by.y = 'FLDNAME', all.x = T, no.dups = F)
lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
lm_result <- lm_result[!lm_result$Structure=='',]
lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
lm_result_total <- lm_result
write.csv(lm_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_MRI_crosssectional_total.csv', row.names = F)


#### AV45 #####
data <- data_all
AV45 <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA_AV45_longitudinal.csv')
data <- merge(data, AV45[,c(1,45:213)], by = 'id', all.x = T)
ctxl <- colnames(data)[grep('CTX_LH',colnames(data))]
ctxr <- colnames(data)[grep('CTX_RH',colnames(data))]
RT <- colnames(data)[grep('RIGHT_',colnames(data))]
LT <- colnames(data)[grep('LEFT_',colnames(data))]
SCT <- c(RT,LT)
CTX <- c(ctxl, ctxr)
imaging_match <- read.csv('E:/My works/Glymphatic/5_ADv_progression/Brain_region_formatting.csv')

library(car)
lmOutput <- function(fit){
  sum <-summary(fit)
  beta <-sum$coefficients[2,1]
  se <-sum$coefficients[2,2]
  tval <-sum$coefficients[2,3]
  PValue <-sum$coefficients[2,4]
  Biomarker.y <- i
  Analytes.x <- n
  N <- nrow(DATA)
  Output <- cbind(Analytes.x, Biomarker.y, beta, PValue, se, tval, N)
  Output
}
# (1) Total
lm_result <- data.frame()
for (i in CTX) {
  for (n in Pt_ls){
    DATA <- data
    DATA<-DATA[complete.cases(DATA[[n]]),]
    DATA<-DATA[complete.cases(DATA[[i]]),]
    fit <-lm(scale(log(DATA[[i]])) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
    Row <- lmOutput(fit)
    lm_result <- rbind(lm_result, Row)
  }
}
for (i in SCT) {
  for (n in Pt_ls){
    DATA <- data
    DATA<-DATA[complete.cases(DATA[[n]]),]
    DATA<-DATA[complete.cases(DATA[[i]]),]
    fit <-lm(scale(DATA[[i]]) ~ scale(log(DATA[[n]])) + AGE + PTGENDER + PTEDUCAT + APOEcarrier, data = DATA)
    Row <- lmOutput(fit)
    lm_result <- rbind(lm_result, Row)
  }
}
lm_result$Sample <- 'Total'
lm_result <- merge(imaging_match, lm_result, by.y = 'Biomarker.y', by.x = 'Imaging', all.x = T, no.dups = F)
lm_result <- merge(Pt_select, lm_result, by.x = 'Analytes', by.y = 'Analytes.x', all.y = T, no.dups = F)
lm_result <- lm_result[!lm_result$Structure=='',]
lm_result$P_Bon <- p.adjust(lm_result$PValue, method = 'bonferroni')
lm_result$P_FDR <- p.adjust(lm_result$PValue, method = 'fdr')
lm_result_total <- lm_result
write.csv(lm_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_AV45_crosssectional_total.csv', row.names = F)


#----- (Longitudinal) LME ------
#### data preparation #####
data0 <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA4_ADv_longitudinal.csv')
data0 <- data0[!is.na(data0$Pt_CompVis),]
data0$CSF_ABETA42_40 <- data0$CSF_ABETA42/data0$CSF_ABETA40
data0$CSF_PTAU181_ABETA40 <- data0$CSF_PTAU181/data0$CSF_ABETA40
names(data0)
data0$DX_AT_status <- paste(data0$DX, data0$DX_AT_status)

data <- read.csv('E:/My works/Proteomics_AD/DATA/CSF_SOMAscan_covariates_CSF707.csv')
Pt <- data[,c(1,9:7016)] #proteomics data
names(Pt[,c(1:30)])
data_progression <- read.csv('E:/My works/Proteomics_AD/3_Predictive/DATA/DATA_progression_status.csv')
names(data_progression)
data_all <- merge(data_progression, Pt, by = 'id')
data_all <- merge(data0, subset(data_all, select = -id), by = 'RID', no.dups = F, all.x = T)

DXs <- colnames(data_all)[grep('DX',colnames(data_all))]
STA <- colnames(data_all)[grep('status',colnames(data_all))]
ls <- c(DXs,STA,'PTGENDER','APOE4','APOEcarrier')
for (i in ls) {
  data_all[[i]]=as.factor(data_all[[i]])
}
summary(data_all[,c(1:40)])

# selected proteins
Pt_select <- read.csv('E:/My works/Proteomics_AD/12_Pt_select/bonferroni_significant_between_group_proteins_272.csv')
Pt_select <- Pt_select[!is.na(Pt_select$MARK),]
summary(Pt_select)
Pt_ls <- Pt_select$Analytes #protein list


#----- LME analysis #####
library(lme4)
library(lmerTest)

ADV_ls <- c("CSF_ABETA40","CSF_ABETA42","CSF_TTAU","CSF_PTAU181",'CSF_PGRN_MSDCRT','CSF_sTREM2_MSDCRT',
            # 'CSFNFL',# without longitudinal data
            "PLASMAPTAU181","PLASMA_NFL","PLASMA_AB42",
            "PET_AV45_SUVR",
            "PET_FDG_SUVR")
Cog_ls <- c("mPACCdigit","mPACCtrailsB","ADNI_MEM","ADNI_EF")
AV1451 <- c("PET_AV1451_SUVR","BRAAK1_SUVR","BRAAK34_SUVR","BRAAK56_SUVR")
MRI_ls <- c("Hippocampus","Temporal_Meta")

LME_Output <- function(fit){
  sum <- summary(fit)
  beta <- sum$coefficients[dim(sum$coefficients)[1],1]
  se <- sum$coefficients[dim(sum$coefficients)[1],2]
  pval <- sum$coefficients[dim(sum$coefficients)[1],5]
  ADV <- i
  Pt <- n
  n_RID <- sum$ngrps
  n_obs <- sum$devcomp$dims[[2]]
  Output <- cbind(ADV, Pt, beta, se, pval, n_RID, n_obs)
  Output
}
data_withfl <- function(ADvar){
  D1 <- DATA[!is.na(DATA[[ADvar]]),]
  a <- D1[D1$M >0,]
  a$TEST <- 1
  a <- a[!duplicated(a$RID),]
  b <- merge(a[,c('RID','TEST')], D1, by = 'RID', all.y = T, no.dups = F)
  c <- b[!is.na(b$TEST),]
  c
}
data_AV1451 <- function(ADvar){
  D1 <- DATA[!is.na(DATA[[ADvar]]),]
  D1 <- D1[order(D1$RID, D1$Month, decreasing = T),]
  D1 <- D1 %>% group_by(RID) %>% mutate(Month_rank = order(Month, decreasing = T))
  D1 <- D1[order(D1$RID, D1$Month, decreasing = F),]
  D1 <- D1 %>% group_by(RID) %>% mutate(Month_rank2 = order(Month, decreasing = F))
  summary(D1$Month_rank)
  D2 <- D1[D1$Month_rank == 1,]
  summary(D2$Month_rank2)
  a <- D2[D2$Month_rank2 >1,]
  a$TEST <- 1
  a <- a[!duplicated(a$RID),]
  b <- merge(a[,c('RID','TEST')], D1, by = 'RID', all.y = T, no.dups = F)
  c <- b[!is.na(b$TEST),]
  DATA <- c
}
{
  # (1) Total
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'Total'
  lme_result$P_Bon <- p.adjust(lme_result$pval, method = 'bonferroni')
  lme_result$P_FDR <- p.adjust(lme_result$pval, method = 'fdr')
  lme_result <- merge(Pt_select, lme_result, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
  lme_result_total <- lme_result
  write.csv(lme_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_total.csv', row.names = F)
  
  # (2) CN/MCI/Dementia
  # CN
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'CN']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'CN']
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'CN']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'CN']
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'CN'
  lme_result_CN <- lme_result
  write.csv(lme_result_CN, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_CN.csv', row.names = F)
  
  # MCI
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'MCI']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'MCI']
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'MCI']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'MCI']
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'MCI'
  lme_result_MCI <- lme_result
  write.csv(lme_result_MCI, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_MCI.csv', row.names = F)
  
  # ZDementia
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'ZDementia']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'ZDementia']
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'ZDementia']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DXBLYES == 'ZDementia']
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'ZDementia'
  lme_result_ZDementia <- lme_result
  write.csv(lme_result_ZDementia, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_ZDementia.csv', row.names = F)
  
  lme_result_S3 <- rbind(lme_result_CN, lme_result_MCI, lme_result_ZDementia)
  lme_result_S3$P_Bon <- p.adjust(lme_result_S3$pval, method = 'bonferroni')
  lme_result_S3$P_FDR <- p.adjust(lme_result_S3$pval, method = 'fdr')
  lme_result_S3 <- merge(Pt_select, lme_result_S3, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
  write.csv(lme_result_S3, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_S3.csv', row.names = F)
  
  # (3) CSF Ab42
  # A-
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A-']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A-']
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A-']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A-']
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'A-'
  lme_result_An <- lme_result
  write.csv(lme_result_An, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_An.csv', row.names = F)
  
  # A+
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A+']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A+']
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A+']
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$CSF_A_status_1098_bl == 'A+']
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'A+'
  lme_result_Ap <- lme_result
  write.csv(lme_result_Ap, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_Ap.csv', row.names = F)
  
  lme_result_A <- rbind(lme_result_Ap, lm_result_An)
  lme_result_A$P_Bon <- p.adjust(lme_result_A$pval, method = 'bonferroni')
  lme_result_A$P_FDR <- p.adjust(lme_result_A$pval, method = 'fdr')
  lme_result_A <- merge(Pt_select, lme_result_A, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
  write.csv(lme_result_A, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_A.csv', row.names = F)
  
  
  # (4) AD continuum
  lme_result <- data.frame()
  for (i in ADV_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DX_AT_status_bl == 'CN A-T-' | data_all$DX_AT_status_bl == 'CN A+T-' | data_all$DX_AT_status_bl == 'CN A+T+' | data_all$DX_AT_status_bl == 'MCI A+T+' | data_all$DX_AT_status_bl == 'ZDementia A+T+',]
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in Cog_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DX_AT_status_bl == 'CN A-T-' | data_all$DX_AT_status_bl == 'CN A+T-' | data_all$DX_AT_status_bl == 'CN A+T+' | data_all$DX_AT_status_bl == 'MCI A+T+' | data_all$DX_AT_status_bl == 'ZDementia A+T+',]
      dta <- data_withfl(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in MRI_ls) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DX_AT_status_bl == 'CN A-T-' | data_all$DX_AT_status_bl == 'CN A+T-' | data_all$DX_AT_status_bl == 'CN A+T+' | data_all$DX_AT_status_bl == 'MCI A+T+' | data_all$DX_AT_status_bl == 'ZDementia A+T+',]
      dta <- data_withfl(i)
      fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  for (i in AV1451) {
    for (n in Pt_ls) {
      DATA <- data_all[data_all$DX_AT_status_bl == 'CN A-T-' | data_all$DX_AT_status_bl == 'CN A+T-' | data_all$DX_AT_status_bl == 'CN A+T+' | data_all$DX_AT_status_bl == 'MCI A+T+' | data_all$DX_AT_status_bl == 'ZDementia A+T+',]
      dta <- data_AV1451(i)
      fit = lmer(scale(dta[[i]]) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
      Row <- LME_Output(fit)
      lme_result <- rbind(lme_result, Row)
    }
  }
  lme_result$Sample <- 'AD continuum'
  lme_result_ADC <- lme_result
  lme_result_ADC$P_Bon <- p.adjust(lme_result_ADC$pval, method = 'bonferroni')
  lme_result_ADC$P_FDR <- p.adjust(lme_result_ADC$pval, method = 'fdr')
  lme_result_ADC <- merge(Pt_select, lme_result_ADC, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
  write.csv(lme_result_ADC, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal_ADC.csv', row.names = F)
}

# merge
lme_result <- rbind(lme_result_total, lme_result_A, lme_result_S3, lme_result_ADC)
lme_result <- merge(Pt_select, lme_result, by = 'Analytes', all.y = T, no.dups = F)
write.csv(lme_result, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_ADv_longitudinal.csv', row.names = F)


#### MRI #####
data <- data_all
MRI  <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA_MRI_longitudinal.csv')
data <- merge(data, MRI[,c(1,65:151)], by = 'id', all.x = T)
CTX <- colnames(data)[grep('CV',colnames(data))]
CTX <- CTX[-c(1,2)]
SCT <- colnames(data)[grep('SV',colnames(data))]
SCT <- SCT[-c(1,2)]
MRI_ls <- c(CTX,SCT)
imaging_match <- read.csv('E:/My works/Glymphatic/2_AD_variable_cross_sectional/Brain_region_formatting_MRI.csv')

library(lme4)
library(lmerTest)
# (1)Total
LME_Output <- function(fit){
  sum <- summary(fit)
  beta <- sum$coefficients[dim(sum$coefficients)[1],1]
  se <- sum$coefficients[dim(sum$coefficients)[1],2]
  pval <- sum$coefficients[dim(sum$coefficients)[1],5]
  ADV <- i
  Pt <- n
  n_RID <- sum$ngrps
  n_obs <- sum$devcomp$dims[[2]]
  Output <- cbind(ADV, Pt, beta, se, pval, n_RID, n_obs)
  Output
}
data_withfl <- function(ADvar){
  D1 <- DATA[!is.na(DATA[[ADvar]]),]
  a <- D1[D1$M >0,]
  a$TEST <- 1
  a <- a[!duplicated(a$RID),]
  b <- merge(a[,c('RID','TEST')], D1, by = 'RID', all.y = T, no.dups = F)
  c <- b[!is.na(b$TEST),]
  c
}
lme_result <- data.frame()
for (i in CTX) {
  for (n in Pt_ls) {
    DATA <- data
    dta <- data_withfl(i)
    fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
    Row <- LME_Output(fit)
    lme_result <- rbind(lme_result, Row)
  }
}
for (i in SCT) {
  for (n in Pt_ls) {
    DATA <- data
    dta <- data_withfl(i)
    fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + ICV + (1 + M|RID), data = dta)
    Row <- LME_Output(fit)
    lme_result <- rbind(lme_result, Row)
  }
}
lme_result$Sample <- 'Total'
lme_result <- merge(imaging_match, lme_result, by.y = 'ADV', by.x = 'FLDNAME', all.x = T, no.dups = F)
lme_result <- merge(Pt_select, lme_result, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
lme_result <- lme_result[!lme_result$Structure=='',]
lme_result$P_Bon <- p.adjust(lme_result$pval, method = 'bonferroni')
lme_result$P_FDR <- p.adjust(lme_result$pval, method = 'fdr')
lme_result_total <- lme_result
write.csv(lme_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_MRI_longitudinal_total.csv', row.names = F)


#### AV45 #####
data <- data_all
AV45 <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA_AV45_longitudinal.csv')
data <- merge(data, AV45[,c(1,45:213)], by = 'id', all.x = T)
ctxl <- colnames(data)[grep('CTX_LH',colnames(data))]
ctxr <- colnames(data)[grep('CTX_RH',colnames(data))]
RT <- colnames(data)[grep('RIGHT_',colnames(data))]
LT <- colnames(data)[grep('LEFT_',colnames(data))]
SCT <- c(RT,LT)
CTX <- c(ctxl, ctxr)
imaging_match <- read.csv('E:/My works/Glymphatic/5_ADv_progression/Brain_region_formatting.csv')

library(lme4)
library(lmerTest)
# (1)Total
LME_Output <- function(fit){
  sum <- summary(fit)
  beta <- sum$coefficients[dim(sum$coefficients)[1],1]
  se <- sum$coefficients[dim(sum$coefficients)[1],2]
  pval <- sum$coefficients[dim(sum$coefficients)[1],5]
  ADV <- i
  Pt <- n
  n_RID <- sum$ngrps
  n_obs <- sum$devcomp$dims[[2]]
  Output <- cbind(ADV, Pt, beta, se, pval, n_RID, n_obs)
  Output
}
data_withfl <- function(ADvar){
  D1 <- DATA[!is.na(DATA[[ADvar]]),]
  a <- D1[D1$M >0,]
  a$TEST <- 1
  a <- a[!duplicated(a$RID),]
  b <- merge(a[,c('RID','TEST')], D1, by = 'RID', all.y = T, no.dups = F)
  c <- b[!is.na(b$TEST),]
  c
}
lme_result <- data.frame()
for (i in CTX) {
  for (n in Pt_ls) {
    DATA <- data
    dta <- data_withfl(i)
    fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
    Row <- LME_Output(fit)
    lme_result <- rbind(lme_result, Row)
  }
}
for (i in SCT) {
  for (n in Pt_ls) {
    DATA <- data
    dta <- data_withfl(i)
    fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
    Row <- LME_Output(fit)
    lme_result <- rbind(lme_result, Row)
  }
}
lme_result$Sample <- 'Total'
lme_result <- merge(imaging_match, lme_result, by.y = 'ADV', by.x = 'Imaging', all.x = T, no.dups = F)
lme_result <- merge(Pt_select, lme_result, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
lme_result <- lme_result[!lme_result$Structure=='',]
lme_result$P_Bon <- p.adjust(lme_result$pval, method = 'bonferroni')
lme_result$P_FDR <- p.adjust(lme_result$pval, method = 'fdr')
lme_result_total <- lme_result
write.csv(lme_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_AV45_longitudinal_total.csv', row.names = F)


#### AV1451 #####
data <- data_all
AV1451 <- read.csv('E:/My works/Proteomics_AD/4_ADv_association/DATA/DATA_AV1451_longitudinal.csv')
ctxl <- colnames(AV1451)[grep('CTX_LH',colnames(AV1451))]
ctxr <- colnames(AV1451)[grep('CTX_RH',colnames(AV1451))]
RT <- colnames(AV1451)[grep('RIGHT_',colnames(AV1451))]
LT <- colnames(AV1451)[grep('LEFT_',colnames(AV1451))]
SCT <- c(RT,LT)
CTX <- c(ctxl, ctxr)
data <- merge(data, AV1451[,c('id',SCT,CTX)], by = 'id', all.x = T)
imaging_match <- read.csv('E:/My works/Glymphatic/5_ADv_progression/Brain_region_formatting.csv')

library(lme4)
library(lmerTest)
# (1)Total
LME_Output <- function(fit){
  sum <- summary(fit)
  beta <- sum$coefficients[dim(sum$coefficients)[1],1]
  se <- sum$coefficients[dim(sum$coefficients)[1],2]
  pval <- sum$coefficients[dim(sum$coefficients)[1],5]
  ADV <- i
  Pt <- n
  n_RID <- sum$ngrps
  n_obs <- sum$devcomp$dims[[2]]
  Output <- cbind(ADV, Pt, beta, se, pval, n_RID, n_obs)
  Output
}
data_AV1451 <- function(ADvar){
  D1 <- DATA[!is.na(DATA[[i]]),]
  D1 <- D1[order(D1$RID, D1$Month, decreasing = T),]
  D1 <- D1 %>% group_by(RID) %>% mutate(Month_rank = order(Month, decreasing = T))
  D1 <- D1[order(D1$RID, D1$Month, decreasing = F),]
  D1 <- D1 %>% group_by(RID) %>% mutate(Month_rank2 = order(Month, decreasing = F))
  summary(D1$Month_rank)
  D2 <- D1[D1$Month_rank == 1,]
  summary(D2$Month_rank2)
  a <- D2[D2$Month_rank2 >1,]
  a$TEST <- 1
  a <- a[!duplicated(a$RID),]
  b <- merge(a[,c('RID','TEST')], D1, by = 'RID', all.y = T, no.dups = F)
  c <- b[!is.na(b$TEST),]
  DATA <- c
}
lme_result <- data.frame()
for (i in CTX) {
  for (n in Pt_ls) {
    DATA <- data
    dta <- data_AV1451(i)
    fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
    Row <- LME_Output(fit)
    lme_result <- rbind(lme_result, Row)
  }
}
for (i in SCT) {
  for (n in Pt_ls) {
    DATA <- data
    dta <- data_AV1451(i)
    fit = lmer(scale(log(dta[[i]])) ~ scale(log(dta[[n]]))*M + AGE + PTGENDER + PTEDUCAT + APOEcarrier + (1 + M|RID), data = dta)
    Row <- LME_Output(fit)
    lme_result <- rbind(lme_result, Row)
  }
}
lme_result$Sample <- 'Total'
lme_result <- merge(imaging_match, lme_result, by.y = 'ADV', by.x = 'Imaging', all.x = T, no.dups = F)
lme_result <- merge(Pt_select, lme_result, by.x = 'Analytes', by.y = 'Pt', all.y = T, no.dups = F)
lme_result <- lme_result[!lme_result$Structure=='',]
lme_result$P_Bon <- p.adjust(lme_result$pval, method = 'bonferroni')
lme_result$P_FDR <- p.adjust(lme_result$pval, method = 'fdr')
lme_result_total <- lme_result
write.csv(lme_result_total, 'E:/My works/Proteomics_AD/4_ADv_association/Results/Results_Pt_AV45_longitudinal_total.csv', row.names = F)











#run plink for adni csf proteins
#
#
for i in {3..274}
do 
/mnt/storage/home1/Huashan3/yl_data/plink2 --bfile /mnt/storage/home1/Huashan3/csd_data/GWAS/GWAS_adni_csf_pro/data/ADNI_merge --out /mnt/storage/home1/Huashan3/csd_data/GWAS/GWAS_adni_csf_pro/out/adni_csf_pro --geno 0.05 --mind 0.05 --maf 0.01 --hwe 1e-6 --glm hide-covar cols=+a1freq --vif 1000  --covar /mnt/storage/home1/Huashan3/csd_data/GWAS/GWAS_adni_csf_pro/data/covar_272pro.txt --covar-col-nums 3,4,5,6,7,8,9,10,11,12,13,14 --covar-variance-standardize --pheno /mnt/storage/home1/Huashan3/csd_data/GWAS/GWAS_adni_csf_pro/data/pheno_272pro.txt --pheno-col-nums ${i}
done







#run Mendelian randomization
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(simex)
library(readxl)
setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916")
#0.0 create folders
#dir.create("exposure_resources")
#dir.create("save_harmonise_dat")
#dir.create("save_report()")
#dir.create("save_result_combine_individ")

#===========PART ONE FORWARD MR ===========#
#===========PART ONE FORWARD MR ===========####
#===========PART ONE FORWARD MR ===========#

setwd("F:/Project/Proj_ADNI_proteomics/GWAS_adni_csf_pro/out")
ls_files <- list.files(pattern = '*.glm.linear.gz')
##adni csf protein clump##
for (file in ls_files) {
  analytes <- substr(substr(file,14,nchar(file)),1,nchar(substr(file,14,nchar(file)))-14)
  #1.read exposure and output clumped files
  setwd("F:/Project/Proj_ADNI_proteomics/GWAS_adni_csf_pro/out")
  expo<- read_exposure_data(file,
                            snp_col  = 'ID',
                            beta_col = 'BETA',
                            se_col = 'SE',
                            eaf_col = 'A1_FREQ',
                            effect_allele_col = 'ALT',
                            other_allele_col = 'REF',
                            pval_col = 'P',
                            samplesize_col = 'OBS_CT',
                            sep = "\t",
                            chr_col = "CHROM",
                            pos_col = "POS",
                            phenotype_col = "pheno")
  
  expo<- rename(expo,rsid=SNP,pval=pval.exposure)#local clump 
  expo_clump <- ld_clump(expo,clump_p = 1E-5,
                         bfile="E:/TEMP_CSD/PROJECT_GWA_hypothalamus/20221027MR_hypo_to_traits/local_clump/EUR",      
                         plink_bin = genetics.binaRies::get_plink_binary())
  expo_clump <- rename(expo_clump,SNP=rsid,pval.exposure=pval)
  setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/exposure_resources")
  write.table(expo_clump,paste0(analytes,'_clump.txt'),row.names=F,quote=F)
}  

#===============start analysis===============#

suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(Cairo))

setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/data/FinnGen_R9_AD")
ls.file.outcome <- list.files(pattern = 'finngen_R9_*')
tbl_ad_data.csv <- read.csv('tbl_ad_data.csv')

setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/exposure_resources")
ls.file.expo.clump <- list.files(pattern = '*_clump.txt')


#if halrted unexpectedly, read result_combine_run.csv file
#setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_result_combine_individ")
#result_combine_f <- read.csv("result_combine_run.finngen_R9_G6_ALZHEIMER.gz.csv")
result_combine_f <- data.frame()
df_norun <- data.frame()
for (j in ls.file.outcome[4]) {
  for ( i in ls.file.expo.clump[253:272]) {
    analytes <- substr(i,1,nchar(i)-10)
    #1.read exposure####
    setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/exposure_resources")
    expo_clump <- fread(i,data.table = F)
    print('COMPLETE: 1.read exposure')
    if (!is.null(expo_clump)) { 
      #2.read outcome####
      setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/data/FinnGen_R9_AD")
      outcome_dat <- read_outcome_data(j,snps=expo_clump$SNP,
                                       snp_col  = 'rsids',
                                       beta_col = 'beta',
                                       se_col = 'sebeta',
                                       pval_col = 'pval', 
                                       #samplesize_col = 'N',
                                       eaf_col = 'af_alt',
                                       effect_allele_col = 'alt',
                                       other_allele_col = 'ref',
                                       sep = '\t'
      )
      
      print('COMPLETE: 2.read outcome')
      if (!is.null(outcome_dat)) { 
        #3.harmonise exposure and outcome####
        trait1.trait2.dat  <- harmonise_data(expo_clump , outcome_dat)
        #trait1.trait2.dat$exposure <- substr(i,1,nchar(i)-14)
        trait1.trait2.dat$outcome <- j
        trait1.trait2.dat$samplesize.outcome <- filter(tbl_ad_data.csv,filename==j)$num_cases+filter(tbl_ad_data.csv,filename==j)$num_controls
        
        setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_harmonise_dat")
        write.csv(trait1.trait2.dat,paste0('harmonised.',analytes,'.',j,'.csv'),row.names = F)
        print('COMPLETE: 3.harmonise exposure and outcome')
        #4.MR analysis and generate OR####
        trait1.trait2.results <- mr(trait1.trait2.dat)
        trait1.trait2.results.withOR = generate_odds_ratios(trait1.trait2.results)
        print(paste0('snps enter mr: ',nrow(filter(trait1.trait2.dat,mr_keep=='TRUE'))))
        print('COMPLETE: 4.MR analysis and generate OR')
        
        if (nrow(filter(trait1.trait2.dat,mr_keep=='TRUE'))>1) { 
          #5.turn the results into one row####
          result_combine <- filter(trait1.trait2.results.withOR, method=='Inverse variance weighted')
          result_MREGGER <- filter(trait1.trait2.results.withOR, method=='MR Egger')
          result_WMe <- filter(trait1.trait2.results.withOR, method=='Weighted median')
          result_WMo<- filter(trait1.trait2.results.withOR, method=='Weighted mode')
          result_SM <- filter(trait1.trait2.results.withOR, method=='Simple mode')
          result_combine <- left_join(result_combine,result_MREGGER,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMe,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMo,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_SM,by=c("id.exposure","id.outcome","outcome","exposure"))
          print('COMPLETE: 5.turn the results into one row, snp_mr_keep>=2')
          #6.MR-EGGER INTERCEPT####
          trait1.trait2.dat.mr_pleiotropy_test = mr_pleiotropy_test(trait1.trait2.dat)
          result_combine$egger_intercept <- trait1.trait2.dat.mr_pleiotropy_test$egger_intercept
          result_combine$egger_intercept.pval <-  trait1.trait2.dat.mr_pleiotropy_test$pval
          print('COMPLETE: 6.MR-EGGER INTERCEPT,snp_mr_keep>=2')
          
          if (length(filter(trait1.trait2.dat,mr_keep=='TRUE')$SNP)>2) { 
            #7.MR-EGGER SIMEX####
            BetaYG = filter(trait1.trait2.dat,mr_keep=='TRUE')$beta.outcome
            BetaXG = filter(trait1.trait2.dat,mr_keep=='TRUE')$beta.exposure
            seBetaYG = filter(trait1.trait2.dat,mr_keep=='TRUE')$se.outcome
            seBetaXG = filter(trait1.trait2.dat,mr_keep=='TRUE')$se.exposure
            Fit2 = lm(BetaYG~BetaXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)
            mod.sim <- simex(Fit2,B=1000,measurement.error = seBetaXG,SIMEXvariable="BetaXG",fitting.method ="quad",asymptotic="FALSE")
            summary(mod.sim)
            simex.beta = summary(mod.sim)[[1]]$jackknife[2]
            simex.se = summary(mod.sim)[[1]]$jackknife[4]
            simex.p = summary(mod.sim)[[1]]$jackknife[8]
            result_combine$isq = Isq(trait1.trait2.results$b[trait1.trait2.results$method == 'MR Egger'],trait1.trait2.results$se[trait1.trait2.results$method == 'MR Egger'])
            result_combine$simex.beta <- simex.beta 
            result_combine$simex.se <- simex.se
            result_combine$simex.p <- simex.p
            print('COMPLETE: 7.MR-EGGER SIMEX,snp_mr_keep)>2')
            
            #8.HETERO####
            trait1.trait2.dat.mr_heterogeneity = mr_heterogeneity(trait1.trait2.dat)
            result_combine$MR_Egger.Q <- trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
            result_combine$MR_Egger.Q_df <- trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
            result_combine$MR_Egger.Q_pval <- trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
            result_combine$Inverse_variance_weighted.Q <- trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_df <- trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_pval <- trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            print('COMPLETE: 8.HETERO,snp_mr_keep>2')
            
          } else {
            #7.MR-EGGER SIMEX####
            result_combine$isq = NA
            result_combine$isq <- as.numeric(result_combine$isq )
            result_combine$simex.beta <- NA
            result_combine$simex.beta <- as.numeric(result_combine$simex.beta)
            result_combine$simex.se <- NA
            result_combine$simex.se <- as.numeric(result_combine$simex.se )
            result_combine$simex.p <- NA
            result_combine$simex.p <- as.numeric(result_combine$simex.p)
            print('COMPLETE: 7.MR-EGGER SIMEX, snp_mr_keep=2')
            
            #8.HETERO####
            trait1.trait2.dat.mr_heterogeneity = mr_heterogeneity(trait1.trait2.dat)
            result_combine$MR_Egger.Q <- NA
            result_combine$MR_Egger.Q <- as.numeric(result_combine$MR_Egger.Q)
            result_combine$MR_Egger.Q_df <- NA
            result_combine$MR_Egger.Q_df <- as.numeric(result_combine$MR_Egger.Q_df)
            result_combine$MR_Egger.Q_pval <- NA
            result_combine$MR_Egger.Q_pval <- as.numeric(result_combine$MR_Egger.Q_pval) 
            result_combine$Inverse_variance_weighted.Q <- trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_df <- trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_pval <- trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            print('COMPLETE: 8.HETERO,snp_mr_keep=2')
          }
          
          #generate report and its path
          setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_report()")
          newpath <- paste0(analytes,"_",j)
          dir.create(newpath)
          mr_report(trait1.trait2.dat,output_path = newpath)  
          print('COMPLETE: generate report and its path,snp_mr_keep>=2')
          
          
        }else{
          #5.turn the results into one row####
          result_combine <- filter(trait1.trait2.results.withOR, method=='Wald ratio')
          result_MREGGER <- filter(trait1.trait2.results.withOR, method=='MR Egger')
          result_WMe <- filter(trait1.trait2.results.withOR, method=='Weighted median')
          result_WMo<- filter(trait1.trait2.results.withOR, method=='Weighted mode')
          result_SM <- filter(trait1.trait2.results.withOR, method=='Simple mode')
          result_combine <- left_join(result_combine,result_MREGGER,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMe,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMo,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_SM,by=c("id.exposure","id.outcome","outcome","exposure"))
          print('COMPLETE: 5.turn the results into one row,snp_mr_keep=1')
          
          #6.MR-EGGER INTERCEPT####
          result_combine$egger_intercept <- NA
          result_combine$egger_intercept <- as.numeric(result_combine$egger_intercept)
          result_combine$egger_intercept.pval <-  NA
          result_combine$egger_intercept.pval <- as.numeric(result_combine$egger_intercept.pval)
          print('COMPLETE: 6.MR-EGGER INTERCEPT,snp_mr_keep=1')
          
          #7.MR-EGGER SIMEX####
          result_combine$isq = NA
          result_combine$isq <- as.numeric(result_combine$isq )
          result_combine$simex.beta <- NA
          result_combine$simex.beta <- as.numeric(result_combine$simex.beta)
          result_combine$simex.se <- NA
          result_combine$simex.se <- as.numeric(result_combine$simex.se )
          result_combine$simex.p <- NA
          result_combine$simex.p <- as.numeric(result_combine$simex.p)
          print('COMPLETE: 7.MR-EGGER SIMEX,snp_mr_keep=1')
          
          #8.HETERO####
          result_combine$MR_Egger.Q <- NA
          result_combine$MR_Egger.Q <- as.numeric(result_combine$MR_Egger.Q)
          result_combine$MR_Egger.Q_df <- NA
          result_combine$MR_Egger.Q_df <- as.numeric(result_combine$MR_Egger.Q_df)
          result_combine$MR_Egger.Q_pval <- NA
          result_combine$MR_Egger.Q_pval <- as.numeric(result_combine$MR_Egger.Q_pval) 
          result_combine$Inverse_variance_weighted.Q <- NA
          result_combine$Inverse_variance_weighted.Q <- as.numeric(result_combine$Inverse_variance_weighted.Q)
          result_combine$Inverse_variance_weighted.Q_df <- NA
          result_combine$Inverse_variance_weighted.Q_df <- as.numeric(result_combine$Inverse_variance_weighted.Q_df)
          result_combine$Inverse_variance_weighted.Q_pval <- NA
          result_combine$Inverse_variance_weighted.Q_pval <- as.numeric(result_combine$Inverse_variance_weighted.Q_pval)
          print('COMPLETE: 8.HETERO,snp_mr_keep=1')
          
          setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_report()")
          newpath <- paste0(analytes,"_",j)
          dir.create(newpath)
          print('COMPLETE: generate path,snp_mr_keep=1')
          
          
        }
        
        #9.other statistics####
        lor <- trait1.trait2.dat$beta.exposure
        af <- trait1.trait2.dat$eaf.exposure
        numbercase=	filter(tbl_ad_data.csv,filename==j)$num_cases
        numbercontrol=filter(tbl_ad_data.csv,filename==j)$num_controls
        ADprevalence=0.05
        
        trait1.trait2.mr_steiger = mr_steiger2(r_exp = get_r_from_pn(trait1.trait2.dat$pval.exposure,trait1.trait2.dat$samplesize.exposure), 
                                               r_out = get_r_from_lor(lor,af,numbercase,numbercontrol,ADprevalence), 
                                               n_exp = trait1.trait2.dat$samplesize.exposure, 
                                               n_out = trait1.trait2.dat$samplesize.outcome)  
        result_combine$steigertest_P <- trait1.trait2.mr_steiger$steiger_test
        result_combine$causal_dir <- trait1.trait2.mr_steiger$correct_causal_direction
        result_combine$steigertest_P.adj <- trait1.trait2.mr_steiger$steiger_test_adj
        result_combine$exp.R2 = trait1.trait2.mr_steiger$r2_exp
        exp.R2=trait1.trait2.mr_steiger$r2_exp
        result_combine$F.stat = (trait1.trait2.dat$samplesize.exposure[1] - dim(trait1.trait2.dat)[1] - 1) / (dim(trait1.trait2.dat)[1]) * exp.R2 / ( 1 - exp.R2 ) 
        print('COMPLETE: 9.other statistics')
        
        #10.power calculate  ####
        ratio = numbercontrol/ numbercase 
        n = numbercase + numbercontrol     
        OR = trait1.trait2.results.withOR$or[trait1.trait2.results.withOR$method == 'Inverse variance weighted'|trait1.trait2.results.withOR$method == 'Wald ratio']
        sig = 0.05
        rsq = exp.R2
        result_combine$power <-  pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*OR-qnorm(1-sig/2))
        print('COMPLETE: 10.power calculate ')
        
        
        setwd(paste0("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_report()/",newpath))
        #11.single_SNP analysis####
        trait1.trait2.single_snp_analysis <- mr_singlesnp(trait1.trait2.dat)
        trait1.trait2.single_snp_analysis.withOR <- generate_odds_ratios(trait1.trait2.single_snp_analysis)
        write.csv(trait1.trait2.single_snp_analysis.withOR,paste0('SSA.',analytes,'.',j,'.csv'),row.names = FALSE)
        print('COMPLETE: 11.single_SNP analysis')
        
        #12.leave-one-out analysis####
        trait1.trait2.dat.mr_leaveoneout <- mr_leaveoneout(trait1.trait2.dat)
        trait1.trait2.dat.mr_leaveoneout.withOR <- generate_odds_ratios(trait1.trait2.dat.mr_leaveoneout)
        write.csv(trait1.trait2.dat.mr_leaveoneout.withOR ,paste0('LOO.',analytes,'.',j,'.csv'),row.names = FALSE)
        print('COMPLETE:12.leave-one-out analysis')
        
        
        #combine results
        setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_result_combine_individ")
        write.csv(result_combine_f,'result_combine_run.csv',row.names=FALSE)
        write.csv(result_combine,paste0('result_combine_I.',analytes,'.',j,'.csv'),row.names = F)
        print('COMPLETE:Individual file output')
        
        result_combine_f <- bind_rows(result_combine_f,result_combine)
        write.csv(result_combine_f,paste0('result_combine_run.',j,'.csv'),row.names=FALSE)
        print('COMPLETE:Merged file output')
        
      } else {
        reason.norun <- 'nosigsnp_in_outcome'
        df_norun_add <- data.frame(analytes=analytes,outcome=j,reason.norun =reason.norun )
        df_norun <- rbind.data.frame(df_norun,df_norun_add)
        #output the record
        setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_result_combine_individ")
        write.csv(df_norun,paste0('nosigsnp_ls.',j,'.csv'),row.names=FALSE)
      }
    } else {
      reason.norun <- 'nosigsnp_in_exposure'
      df_norun_add <- data.frame(analytes=analytes,outcome=j,
                                 reason.norun=reason.norun )
      df_norun <- rbind.data.frame(df_norun,df_norun_add)
      #output the record
      setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0916/save_result_combine_individ")
      write.csv(df_norun,paste0('nosigsnp_ls.',j,'.csv'),row.names=FALSE)
      
    }
  }
}
print('COMPLETE: THE WHOLE LOOP')
write.csv(result_combine_f,'result_combine_run.csv',row.names=FALSE)

#===========PART TWO REVERSE MR ===========#
#===========PART TWO REVERSE MR ===========####
#===========PART TWO REVERSE MR ===========#
#dir.create("exposure_resources")
#dir.create("save_harmonise_dat")
#dir.create("save_report()")
#dir.create("save_result_combine_individ")



setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/data/FinnGen_R9_AD")
ls.file.disease <- list.files(pattern = 'finngen_R9_*')
tbl_ad_data.csv <- read.csv('tbl_ad_data.csv')


for (j in ls.file.disease[1] ) {
  #1.read exposure and output clumped files
  setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/data/FinnGen_R9_AD")
  expo<- read_exposure_data(j,
                            snp_col  = 'rsids',
                            beta_col = 'beta',
                            se_col = 'sebeta',
                            eaf_col = 'af_alt',
                            effect_allele_col = 'alt',
                            other_allele_col = 'ref',
                            pval_col = 'pval', 
                            #samplesize_col = 'OBS_CT',
                            sep = "\t"
  )
  print('COMPLETE: read exposure data')
  expo<- rename(expo,rsid=SNP,pval=pval.exposure)#local clump
  expo_clump <- ld_clump(expo,clump_p = 5E-8,
                         bfile="E:/TEMP_CSD/PROJECT_GWA_hypothalamus/20221027MR_hypo_to_traits/local_clump/EUR",      
                         plink_bin = genetics.binaRies::get_plink_binary())
  print(paste0('COMPLETE: local clump: ',nrow(expo_clump),' snps'))
  expo_clump <- rename(expo_clump,SNP=rsid,pval.exposure=pval)
  setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/exposure_resources")
  write.table(expo_clump,paste0(j,'_clump.txt'),row.names=F,quote=F)
}  



#===============start analysis===============#
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(Cairo))

setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/exposure_resources")
ls_files_disease_clump <- list.files(pattern = '*clump.txt')

setwd("F:/Project/Proj_ADNI_proteomics/GWAS_adni_csf_pro/out")
ls_files <- list.files(pattern = '*.glm.linear.gz')
#ls_files_target <- ls_files[grep("X4179.57|X5694.57|X12853.112|X10479.18|X3216.2|X10980.11|X10351.51|X13388.57",ls_files)]
#ls_files_target_tp <- ls_files[!( ls_files %in% ls_files_target)]
result_combine_f <- data.frame()
df_norun <- data.frame()
for (j in ls_files_disease_clump[1]) {
  for (i in ls_files) {
    analytes <- substr(substr(i,14,nchar(i)),1,nchar(substr(i,14,nchar(i)))-14)
    #1.read exposure
    setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/exposure_resources")
    expo_clump <- fread(j,data.table = F)
    print('COMPLETE: 1.read exposure')
    if (!is.null(expo_clump)) { 
      #2.read outcome
      setwd("F:/Project/Proj_ADNI_proteomics/GWAS_adni_csf_pro/out")
      outcome_dat <- fread(i,data.table = F)
      outcome_dat <- filter(outcome_dat,ID %in% expo_clump$SNP)
      if (nrow(outcome_dat)!=0) { 
        outcome_dat<- format_data(outcome_dat,snp_col  = 'ID',type = 'outcome',
                                  beta_col = 'BETA',se_col = 'SE',eaf_col = 'A1_FREQ',phenotype_col = "pheno",
                                  effect_allele_col = 'ALT',other_allele_col = 'REF',pval_col = 'P',
                                  samplesize_col = 'OBS_CT',chr_col = "CHROM",pos_col = "POS")
        print('COMPLETE: 2.read outcome')
      }
      if (nrow(outcome_dat)!=0) { 
        #3.harmonise exposure and outcome
        trait1.trait2.dat  <- harmonise_data(expo_clump , outcome_dat)
        trait1.trait2.dat$exposure <- j
        trait1.trait2.dat$samplesize.exposure <- filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_cases+filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_controls
        
        setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_harmonise_dat")
        write.csv(trait1.trait2.dat,paste0('harmonised.',j,'.',analytes,'.csv'),row.names = F)
        print('COMPLETE: 3.harmonise exposure and outcome')
        #4.MR analysis and generate OR
        trait1.trait2.results <- mr(trait1.trait2.dat)
        trait1.trait2.results.withOR = generate_odds_ratios(trait1.trait2.results)
        print(paste0('snps enter mr: ',nrow(filter(trait1.trait2.dat,mr_keep=='TRUE'))))
        print('COMPLETE: 4.MR analysis and generate OR')
        
        if (nrow(filter(trait1.trait2.dat,mr_keep=='TRUE'))>1) { 
          #5.turn the results into one row
          result_combine <- filter(trait1.trait2.results.withOR, method=='Inverse variance weighted')
          result_MREGGER <- filter(trait1.trait2.results.withOR, method=='MR Egger')
          result_WMe <- filter(trait1.trait2.results.withOR, method=='Weighted median')
          result_WMo<- filter(trait1.trait2.results.withOR, method=='Weighted mode')
          result_SM <- filter(trait1.trait2.results.withOR, method=='Simple mode')
          result_combine <- left_join(result_combine,result_MREGGER,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMe,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMo,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_SM,by=c("id.exposure","id.outcome","outcome","exposure"))
          print('COMPLETE: 5.turn the results into one row, snp_mr_keep>=2')
          #6.MR-EGGER INTERCEPT
          trait1.trait2.dat.mr_pleiotropy_test = mr_pleiotropy_test(trait1.trait2.dat)
          result_combine$egger_intercept <- trait1.trait2.dat.mr_pleiotropy_test$egger_intercept
          result_combine$egger_intercept.pval <-  trait1.trait2.dat.mr_pleiotropy_test$pval
          print('COMPLETE: 6.MR-EGGER INTERCEPT,snp_mr_keep>=2')
          
          if (length(filter(trait1.trait2.dat,mr_keep=='TRUE')$SNP)>2) { 
            #7.MR-EGGER SIMEX
            BetaYG = filter(trait1.trait2.dat,mr_keep=='TRUE')$beta.outcome
            BetaXG = filter(trait1.trait2.dat,mr_keep=='TRUE')$beta.exposure
            seBetaYG = filter(trait1.trait2.dat,mr_keep=='TRUE')$se.outcome
            seBetaXG = filter(trait1.trait2.dat,mr_keep=='TRUE')$se.exposure
            Fit2 = lm(BetaYG~BetaXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)
            mod.sim <- simex(Fit2,B=1000,measurement.error = seBetaXG,SIMEXvariable="BetaXG",fitting.method ="quad",asymptotic="FALSE")
            summary(mod.sim)
            simex.beta = summary(mod.sim)[[1]]$jackknife[2]
            simex.se = summary(mod.sim)[[1]]$jackknife[4]
            simex.p = summary(mod.sim)[[1]]$jackknife[8]
            result_combine$isq = Isq(trait1.trait2.results$b[trait1.trait2.results$method == 'MR Egger'],trait1.trait2.results$se[trait1.trait2.results$method == 'MR Egger'])
            result_combine$simex.beta <- simex.beta 
            result_combine$simex.se <- simex.se
            result_combine$simex.p <- simex.p
            print('COMPLETE: 7.MR-EGGER SIMEX,snp_mr_keep)>2')
            
            #8.HETERO
            trait1.trait2.dat.mr_heterogeneity = mr_heterogeneity(trait1.trait2.dat)
            result_combine$MR_Egger.Q <- trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
            result_combine$MR_Egger.Q_df <- trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
            result_combine$MR_Egger.Q_pval <- trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'MR Egger']
            result_combine$Inverse_variance_weighted.Q <- trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_df <- trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_pval <- trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            print('COMPLETE: 8.HETERO,snp_mr_keep>2')
            
          } else {
            #7.MR-EGGER SIMEX
            result_combine$isq = NA
            result_combine$isq <- as.numeric(result_combine$isq )
            result_combine$simex.beta <- NA
            result_combine$simex.beta <- as.numeric(result_combine$simex.beta)
            result_combine$simex.se <- NA
            result_combine$simex.se <- as.numeric(result_combine$simex.se )
            result_combine$simex.p <- NA
            result_combine$simex.p <- as.numeric(result_combine$simex.p)
            print('COMPLETE: 7.MR-EGGER SIMEX, snp_mr_keep=2')
            
            #8.HETERO
            trait1.trait2.dat.mr_heterogeneity = mr_heterogeneity(trait1.trait2.dat)
            result_combine$MR_Egger.Q <- NA
            result_combine$MR_Egger.Q <- as.numeric(result_combine$MR_Egger.Q)
            result_combine$MR_Egger.Q_df <- NA
            result_combine$MR_Egger.Q_df <- as.numeric(result_combine$MR_Egger.Q_df)
            result_combine$MR_Egger.Q_pval <- NA
            result_combine$MR_Egger.Q_pval <- as.numeric(result_combine$MR_Egger.Q_pval) 
            result_combine$Inverse_variance_weighted.Q <- trait1.trait2.dat.mr_heterogeneity$Q[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_df <- trait1.trait2.dat.mr_heterogeneity$Q_df[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            result_combine$Inverse_variance_weighted.Q_pval <- trait1.trait2.dat.mr_heterogeneity$Q_pval[trait1.trait2.dat.mr_heterogeneity$method == 'Inverse variance weighted']
            print('COMPLETE: 8.HETERO,snp_mr_keep=2')
          }
          
          #generate report and its path
          setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_report()")
          newpath <- paste0(j,"_",analytes)
          dir.create(newpath)
          mr_report(trait1.trait2.dat,output_path = newpath)  
          print('COMPLETE: generate report and its path,snp_mr_keep>=2')
          
          
        }else{
          #5.turn the results into one row
          result_combine <- filter(trait1.trait2.results.withOR, method=='Wald ratio')
          result_MREGGER <- filter(trait1.trait2.results.withOR, method=='MR Egger')
          result_WMe <- filter(trait1.trait2.results.withOR, method=='Weighted median')
          result_WMo<- filter(trait1.trait2.results.withOR, method=='Weighted mode')
          result_SM <- filter(trait1.trait2.results.withOR, method=='Simple mode')
          result_combine <- left_join(result_combine,result_MREGGER,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMe,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_WMo,by=c("id.exposure","id.outcome","outcome","exposure"))
          result_combine <- left_join(result_combine,result_SM,by=c("id.exposure","id.outcome","outcome","exposure"))
          print('COMPLETE: 5.turn the results into one row,snp_mr_keep=1')
          
          #6.MR-EGGER INTERCEPT
          result_combine$egger_intercept <- NA
          result_combine$egger_intercept <- as.numeric(result_combine$egger_intercept)
          result_combine$egger_intercept.pval <-  NA
          result_combine$egger_intercept.pval <- as.numeric(result_combine$egger_intercept.pval)
          print('COMPLETE: 6.MR-EGGER INTERCEPT,snp_mr_keep=1')
          
          #7.MR-EGGER SIMEX
          result_combine$isq = NA
          result_combine$isq <- as.numeric(result_combine$isq )
          result_combine$simex.beta <- NA
          result_combine$simex.beta <- as.numeric(result_combine$simex.beta)
          result_combine$simex.se <- NA
          result_combine$simex.se <- as.numeric(result_combine$simex.se )
          result_combine$simex.p <- NA
          result_combine$simex.p <- as.numeric(result_combine$simex.p)
          print('COMPLETE: 7.MR-EGGER SIMEX,snp_mr_keep=1')
          
          #8.HETERO
          result_combine$MR_Egger.Q <- NA
          result_combine$MR_Egger.Q <- as.numeric(result_combine$MR_Egger.Q)
          result_combine$MR_Egger.Q_df <- NA
          result_combine$MR_Egger.Q_df <- as.numeric(result_combine$MR_Egger.Q_df)
          result_combine$MR_Egger.Q_pval <- NA
          result_combine$MR_Egger.Q_pval <- as.numeric(result_combine$MR_Egger.Q_pval) 
          result_combine$Inverse_variance_weighted.Q <- NA
          result_combine$Inverse_variance_weighted.Q <- as.numeric(result_combine$Inverse_variance_weighted.Q)
          result_combine$Inverse_variance_weighted.Q_df <- NA
          result_combine$Inverse_variance_weighted.Q_df <- as.numeric(result_combine$Inverse_variance_weighted.Q_df)
          result_combine$Inverse_variance_weighted.Q_pval <- NA
          result_combine$Inverse_variance_weighted.Q_pval <- as.numeric(result_combine$Inverse_variance_weighted.Q_pval)
          print('COMPLETE: 8.HETERO,snp_mr_keep=1')
          
          setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_report()")
          newpath <- paste0(j,"_",analytes)
          dir.create(newpath)
          print('COMPLETE: generate path,snp_mr_keep=1')
          
          
        }
        
        #9.other statistics
        lor <- trait1.trait2.dat$beta.exposure
        af <- trait1.trait2.dat$eaf.outcome
        trait1.trait2.dat$samplesize.exposure <- filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_cases+filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_controls
        
        numbercase=	filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_cases
        numbercontrol=filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_controls
        ADprevalence=0.05
        
        trait1.trait2.mr_steiger = mr_steiger2(r_out = get_r_from_pn(trait1.trait2.dat$pval.outcome,trait1.trait2.dat$samplesize.outcome), 
                                               r_exp = get_r_from_lor(lor,af,numbercase,numbercontrol,ADprevalence), 
                                               n_exp = trait1.trait2.dat$samplesize.exposure, 
                                               n_out = trait1.trait2.dat$samplesize.outcome)  
        result_combine$steigertest_P <- trait1.trait2.mr_steiger$steiger_test
        result_combine$causal_dir <- trait1.trait2.mr_steiger$correct_causal_direction
        result_combine$steigertest_P.adj <- trait1.trait2.mr_steiger$steiger_test_adj
        result_combine$exp.R2 = trait1.trait2.mr_steiger$r2_exp
        exp.R2=trait1.trait2.mr_steiger$r2_exp
        result_combine$F.stat = (trait1.trait2.dat$samplesize.exposure[1] - dim(trait1.trait2.dat)[1] - 1) / (dim(trait1.trait2.dat)[1]) * exp.R2 / ( 1 - exp.R2 ) 
        print('COMPLETE: 9.other statistics')
        
        #10.power calculate  
        n = filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_cases+filter(tbl_ad_data.csv,filename==substr(j,1,nchar(j)-10))$num_controls     
        OR = trait1.trait2.results.withOR$or[trait1.trait2.results.withOR$method == 'Inverse variance weighted'|trait1.trait2.results.withOR$method == 'Wald ratio']
        sig = 0.05
        rsq = exp.R2
        result_combine$power <-  pnorm(sqrt(n*rsq)*OR-qnorm(1-sig/2))
        print('COMPLETE: 10.power calculate ')
        
        
        setwd(paste0("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_report()/",newpath))
        #11.single_SNP analysis
        trait1.trait2.single_snp_analysis <- mr_singlesnp(trait1.trait2.dat)
        trait1.trait2.single_snp_analysis.withOR <- generate_odds_ratios(trait1.trait2.single_snp_analysis)
        write.csv(trait1.trait2.single_snp_analysis.withOR,paste0('SSA.',j,'.',analytes,'.csv'),row.names = FALSE)
        print('COMPLETE: 11.single_SNP analysis')
        
        #12.leave-one-out analysis
        trait1.trait2.dat.mr_leaveoneout <- mr_leaveoneout(trait1.trait2.dat)
        trait1.trait2.dat.mr_leaveoneout.withOR <- generate_odds_ratios(trait1.trait2.dat.mr_leaveoneout)
        write.csv(trait1.trait2.dat.mr_leaveoneout.withOR ,paste0('LOO.',j,'.',analytes,'.csv'),row.names = FALSE)
        print('COMPLETE:12.leave-one-out analysis')
        
        
        #combine results
        setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_result_combine_individ")
        write.csv(result_combine_f,'result_combine_run.csv',row.names=FALSE)
        write.csv(result_combine,paste0('result_combine_I.',j,'.',analytes,'.csv'),row.names = F)
        print('COMPLETE:Individual file output')
        
        result_combine_f <- bind_rows(result_combine_f,result_combine)
        write.csv(result_combine_f,paste0('result_combine_run.',j,'.csv'),row.names=FALSE)
        print('COMPLETE:Merged file output')
        
      } else {
        reason.norun <- 'nosigsnp_in_outcome'
        df_norun_add <- data.frame(disease=j,analytes=analytes,
                                   reason.norun =reason.norun )
        df_norun <- rbind.data.frame(df_norun,df_norun_add)
        #output the record
        setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_result_combine_individ")
        write.csv(df_norun,paste0('nosigsnp_ls.',j,'.csv'),row.names=FALSE)
      }
    } else {
      reason.norun <- 'nosigsnp_in_exposure'
      df_norun_add <- data.frame(disease=j,analytes=analytes,
                                 reason.norun=reason.norun )
      df_norun <- rbind.data.frame(df_norun,df_norun_add)
      #output the record
      setwd("F:/Project/Proj_ADNI_proteomics/MR_analysis/analysis0919_reverseMR/save_result_combine_individ")
      write.csv(df_norun,paste0('nosigsnp_ls.',j,'.csv'),row.names=FALSE)
      
    }
  }
}

