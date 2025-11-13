# load packages
 require(Matching)
 require(lme4)
### load data set
# insert here your data
data <- read.csv(...)

# covariates: all variables in SDO (hospital discharge sheet) for all sardinian hospitals
# outcome variable: apgar score
# treatment variable: caesarean section (names cs; 0=no; 1=yes)

###### generate target subsample (=non complicated pregancies)

# eg (BMJ 2010)  15 and 44 years who had a singleton birth, and to NHS trusts whose obstetric units had more than 1000 deliveries in the 12 month period;

# eg (Health Affairs 2013) nulliparous, term (>=37 weeks) , singleton, and vertex births

# eg (Arch Pediatrics 2000) preterm singletons in any position

# subset of non-complicated pregnancies ( same as health_affairs_2013)

dha<-data[data$Vitalita==1 & data$Presenza.malformazione==2 & data$momage>15 & data$momage<44 & data$Genere.parto==1 & data$Numero.parti.precedenti==0 & data$Presentazione==1 & data$Eta.gestazionale>30,]# = health_affairs_2013 + preterm (tranne severe preterm)#14757 



##### CHOOSE SUBSAMPLE AND ORDER IT (IMPORTANT FOR CORRECT WITHIN MATCHING)


d<-dha

tb2<-table(factor(d$Presidio))
d$ord.Presidio<-factor(d$Presidio,levels=names(tb2[order(tb2, decreasing = TRUE)]))
d<-d[order(d$ord.Presidio),] 


#### ICC calculation

 fiticc <- glmer(formula = lowapg ~ momage_jbm  
   + m_educ_2 + m_educ_3 + m_educ_miss  
   + babyweight25 + babyweight40 
   + labind
   + preterm + latenorm
   + Decorso.gravidanza + (1 | Presidio), data=d, family = binomial(link="logit"))
   
   icc<-as.numeric(summary(fiticc)$varcor)/((pi^2)/3+(as.numeric(summary(fiticc)$varcor)))
    
    fiticc <- glmer(formula = lowapg ~ momage_jbm  
   + m_educ_2 + m_educ_3 + m_educ_miss  
   + babyweight25 + babyweight40 
   +labind
   + preterm + latenorm
   + Decorso.gravidanza + (1 | Presidio), data=d[d$cs==1,], family = binomial(link="logit"))
   
     icct1<-as.numeric(summary(fiticc)$varcor)/((pi^2/3)+(as.numeric(summary(fiticc)$varcor)))
   
       fiticc <- glmer(formula = lowapg ~ momage_jbm  
   + m_educ_2 + m_educ_3 + m_educ_miss  
   + babyweight25 + babyweight40 
   +labind
   + preterm + latenorm
   + Decorso.gravidanza + (1 | Presidio), data=d[d$cs==0,], family = binomial(link="logit"))
   

  icct0<-as.numeric(summary(fiticc)$varcor)/((pi^2/3)+(as.numeric(summary(fiticc)$varcor)))


######
# propensity score model: logit (used in matching approaches ABC below)
######


# logit ps
    fit1<-glm(formula =
    cs ~ momage_jbm  
   + m_educ_2 + m_educ_3  
   + babyweight25 + babyweight40 
   + labind
   + Eta.gestazionale
   + Decorso.gravidanza 
   , data=d, family = binomial(link = "logit"))
 d$log.ps<-fitted(fit1)
 


###### "MATCH" BEFORE (all obs...)

attbef<-mean(d$lowapg[d$cs==1])-mean(d$lowapg[d$cs==0])

p0<-(length(d$lowapg[d$cs==1])*mean(d$lowapg[d$cs==1])+length(d$lowapg[d$cs==0])*mean(d$lowapg[d$cs==0]))/dim(d)[1]

sebef<-sqrt((p0*(1-p0)*(1/length(d$lowapg[d$cs==1])+1/length(d$lowapg[d$cs==0]))))

bal = MatchBalance(cs ~ momage20 + momage2024 + momage2529+momage3035+momage35    + m_educ_1 + m_educ_2 + m_educ_3 + m_educ_miss
 + babyweight25 + babyweight2540 + babyweight40 +labind
   + Eta.gestazionale
   + Decorso.gravidanza
, data=d,ks = FALSE, nboots=0, print.level=0) 

asam_bef    <-vector();for(i in 1:15){asam_bef[[i]] <- bal$BeforeMatching[[i]]$sdiff}

mean(abs(asam_bef))

####### a) pooled matching


rr = Match(Y=d$lowapg, Tr=d$cs, X=d$log.ps, caliper=0.25, M=1, replace=TRUE,ties=TRUE)

#dropped units by hosp
mndrop<-table(d[rr$index.dropped,]$ord.Presidio)

#matched dataset
 md<-rbind(d[rr$index.treated,],d[rr$index.control,])

#adjusted att and se (after pooled match)
attaftpm<-rr$est
 
#sepm<-rr$se.standard
m0<-lm(formula = lowapg ~ cs, data=md, weights=c(rr$weights,rr$weights))#summary(m0)
m0.vcovCL<-cluster.vcov(m0, md$Presidio)
sepm<-coeftest(m0, m0.vcovCL)[4]

# average number of replicates of a control in md
ff<-(aggregate(rr$weights ~ rr$index.control , data = cbind(rr$weights,rr$index.control) , sum))

arA<-sum(ff[,"rr$weights"][which(ff[,"rr$weights"]>1)]-1)/
         length(unique(ff[,"rr$index.control"]))#in general
#arA<-(sum(rr$weights)-length(unique(rr$index.control)))/length(unique(rr$index.control))# reduction when TIES=FALSE


# balance after pooled matching:
bal_apm = MatchBalance(cs ~ momage20 + momage2024 + momage2529+momage3035+momage35 +
m_educ_1 + m_educ_2 + m_educ_3 + m_educ_miss
 + babyweight25 + babyweight2540 + babyweight40 +labind
   + Eta.gestazionale
   + Decorso.gravidanza
, data=d, match.out=rr, ks = FALSE, nboots=0, print.level=0) 


asam_aftpm    <-vector();for(i in 1:15){asam_aftpm[[i]] <- bal_apm$AfterMatching[[i]]$sdiff}

mean(abs(asam_aftpm))


####### b) matching within hospital


cumfreq<-c(0,cumsum(table(d$ord.Presidio)))

#m.data.with<-list();um.treat.with<-list()

inddrop<-vector(); indtreat<-vector(); indcontr<-vector(); weights<-vector(); indatt<-vector();indse<-vector()

for(i in 1:(nlevels(d$ord.Presidio)-1))
   {	
	        rr <- Match(
Y      =d[d$ord.Presidio==levels(d$ord.Presidio)[i],]$lowapg,
Tr     =d[d$ord.Presidio==levels(d$ord.Presidio)[i],]$cs,
X      =d[d$ord.Presidio==levels(d$ord.Presidio)[i],]$log.ps,
caliper=0.25*sd(d$log.ps)/sd(d[d$ord.Presidio==levels(d$ord.Presidio)[i],]$log.ps), M=1, replace=TRUE,ties=TRUE)


inddrop <-c(inddrop,cumfreq[i]+rr$index.dropped)
indtreat<-c(indtreat,cumfreq[i]+rr$index.treated)
indcontr<-c(indcontr,cumfreq[i]+rr$index.control)
indatt<-c(indatt,rr$est)
indse<-c(indse,rr$se.standard)
weights<-c(weights,rr$weights)	
	}


mwndrop<-table(d[inddrop,]$ord.Presidio)

mwd<-rbind(d[indtreat,],d[indcontr,])

attaftwm<-sum((mwd$lowapg[mwd$cs==1]-mwd$lowapg[mwd$cs==0])*weights/sum(weights))#mean(mwd$lowapg[mwd$cs==1]-mwd$lowapg[mwd$cs==0])#only if TIES=FALSE

#sewm<-sqrt(sum(indse*as.numeric(table(d$ord.Presidio[d$cs==1])[1:19]-mwndrop[1:19])/sum(table(d$ord.Presidio[d$cs==1])[1:19]-mwndrop[1:19])^2))

m0<-lm(formula = lowapg ~ cs, data=mwd, weights=c(weights,weights))#summary(m0)
m0.vcovCL<-cluster.vcov(m0, mwd$Presidio)
sewm<-coeftest(m0, m0.vcovCL)[4]

#average number of replicates
ff<-(aggregate(weights ~ indcontr , data = cbind(weights,indcontr) , sum))

arB<-sum(ff[,"weights"][which(ff[,"weights"]>1)]-1)/
         length(unique(ff[,"indcontr"]))#in general



# balance after within hospital matching 

bal_amw = MatchBalance(cs ~ momage20 + momage2024 + momage2529+momage3035+momage35 + m_educ_1 + m_educ_2 + m_educ_3 + m_educ_miss
 + babyweight25 + babyweight2540 + babyweight40 +labind
   + Eta.gestazionale + Decorso.gravidanza
   , data=mwd,ks = FALSE, weights=c(weights,weights), nboots=0, print.level=0) 


asam_aftmw    <-vector();for(i in 1:15){asam_aftmw[[i]] <- bal_amw$BeforeMatching[[i]]$sdiff}

mean(abs(asam_aftmw))


####### 
####### c) matching within hospital and remaining units between hospital (preferential within)


umwd<-rbind(d[inddrop,],d[which(d$cs==0),])#treated not matched within and controls

rrc <- Match(Y=umwd$lowapg,Tr=umwd$cs, X=umwd$log.ps,
caliper=0.25*sd(d$log.ps)/sd(umwd$log.ps), M=1, replace=TRUE,ties=TRUE)

mbdata<-rbind(umwd[rrc$index.treated,],umwd[rrc$index.control,])#treated matched between (but not within) and matched controls

mwbd<-rbind(mwd,mbdata) # treated matched within or between (and matched controls)

#dropped after within between match
mwbndrop<-table(d[rrc$index.dropped,]$ord.Presidio)

#att after within between match
attaftmwb<-sum((mwbd$lowapg[mwbd$cs==1]-mwbd$lowapg[mwbd$cs==0])*c(weights,rrb$weights)/sum(c(weights,rrb$weights)))#mean(mwbd$sim.lowapg[mwbd$sim.cs==1]-mwbd$sim.lowapg[mwbd$sim.cs==0])#only if TIES=FALSE

#sepw<-sqrt(sum(c(indse,rrb$se.standard)*c(as.numeric(table(d$ord.Presidio[d$cs==1])[1:19]-mwndrop[1:19]),sum(mbdata$cs))/sum(c(as.numeric(table(d$ord.Presidio[d$cs==1])[1:19]-mwndrop[1:19]),sum(mbdata$cs)))^2))

m0<-lm(formula = lowapg ~ cs, data=mwbd, weights=rep(c(weights,rrb$weights),2))#summary(m0)
m0.vcovCL<-cluster.vcov(m0, mwbd$Presidio)
sepw<-coeftest(m0, m0.vcovCL)[4]

# average number of replicates
indcontrc<-c(indcontr,rrc$index.control)
weightsc<-c(weights,rrc$weights)

ff<-(aggregate(weightsc ~ indcontrc , data = cbind(weightsc,indcontrc) , sum))

arC<-sum(ff[,"weightsc"][which(ff[,"weightsc"]>1)]-1)/
         length(unique(ff[,"indcontrc"]))#in general


####### balance after within - between matching 


bal_amwb = MatchBalance(cs ~ momage20 + momage2024 + momage2529+momage3035+momage35 +
m_educ_1 + m_educ_2 + m_educ_3 + m_educ_miss
 + babyweight25 + babyweight2540 + babyweight40 +labind
   + Eta.gestazionale
   + Decorso.gravidanza
, data=mwbd,ks = FALSE, weights=c(weights,weights,rrb$weights,rrb$weights), nboots=0, print.level=0)

asam_aftmwb    <-vector();for(i in 1:15){asam_aftmwb[[i]] <- bal_amwb$BeforeMatching[[i]]$sdiff}
mean(abs(asam_aftmwb))

###### different ps model: dummy logit

 # dummy logit ps
   
    fit2<-glm(formula = cs ~ momage_jbm  
   + m_educ_2 + m_educ_3 + m_educ_miss  
   + babyweight25 + babyweight40 
   + labind
   + Eta.gestazionale
   #+ multip.nocs 
   + Decorso.gravidanza 
   + ord.Presidio -1 # dummies, -1 per togliere intercetta
   , data=d, family = binomial(link = "logit"))
    
    d$dlog.ps<-fitted(fit2)


####### a) pooled matching


rrD = Match(Y=d$lowapg, Tr=d$cs, X=d$dlog.ps, caliper=0.25, M=1, replace=TRUE,ties=TRUE)

#dropped units by hosp
mndropD<-table(d[rrD$index.dropped,]$ord.Presidio)

#matched dataset
mdD<-rbind(d[rrD$index.treated,],d[rrD$index.control,])


#adjusted att and se for D
attaftD<-rrD$est 
#seD<-rrD$se.standard
m0<-lm(formula = lowapg ~ cs, data=mdD, weights=c(rrD$weights,rrD$weights)) #summary(m0)
m0.vcovCL<-cluster.vcov(m0, mdD$Presidio)
seD<-coeftest(m0, m0.vcovCL)[4]

#average number of replicates

ff<-(aggregate(rrD$weights ~ rrD$index.control , data = cbind(rrD$weights,rrD$index.control) , sum))

arD<-sum(ff[,"rrD$weights"][which(ff[,"rrD$weights"]>1)]-1)/
         length(unique(ff[,"rrD$index.control"]))#in general



# balance after pooled matching:


bal_apmD = MatchBalance(cs ~ momage20 + momage2024 + momage2529+momage3035+momage35 +
m_educ_1 + m_educ_2 + m_educ_3 + m_educ_miss
 + babyweight25 + babyweight2540 + babyweight40 +labind
   + Eta.gestazionale
   + Decorso.gravidanza
, data=d, match.out=rrD, ks = FALSE, nboots=0, print.level=0) 


asam_aftD    <-vector();for(i in 1:15){asam_aftD[[i]] <- bal_apmD$AfterMatching[[i]]$sdiff}

mean(abs(asam_aftD))

###### choose ps model: multilevel logit (approccio E)

    #multilevel ps

 fit3 <- glmer(formula = cs ~ momage_jbm  
   + m_educ_2 + m_educ_3 + m_educ_miss  
   + babyweight25 + babyweight40 
   +labind
   + Eta.gestazionale
   #+ preterm + latenorm 
   + Decorso.gravidanza + (1 | ord.Presidio), data=d, family = binomial(link="logit"))

  
    d$ml.ps<-fitted(fit3)



####### matching using multilevel ps


rrE = Match(Y=d$lowapg, Tr=d$cs, X=d$ml.ps, caliper=0.25, M=1, replace=TRUE,ties=TRUE)

#dropped units by hosp
mndropE<-table(d[rrE$index.dropped,]$ord.Presidio)
#matched dataset
mdE<-rbind(d[rrE$index.treated,],d[rrE$index.control,])

#adjusted att and se (after pooled match)
attaftE<-rrE$est 
#seE<-rrE$se.standard
m0<-lm(formula = lowapg ~ cs, data=mdE, weights=c(rrE$weights,rrE$weights)) #summary(m0)
m0.vcovCL<-cluster.vcov(m0, mdE$Presidio)
seE<-coeftest(m0, m0.vcovCL)[4]

#average number of replicates

ff<-(aggregate(rrE$weights ~ rrE$index.control , data = cbind(rrE$weights,rrE$index.control) , sum))

arE<-sum(ff[,"rrE$weights"][which(ff[,"rrE$weights"]>1)]-1)/
         length(unique(ff[,"rrE$index.control"]))#in general



# balance after pooled matching:


bal_apmE = MatchBalance(cs ~ momage20 + momage2024 + momage2529+momage3035+momage35 +
m_educ_1 + m_educ_2 + m_educ_3 + m_educ_miss
 + babyweight25 + babyweight2540 + babyweight40 +labind
   + Eta.gestazionale
   + Decorso.gravidanza 
, data=d, match.out=rrE, ks = FALSE, nboots=0, print.level=0) 


asam_aftE    <-vector();for(i in 1:15){asam_aftE[[i]] <- bal_apmE$AfterMatching[[i]]$sdiff}

mean(abs(asam_aftE))





    