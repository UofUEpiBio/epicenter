
## ----packages, message=FALSE, warning=FALSE---------------------------------------------------------------------------------
library(network)
library(tidyverse)
library(data.table)
library(ergm)
library(ergm.multi)
library(netplot)
library(ggplot2)
library(gridExtra)
library(Rglpk)
library('intergraph')

sessionInfo()

## ----read in -clean--networks -------------------------------------------------------------------------------------------------------
load('network99_f.RData') # a list of 99 networks to fit ERGMs
class(networks) # a list of network objects
load('big_net_mat_f.RData') # a large network (all 99 networks combined) for summaries.
class(bignet_mat) # a network object

bignet_mat # 4212 nodes, 813 HCWs, 3399 residents
813/4212 # 19.3% HCWs -- they have b2degree of 0. 

# Get observed sufficient statistics --- to use for assessing model fit later
# count of b1degree
(sum.b1degree=summary(Networks(networks) ~ b1degree(0:50))) # b1 degree goes as far as 42 
sum(sum.b1degree) #4212
# count of b2degree
sum.b2degree=summary(Networks(networks) ~ b2degree(0:10)) 
sum.b2degree
#Geodesic distance
sum.gdist=ergm.geodistdist(Networks(networks))
sum.gdist[1:7]
#keep only non zero values
sum.gdist.gt0=sum.gdist[sum.gdist!=0]
sum.gdist.gt0

# Calculate proportions
p1.obs=sum.b1degree/sum(sum.b1degree);
p2.obs=sum.b2degree/sum(sum.b2degree);
p2.obs
pgd.obs=sum.gdist.gt0/(sum(sum.gdist.gt0))
plot(p1.obs)
plot(p2.obs)


# Look at a small network
# examine individual networks
set.seed(1234) 
neti=networks[[4]] # a pretty small network
plot(neti)
col=c(hcw = "blue", residents = "green")[neti %v% "group"]
#col[neti %v% "dialysis"==1]="red"
nplot(neti, edge.curvature=0,vertex.color = col,
      vertex.label=NULL)
dg1=summary(neti ~ b1degree(0:10)) # all 27 residents have b1degree of 0
# b1 degree -- degrees of HCWs
# b2 degree -- degrees of residents
dg2=summary(neti ~ b2degree(0:10)) # b1 degree are mostly 0.
dg=summary(neti ~ degree(0:10))
cbind(dg1, dg2, dg) # That's correct.
dg
# try fitting a simple ergm on the small network
ef1=ergm(neti ~ edges)
gg=gof(ef1~ degree + esp + distance) # goodness of fit - this works
gg=gof(ef1)
plot(gg)
dg/sum(dg)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Making the two new plots using vertex attributes

# 1) distribution of network sizes (by node type) and 
# 2) correlation matrix of resident status (e.g., # of individuals with MDRO and under wound care)

# Count the number of HCWs and residents in each network
counts=lapply(networks, function(x){data.table(hcw=table(x %v% 'group')[1], 
residents=table(x %v% 'group')[2] )}) %>% rbindlist
head(counts)

#  plot 
dplot=counts[,.N, keyby = .(hcw, residents)]
dplot %>% head
dplot$N=as.factor(dplot$N)
ggplot(dplot, aes(residents, hcw, fill=N))+geom_tile()+scale_fill_manual(values=c('lightgreen', 'green4'))+theme_bw()+ylab('#HCWs')+xlab('#Residents')#+ylim(0, 18)
#ggsave('network_size.pdf', width=9, height=4)

#view nodal attributes of bignet_mat
bignet_mat %v% "group" 
bignet_mat %v% "type"
bignet_mat %v% "diabetes"

#get all the nodal attributes and combine them into a data frame
dd=as.data.frame(cbind(bignet_mat %v% "group", bignet_mat %v% "type", bignet_mat %v% "diabetes", bignet_mat %v% "wound_care", bignet_mat %v% "dialysis",
 bignet_mat %v% "wheelchair", bignet_mat %v% "bedridden", bignet_mat %v% "ventilator", bignet_mat %v% "mdro"))
dim(dd) # 4212 rows, 9 columns
head(dd)
#variable names
names(dd)=c("group", "type", "diabetes", "wound_care", "dialysis", "wheelchair", "bedridden", "ventilator", "mdro")
setDT(dd)
table(dd$mdro) # 183 residents with mdro
dd[group=='residents' & mdro==1, .N, by=.(diabetes, wound_care, dialysis, wheelchair, bedridden, ventilator)] %>% arrange(N)

res=dd[group=='residents'][,group:=NULL]
head(res)
dim(res) # 3399 residents, 8 attributes
mm=as.data.frame(res[,-1]) 
str(mm)
#convert all the attributes into numeric
mm=apply(mm, 2, as.numeric)

#calculate correlation matrix of mm
round(cor(mm),2)
library(reshape2)
mcor=cor(mm)
mcor[lower.tri(mcor)]=NA
mcor[mcor>0.99]=NA # diag
mcor %>% mrange
dplot=melt(mcor, na.rm=TRUE)
head(dplot)

ggplot(data = dplot, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "seagreen3", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.4,0.4), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ xlab('')+ylab('')+ 
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
#ggsave('Correlation_residents_care.pdf', width=7, height=7)

#high correlation between diabetes and dialysis, wheelchair and bedridden (negative!), ventilator and bedridden, mdro and ventilator



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Now load fitted ERGM models (mode 9 and model 10)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('res9.RData')
# Refit model 9
default_control <- list( 
  init=as.numeric(coef(res9)),
  seed    = 44485)

# try a simpler model
res <- ergm(
  Networks(networks) ~ 
    N(~edges, ~log(n))+
    N(~b2concurrent)+
    N(~b2starmix(2, "type"), ~log(n))+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf))

summary(res)

# print the model using textreg
texreg::knitreg( res, single.row = TRUE, scalebox   = .8)


res9 <- ergm(
  Networks(networks) ~ 
    N(~edges, ~log(n))+
    N(~b2concurrent)+
    b2starmix(2, "type")+
    gwb1dsp(0.25, fixed=TRUE)+ #gwb2dsp(0.75, fixed=TRUE)+ 
    b1starmix(2, "ventilator", diff = TRUE) +
    b1starmix(2, "wound_care", diff = TRUE) +
    b1starmix(2, "bedridden", diff = TRUE) +
    b1starmix(2, "diabetes", diff=TRUE) +
    b1starmix(2, "mdro", diff=TRUE)+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))


gofn9=gofN(res9)
gofn9
plot(gofn9)

#pdf('GOFn_model9.pdf')
#plot(gofn9)
#dev.off()
summary(res9)

#sim9=simulate(res9, nsim=1000)
load('sim9_1000.RData')

# Calculate b1dgree, b2degree, and geodesic distance for each simulated network
# b1 degree
getb1dgree=function(neti){
out=summary(neti ~ b1degree(0:50))
return(out/sum(out))}
# b2 degree
getb2dgree=function(neti){
out=summary(neti ~ b2degree(0:10))
return(out/sum(out))}
# geodesic distance
# use ergm.geodistdist( )


# get these statistics for all simulated networks
b1degree.1000=lapply(sim9, getb1dgree) %>% do.call(rbind, .)
b2degree.1000=lapply(sim9, getb2dgree) %>% do.call(rbind, .)
gd.1000=lapply(sim9, ergm.geodistdist) %>% do.call(rbind, .)

save(b1degree.1000, file='b1degree.1000.RData')
#save(b2degree.1000, file='b2degree.1000.RData')
#save(gd.1000, file='gd.1000.RData')


load('b1degree.1000.RData')
load('b2degree.1000.RData')
load('gd.1000.RData')

# Assess GOF using b1degree
# stack all columns in b1degree.1000 into one big column
ddd=NULL
for (i in 1:51){ # b1 degree 0 ~ 50
ddd=rbind(ddd, data.frame(Obs=p1.obs[i], Sim=b1degree.1000[,i], index=i-1))
}
head(ddd)
# for each unique value of index, make a box plot of Sim using ggplot2
p1=ggplot(ddd, aes(x=index, y=Sim, group=index))+geom_boxplot()+xlim(0,10)+geom_point(aes(x=index, y=Obs), color='red')+theme_bw()
pdf('GOF_b1degree.pdf')
p1+ylim(0, 0.03)+scale_x_continuous(name='b1 degree', breaks=0:20, limits=c(0, 20))+
ylab('Proportioin')+geom_text(x=5, y=0.03, label='proportion of b1 degree 0 is 0.806')
dev.off()

# stack all columns in b2degree.1000 into one big column
ddd=NULL
for (i in 1:10){
ddd=rbind(ddd, data.frame(Obs=p2.obs[i], Sim=b2degree.1000[,i], index=i-1))
}
head(ddd)
# for each unique value of index, make a box plot of Sim using ggplot2
p2=ggplot(ddd, aes(x=index, y=Sim, group=index))+geom_boxplot()+xlim(0,10)+geom_point(aes(x=index, y=Obs), color='red')+theme_bw()
pdf('GOF_b2degree.pdf')
p2+geom_point(aes(x=0.07, y=0.193), color='black')+scale_x_continuous(name='b2 degree', breaks=0:10, limits=c(0, 10))+ylab('Proportioin')
dev.off()

# Now geodestic distance
head(gd.1000) 
dim(gd.1000)
#keep only colums that are not all 0
gd.1000=gd.1000[,colSums(gd.1000)!=0]
dim(gd.1000)
head(gd.1000)
gd.1000=gd.1000[,1:20]

p.gd.1000=t(apply(gd.1000, 1, function(x){x/sum(x)}))
dim(p.gd.1000)
rowSums(p.gd.1000)
p.gd.1000=as.data.frame(p.gd.1000)
head(p.gd.1000)
sum(pgd.obs)
pgd.obs=pgd.obs[1:9]
pgd.obs=pgd.obs/sum(pgd.obs)
sum(pgd.obs)
# stack all columns in b2degree.1000 into one big column
ddd=NULL
for (i in 1:9){
ddd=rbind(ddd, data.frame(Obs=pgd.obs[i], Sim=p.pd.1000[,i], dist=i))
}
head(ddd)
# for each unique value of dist, make a box plot of Sim using ggplot2
p3=ggplot(ddd, aes(x=dist, y=Sim, group=dist))+geom_boxplot()+xlim(0,10)+geom_point(aes(x=dist, y=Obs), color='red')+theme_bw()
pdf('GOF_gdist.pdf')
p3
dev.off()



anys=apply(gd.1000,2,mean)
pgd.obs

ddd=NULL
for (i in 1:10){
ddd=rbind(ddd, data.frame(Obs=p2.obs[i], Sim=b2degree.1000[,i], index=i-1))
}
head(ddd)
# for each unique value of index, make a box plot of Sim using ggplot2
p2=ggplot(ddd, aes(x=index, y=Sim, group=index))+geom_boxplot()+xlim(0,10)+geom_point(aes(x=index, y=Obs), color='red')+theme_bw()
pdf('GOF_b2degree.pdf')
p2+geom_point(aes(x=0.07, y=0.193), color='black')+scale_x_continuous(name='b2 degree', breaks=0:10, limits=c(0, 10))+ylab('Proportioin')
dev.off()


length(pgd.obs) #4212
pgd.obs[1:5] 
summary(pgd.obs)

plot(pgd.obs)

dim(gd.1000)
names(gd.1000)

# Now make box plot against observed
dim(b2degree.1000)
head(b2degree.1000)
miqr(b2degree.1000[,5])
p2.obs
p1.obs

apply(b2degree.1000, 2, median)
apply(b1degree.1000, 2,  median)



class(net)
class(net[[2]])
length(net[[1]])

summary(net)
names(net)


texreg::knitreg( list(res7a, res8, res9), single.row = TRUE, scalebox   = .8)
as.numeric(coef(res9))

ini=c( 3.4375035568, -1.0918000923, -0.4256032660, -0.2157707568, -0.3205728296 ,-0.1791093346      ,    -Inf  ,        -Inf,
       -0.6274064325, -0.9046699707, -1.0780538275 , 0.0041520358 , 0.0337860322,  0.0089902734 , 0.0706534435 ,0,0,
       0.0230620999 , 0.0326673457   ,     -Inf   ,       -Inf)


default_control <- list( 
  init=ini,
  seed    = 44485)

res10 <- ergm(
  Networks(networks) ~ 
    N(~edges, ~log(n))+
    N(~b2concurrent)+
    b2starmix(2, "type")+
    gwb1dsp(0.25, fixed=TRUE)+ #gwb2dsp(0.75, fixed=TRUE)+ 
    b1starmix(2, "ventilator", diff = TRUE) +
    b1starmix(2, "wound_care", diff = TRUE) +
    b1starmix(2, "wheelchair", diff = TRUE) +
    b1starmix(2, "diabetes", diff=TRUE) +
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

#save(res10, file='res10.RData')
load('res10.RData')
summary(res10)

gn=gofN(res10)
pdf('gg.pdf')
plot(gn)
dev.off()

#mcmc.diagnostics(res10)

texreg::knitreg( list(res7a, res8, res9, res10), single.row = TRUE, scalebox   = .8)
as.numeric(coef(res9))

as.numeric(coef(res10))
ini=c(3.431193033, -1.075791164, -0.445590184 ,-0.224438208 ,-0.342860801, -0.172101013,    -Inf  ,   -Inf,
      -0.655549717, -0.893647216 ,-1.111844295 , 0.006933418,
  0.043404094 , 0.003345014 , 0.035416426,  0.024568491 , 0.014448774 , 0.014466141 , 0.024784239 ,0,0,        -Inf   ,      -Inf)

default_control <- list( 
  init=ini,
  seed    = 44485)

res11 <- ergm(
  Networks(networks) ~ 
    N(~edges, ~log(n))+
    N(~b2concurrent)+
    b2starmix(2, "type")+
    gwb1dsp(0.25, fixed=TRUE)+ #gwb2dsp(0.75, fixed=TRUE)+ 
    b1starmix(2, "ventilator", diff = TRUE) +
    b1starmix(2, "wound_care", diff = TRUE) +
    b1starmix(2, "wheelchair", diff = TRUE) +
    b1starmix(2, "diabetes", diff=TRUE) +
    b1starmix(2, "dialysis", diff = TRUE) +
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res11, file='res11.RData')
mcmc.diagnostics(res11)
texreg::knitreg( list(res7a, res8, res9, res10, res11), single.row = TRUE, scalebox   = .8)

gof(res11)
