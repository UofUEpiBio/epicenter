## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)
# a few personal functions to make things easy
naton=function(x){y=x;y[is.na(x)]=0; return(y)} # convert NA to 0
nuni=function(x){length(unique(x))} # count unique values
findmiss=function(dat){M=apply(dat,2,function(x){sum(is.na(x))}); M[M>0]} # count missings in all columns
#setwd("/uufs/chpc.utah.edu/common/home/u0942431/school")
setwd("C:/Users/u0942431/Box/Epicenter LTCF Data Analysis/Chong/Analysis_New")
options(width = 170)

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(network)
library(tidyverse)
library(data.table)
library(ergm)
library(ergm.multi)
library(netplot)
library(ggplot2)
library(gridExtra)
library(parallel)
library(Rglpk)


## -----------------------------------------------------------------------------
packageVersion("network")
packageVersion("ergm")
packageVersion("ergm.multi")


## ----read-clean---------------------------------------------------------------
load('network93_f2.RData') # a list of 93 networks to fit ERGMs
load('big_net_mat_f2.RData') # a large network (all 99 networks combined) for summaries.
load('dd_nodal_attr.RData')
dim(dd)
glimpse(dd)
setDT(dd)

networks1=networks[-57]
length(networks1)

dd$NetID=as.numeric(dd$NetID)
dd=dd[NetID!=57]
dd$NetID[dd$NetID>57]=dd$NetID[dd$NetID>57]-1
unique(dd$NetID)

Nurse=dd[, n:=sum(Nurse==1), NetID][n>1, unique(NetID)] %>% as.numeric # requiring at least 2 nurses in a network 
CNA=dd[, n:=sum(CNA==1), NetID][n>1, unique(NetID)] %>% as.numeric
PTOTRT=dd[, n:=sum(PTOTRT==1), NetID][n>1, unique(NetID)] %>% as.numeric
Other=dd[, n:=sum(Other==1), NetID][n>1, unique(NetID)] 
nuni(c(Other, PTOTRT)) # 46

wh=intersect(Nurse, PTOTRT)
length(wh)
dcheck=dd[NetID %in% wh & group=='hcw', .(NetID, Nurse, PTOTRT, Other, CNA)]
dcheck[, count:=sum(CNA=="1")+sum(Other=="1"), NetID]
PTOTRT.nurse=as.numeric(dcheck[count>1, unique(NetID)])
length(PTOTRT.nurse) # all 30 networks are fine!

dialysis=dd[, n:=sum(dialysis=='Yes'), NetID][n>0, unique(NetID)]
wound_care=dd[, n:=sum(wound_care=='Yes'), NetID][n>0, unique(NetID)] 
mdro=dd[, n:=sum(mdro=='Yes'), NetID][n>0, unique(NetID)] 
bedridden=dd[, n:=sum(bedridden=='Yes'), NetID][n>0, unique(NetID)]  # requiring at least 2 nurses in a network 
bedridden
wheelchair=dd[, n:=sum(wheelchair=='Yes'), NetID][n>0, unique(NetID)]  # requiring at least 2 nurses in a network 

dialysis2=dd[, n:=sum(dialysis=='Yes'), NetID][n>1, unique(NetID)]
diabetes2=dd[, n:=sum(diabetes=='Yes'), NetID][n>1, unique(NetID)]
wound_care2=dd[, n:=sum(wound_care=='Yes'), NetID][n>1, unique(NetID)] 
ventilator2=dd[, n:=sum(ventilator=='Yes'), NetID][n>1, unique(NetID)] 
mdro2=dd[, n:=sum(mdro=='Yes'), NetID][n>1, unique(NetID)] 
wheelchair2=dd[, n:=sum(wheelchair=='Yes'), NetID][n>1, unique(NetID)]  # requiring at least 2 nurses in a network 
wheelchair2


## ---- eval=FALSE, echo=FALSE--------------------------------------------------
## # Example of a graph
## outliers <- c(4, 54, 57, 62 ) # networks with only 1 HCW
## set.seed(1231)
## net_figs <- vector("list", 4L)
## for (i in 1:4)
## j=outliers[i]
##   net_figs[[i]] <- nplot(
##     networks[[j]],
##     edge.curvature=0,
##     vertex.color = c(   hcw = "blue", residents = "green"
##       )[networks[[j]] %v% "group"]
##     )
## 
## # Putting them together using the gridExtra package
## gridExtra::grid.arrange(grobs=net_figs, nrow = 2, ncol = 2)


## ---- echo=FALSE--------------------------------------------------------------
log_n <- log(summary(bignet_mat ~ degrange(0, by = "net_id")))
# Summaries model parameters () 
termsum=summary(bignet_mat ~
    edges + 
    Sum("sum" ~ Sum(log_n ~ b1factor("net_id", levels=TRUE), "a"), "log_n" ) +
    b2star(2)+
    b2concurrent+
    b1star(2)+ #not in the model but.
    b2starmix(2, "type1")+ 
    b1dsp(1:9)+ # HCWs sharing 1~9 residents
    b2dsp(1:9)+ # residents sharing 1~9 HCWs
    b1starmix(2, "diabetes", diff = TRUE)+b2factor('diabetes')+
    b1starmix(2, "wound_care", diff = TRUE)+ b2factor('wound_care')+
    b1starmix(2, "dialysis", diff = TRUE) + b2factor('dialysis')+
    b1starmix(2, "wheelchair", diff = TRUE) +b2factor('wheelchair')+
    b1starmix(2, "bedridden", diff = TRUE)+b2factor('bedridden')+
    b1starmix(2, "ventilator", diff = TRUE) +b2factor('ventilator')+
    b1starmix(2, "mdro", diff = TRUE)+b2factor('mdro')
)

## ---- echo=FALSE--------------------------------------------------------------
dfrt=data.table(Ergm.Terms=names(termsum), Count=as.integer(termsum))
print(dfrt)
# b2factor:likely dialysis, mdro, bedridden, wound_care

## -----------------------------------------------------------------------------
default_control.1 <- list(
  seed            = 44485,
  main.method     = "Stochastic-Approximation", # Robbins-Monro Stochastic Approx.
  SA.nsubphases   = 4 * 1, # 1 the default
  SA.burnin       = 1024 * 4, # Twice the default
  SA.interval     = 1024 * 4 # four times... This is being passed to SA instead
)



## ---- cache = TRUE, message=FALSE, warning=FALSE, results='hide'--------------
res0 <- ergm(
  Networks(networks) ~
  N(~edges, ~log(n))+
  offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", c(default_control.1))
  )
#g=gofN(res0)
#plot(g)
#summary(g)["Pearson residuals"]
#summary(g)["Variance of Pearson residuals"] 


## ---- cache = TRUE, message=FALSE, warning=FALSE, results='hide'--------------
res1 <- ergm(
  Networks(networks) ~
  N(~edges, ~log(n))+ 
  N(~b2star(2))+   
  +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", c(default_control.1))
)
save(res1, file='res1_b2star2.RData')
g=gofN(res1)
summary(g)
names(summary(g))
summary(g)["variance of Pearson residuals"]


load('res2_c.RData')
default_control <- list( init= coef(res2.c),  seed    = 44485)

res2.cnew <- ergm(
  Networks(networks1) ~
  N(~edges, ~log(n))+
  N(~gwb1dsp(0.75, fixed=TRUE))+
  +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

#save(res2.c, file='res2_c.RData')
#load('res2_c.RData')
mcmc.diagnostics(res2.cnew)
g2=gofN(res2.cnew)
plot(g2)
save(res2.cnew, file='res2_cnew.RData')

summary(g)["Variance of Pearson residuals"]


coef(res2.c)
ini=c(3, -1, -1, 0, -Inf, -Inf)
default_control <- list( init= ini,  seed    = 44485)

# add b2factor
res5 <- ergm(
  Networks(networks) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2factor('dialysis'), subset=dialysis)+
    +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res5, file='res5_b2factordialysis.RData')



ini=c(3.2, -0.9, -1.2, -0.2, 0, -Inf, -Inf)
default_control <- list( init= ini,  seed    = 44485)

# add b2factor
res5.a <- ergm(
  Networks(networks) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)
    +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res5.a, file='res5_b2factordialysis_woundcare.RData')
mcmc.diagnostics(res5.a)
g=gofN(res5.a)
plot(g)

texreg::knitreg( list('gwb1dsp 0.75' = res2.c, res5, res5.a), single.row = TRUE, scalebox  = .8)

coef(res5.a)
ini=c(3, -0.9, -1.2, -0.2, 0.27, 0, -Inf, -Inf)
default_control <- list( init= ini,  seed    = 44485)

# add b2factor
res5.b <- ergm(
  Networks(networks) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('mdro'), subset=mdro)+
  +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res5.b, file='res5_b2factordialysis_woundcare_mdro.RData')


res5.c <- ergm(
  Networks(networks) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res5.c, file='res5_b2factordialysis_woundcare_bedridden.RData')


texreg::knitreg( list('gwb1dsp 0.75' = res2.c, res5, res5.a, res5.c), single.row = TRUE, scalebox  = .8)



# Now put b2star2 and gwb1dsp together
  #   N(1)~edges     N(log(n))~edges     b2star(2)  N(1)~gwb1dsp.fixed.0.75    offset(b1degree0)   offset(b2degree0)
ini=c(    3.1676520  ,     -0.8936725   ,      0,            -1.2297144      ,       -Inf       ,      -Inf )
default_control <- list( init= ini,  seed    = 44485)

res3 <- ergm(
  Networks(networks) ~
  N(~edges, ~log(n))+
  N(~b2star(2))+
  N(~gwb1dsp(1, fixed=TRUE))+
  +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

mcmc.diagnostics(res3)
g=gofN(res3)
summary(g)["Variance of Pearson residuals"]
plot(g)
load('res3.RData')
#save(res3, file='res3.RData')



texreg::knitreg( list('b2star2'=res1, 'gwb1dsp 0.75' = res2.c,  'Both'=res3), single.row = TRUE, scalebox  = .8)
# use decay 0.75 or 0.5.

# not bad, but not much better. Obviously correlated.

# Check interaction between gwb1dsp and log(n) - use res2.c as the base model
ini=c(    3,     -0.8936725   ,       -1.2297144   ,  0  ,    0,   -Inf       ,      -Inf )
default_control <- list( init= ini,  seed    = 44485)
# this will be the chosen res2.
 
res4 <- ergm(
  Networks(networks) ~
  N(~edges, ~log(n))+
  N(~gwb1dsp(0.75, fixed=TRUE))+
  N(~b2starmix(2, 'Nurse'), subset=Nurse )+ # has to be numeric
  +offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res4, file='res4_gwb1dsp_Nurse.RData')
load('res4_gwb1dsp_Nurse.RData')
mcmc.diagnostics(res4)
g=gofN(res4)
plot(g)
summary(g)["Variance of Pearson residuals"] # difficult to tell whether it has improved.


texreg::knitreg( list( res2.c, res3,  '+Nurse'=res4), single.row = TRUE, scalebox  = .8)
# 
coef(res5.c)
coef(res4)

ini=c(    3,     -0.8936725   ,       -1.2297144   ,  -0.06  ,    -0.08,  -0.22, 0.27,  0.15,  -Inf       ,      -Inf )
default_control <- list( init= ini,  seed    = 44485)
# this will be the chosen res2.
 
PTOTRT=as.numeric(PTOTRT)
Nurse=as.numeric(Nurse)

load('res5_b2factor_Nurse.RData')
default_control <- list( init= coef(res5),  seed    = 44485)
res5new <- ergm(
  Networks(networks1) ~
  N(~edges, ~log(n))+
  N(~gwb1dsp(0.75, fixed=TRUE))+
  N(~b2starmix(2, 'Nurse'), subset=Nurse )+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
  offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
save(res5new, file='res5_b2factor_Nurse_new.RData')
#save(res5, file='res5_b2factor_Nurse.RData')
par(mar = c(1, 1, 1, 1))
#load("res5_Nurse_PT.RData")
mcmc.diagnostics(res5new)
texreg::knitreg( list( res2.cnew, res5new), single.row = TRUE, scalebox  = .8)

g=gofN(res5new)
plot(g)


#load('res5_b2factor_Nurse.RData')

# Model 6, HCW type nurse, b2factors, b1starmix wound care

default_control <- list( init= coef(res6_new),  seed    = 44485)
#setdiff(wound_care, wound_care2)

res6_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care2)+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res6_new, file='res6_new.RData')
load('res6_new.RData')
load('res5_b2factor_Nurse_new.RData')
texreg::knitreg( list(res5new, res6_new), single.row = TRUE, scalebox  = .8)
g=gof(res6_new)
plot(g)


# Model 7, HCW type nurse, b2factors, b1starmix wound care, ventilator,
coef(res6_new)
#N(1)~edges                     N(log(n))~edges             N(1)~gwb1dsp.fixed.0.75          N(1)~b2starmix.2.Nurse.0.0 
ini=c(3.74088063      ,        -1.14954895      ,            -1.02778021         ,          -0.08828405 ,
#N(1)~b2starmix.2.Nurse.0.1          N(1)~b2factor.dialysis.Yes        N(1)~b2factor.wound_care.Yes         N(1)~b2factor.bedridden.Yes 
-0.25500002          ,     -0.21123648           ,   0.28677099             ,     0.20095538 ,
#N(1)~b1starmix.2.wound_care.HCW.No N(1)~b1starmix.2.wound_care.HCW.Yes                   offset(b1degree0)                   offset(b2degree0) 
0.03216966  ,           0.08562665          ,  0,0,         -Inf        ,              -Inf)

default_control <- list( init=ini,  seed    = 44485)

res7_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care2)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res7_new, file='res7_new.RData')
#load('res6_new.RData')
texreg::knitreg( list( res6_new, res7_new), single.row = TRUE, scalebox  = .8)
g=gof(res7_new)
plot(g)

# Model 8, HCW type nurse, b2factors, b1starmix wound care, ventilator, diabetes

coef(res7_new)
#N(1)~edges                     N(log(n))~edges             N(1)~gwb1dsp.fixed.0.75          N(1)~b2starmix.2.Nurse.0.0 
ini=c(3.58001627   ,     -1.10548015       ,      -1.02659913   ,      -0.11069885 ,
#N(1)~b2starmix.2.Nurse.0.1          N(1)~b2factor.dialysis.Yes        N(1)~b2factor.wound_care.Yes         N(1)~b2factor.bedridden.Yes 
-0.29830097     ,    -0.17553659       ,       0.30153264      ,      0.16151667 ,
#N(1)~b1starmix.2.wound_care.HCW.No N(1)~b1starmix.2.wound_care.HCW.Yes  N(1)~b1starmix.2.ventilator.HCW.No N(1)~b1starmix.2.ventilator.HCW.Yes 
0.03100716          ,   0.06761051    ,          0.01851044        ,       0.04759824 , 0,0,
#offset(b1degree0)                   offset(b2degree0) 
-Inf           ,             -Inf )
default_control <- list( init= ini,  seed    = 44485)


res8_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
load('res8_new.RData')
save(res8_new, file='res8_new.RData')

texreg::knitreg( list(  res6_new, res7_new, res8_new), single.row = TRUE, scalebox  = .8)
g=gofN(res8_new)

pdf('gofN_model8.pdf')
plot(g)
dev.off()
mcmc.diagnostics(res8_new)


#N(1)~edges                     N(log(n))~edges             N(1)~gwb1dsp.fixed.0.75          N(1)~b2starmix.2.Nurse.0.0 

ini=c(3.64107433      ,             -1.14903980        ,                 -0.96426243                    ,     -0.10751787 ,
#N(1)~b2starmix.2.Nurse.0.1          N(1)~b2factor.dialysis.Yes        N(1)~b2factor.wound_care.Yes         N(1)~b2factor.bedridden.Yes 
-0.29451332               ,          -0.20798950           ,               0.22684924              ,           0.17172971 ,
#N(1)~b1starmix.2.wound_care.HCW.No N(1)~b1starmix.2.wound_care.HCW.Yes  N(1)~b1starmix.2.ventilator.HCW.No N(1)~b1starmix.2.ventilator.HCW.Yes 
0.01705399              ,            0.04478927         ,                 0.01816528         ,                 0.04513037 ,
#N(1)~b1starmix.2.diabetes.HCW.No   N(1)~b1starmix.2.diabetes.HCW.Yes                   offset(b1degree0)                   offset(b2degree0) 
0.02583105                   ,       0.04475215           ,             0,    -Inf                        ,        -Inf )
default_control <- list( init= ini,  seed    = 44485)


res9_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2factor('mdro'), subset=mdro)+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
#load('res8_new.RData')
save(res9_new, file='res9_new.RData')

texreg::knitreg( list(res8_new, res9_new), single.row = TRUE, scalebox  = .8)
g=gofN(res8_new)
plot(g)
mcmc.diagnostics(res8_new)



res9_new2 <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'mdro', diff=TRUE), subset=mdro2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
#load('res8_new.RData')
save(res9_new, file='res9_new.RData')

texreg::knitreg( list(res8_new, res9_new), single.row = TRUE, scalebox  = .8)



ini=c(3.64107433      ,             -1.14903980        ,                 -0.96426243                    ,     -0.10751787 ,
#N(1)~b2starmix.2.Nurse.0.1          N(1)~b2factor.dialysis.Yes        N(1)~b2factor.wound_care.Yes         N(1)~b2factor.bedridden.Yes 
-0.29451332               ,          -0.20798950           ,               0.22684924              ,           0.17172971 ,
#N(1)~b1starmix.2.wound_care.HCW.No N(1)~b1starmix.2.wound_care.HCW.Yes  N(1)~b1starmix.2.ventilator.HCW.No N(1)~b1starmix.2.ventilator.HCW.Yes 
0.01705399              ,            0.04478927         ,                 0.01816528         ,                 0.04513037 ,
#N(1)~b1starmix.2.diabetes.HCW.No   N(1)~b1starmix.2.diabetes.HCW.Yes                   offset(b1degree0)                   offset(b2degree0) 
0.02583105                   ,       0.04475215           ,             0, 0,   -Inf                        ,        -Inf )
default_control <- list( init= ini,  seed    = 44485)


res10_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2starmix(2, 'PTOTRT'), subset=PTOTRT)+ # has to be numeric
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
#load('res8_new.RData')
save(res10_new, file='res10_new.RData')

texreg::knitreg( list(res8_new, res9_new, res10_new), single.row = TRUE, scalebox  = .8)
g=gofN(res10_new)
pdf('gofN_model10_new.pdf')
plot(g)
mcmc.diagnostics(res10_new)
dev.off()


# replace bedridden with wheelchair

default_control <- list( init= coef(res10_new),  seed    = 44485)

res11_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('wheelchair'), subset=wheelchair)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2starmix(2, 'PTOTRT'), subset=PTOTRT)+ # has to be numeric
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
#load('res8_new.RData')
save(res11_new, file='res11_new.RData')

texreg::knitreg( list(res8_new, res9_new, res10_new, res11_new), single.row = TRUE, scalebox  = .8)

ini=coef(res10_new)[-6]
default_control <- list( init=ini,  seed    = 44485)

res12_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    #N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2starmix(2, 'PTOTRT'), subset=PTOTRT)+ # has to be numeric
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
#load('res8_new.RData')
save(res12_new, file='res12_new.RData')

texreg::knitreg( list(res8_new, res9_new, res10_new, res11_new, res12_new), single.row = TRUE, scalebox  = .8)

ini=coef(res10_new)[-6]
default_control <- list( init=ini,  seed    = 44485)

coef(res12_new)
ini=coef(res12_new)[-c(8,9)]
default_control <- list( init=ini,  seed    = 44485)

res13_new <- ergm(
  Networks(networks1) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    #N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    #N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2starmix(2, 'PTOTRT'), subset=PTOTRT)+ # has to be numeric
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))
#load('res8_new.RData')

save(res13_new, file='res13_new.RData')

texreg::knitreg( list(res8_new, res9_new, res10_new, res11_new, res12_new, res13_new), single.row = TRUE, scalebox  = .8)



# Last, fit model using original 93 networks.

load('dd_nodal_attr.RData')
dim(dd)
setDT(dd)
dd$NetID=as.numeric(dd$NetID)

Nurse=dd[, n:=sum(Nurse==1), NetID][n>1, unique(NetID)] %>% as.numeric # requiring at least 2 nurses in a network 
CNA=dd[, n:=sum(CNA==1), NetID][n>1, unique(NetID)] %>% as.numeric
PTOTRT=dd[, n:=sum(PTOTRT==1), NetID][n>1, unique(NetID)] %>% as.numeric
Other=dd[, n:=sum(Other==1), NetID][n>1, unique(NetID)] 
nuni(c(Other, PTOTRT)) # 46

wh=intersect(Nurse, PTOTRT)
length(wh)
dcheck=dd[NetID %in% wh & group=='hcw', .(NetID, Nurse, PTOTRT, Other, CNA)]
dcheck[, count:=sum(CNA=="1")+sum(Other=="1"), NetID]
PTOTRT.nurse=as.numeric(dcheck[count>1, unique(NetID)])
length(PTOTRT.nurse) # all 30 networks are fine!

dialysis=dd[, n:=sum(dialysis=='Yes'), NetID][n>0, unique(NetID)]
wound_care=dd[, n:=sum(wound_care=='Yes'), NetID][n>0, unique(NetID)] 
mdro=dd[, n:=sum(mdro=='Yes'), NetID][n>0, unique(NetID)] 
bedridden=dd[, n:=sum(bedridden=='Yes'), NetID][n>0, unique(NetID)]  # requiring at least 2 nurses in a network 
bedridden

dialysis2=dd[, n:=sum(dialysis=='Yes'), NetID][n>1, unique(NetID)]
diabetes2=dd[, n:=sum(diabetes=='Yes'), NetID][n>1, unique(NetID)]
wound_care2=dd[, n:=sum(wound_care=='Yes'), NetID][n>1, unique(NetID)] 
ventilator2=dd[, n:=sum(ventilator=='Yes'), NetID][n>1, unique(NetID)] 
mdro2=dd[, n:=sum(mdro=='Yes'), NetID][n>1, unique(NetID)] 
wheelchair2=dd[, n:=sum(wheelchair=='Yes'), NetID][n>1, unique(NetID)]  # requiring at least 2 nurses in a network 
wheelchair2


default_control <- list( init= coef(res10_new),  seed    = 44485)

res10_93 <- ergm(
  Networks(networks) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2starmix(2, 'PTOTRT'), subset=PTOTRT)+ # has to be numeric
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

g=gofN(res10_93)

pdf('mcmc_diagnostics_and_gofN_model10_93.pdf')
mcmc.diagnostics(res10_93)
plot(g)
dev.off()

save(res10_93, file='res10_93.RData')

load('res12_new.RData')
coef(res12_new)
g=gofN(res12_new)

pdf('mcmc_diagnostics_and_gofN_model12_new.pdf')
mcmc.diagnostics(res12_new)
plot(g)
dev.off()



default_control <- list( init= coef(res12_new),  seed    = 44485)

res12_93 <- ergm(
  Networks(networks) ~
    N(~edges, ~log(n))+
    N(~gwb1dsp(0.75, fixed=TRUE))+
    N(~b2starmix(2, 'Nurse'), subset=Nurse)+ # has to be numeric
    #N(~b2factor('dialysis'), subset=dialysis)+
    N(~b2factor('wound_care'), subset=wound_care)+
    N(~b2factor('bedridden'), subset=bedridden)+
    N(~b1starmix(2, 'wound_care', diff=TRUE), subset=wound_care)+
    N(~b1starmix(2, 'ventilator', diff=TRUE), subset=ventilator2)+
    N(~b1starmix(2, 'diabetes', diff=TRUE), subset=diabetes2)+
    N(~b2starmix(2, 'PTOTRT'), subset=PTOTRT)+ # has to be numeric
    offset(b1degree(0))+offset(b2degree(0)),
  offset.coef = c(-Inf, -Inf),
  control = do.call("control.ergm", default_control))

save(res12_93, file='res12_93.RData')


load('res10_93.RData')
texreg::knitreg( list(res10_93,res12_93), single.row = TRUE, scalebox  = .8)
g=gofN(res12_93)

pdf('mcmc_diagnostics_and_gofN_model12_93.pdf')
mcmc.diagnostics(res12_93)
plot(g)
dev.off()


net=networks[[57]]
nplot( net, vertex.color=c(  hcw = "blue", residents = "green"
       )[net %v% "group"])

nplot( net, vertex.color=c(  HCW = "blue", Yes = "red", No="green"
       )[net %v% "wound_care"])
