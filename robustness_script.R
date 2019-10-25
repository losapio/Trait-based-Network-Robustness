# Resistance of plant–plant networks to biodiversity loss and secondary extinctions following simulated environmental changes

# Gianalberto Losapio* & Christian Schöb
# Department of Evolutionary Biology and Environmental Studies, University of Zurich, Winterthurerstrasse 190, CH-8057 Zurich, Switzerland
# Correspondence author: E-mail: gianalberto.losapio@ieu.uzh.ch

# Functional Ecology, doi: 10.1111/1365-2435.12839

library(FactoMineR) # PCA analysis
library(igraph) # network analysis
library(rnetcarto) # network modularity
library(asreml) # mixed effects models
library(pascal) # Wald test
library(nnet) # multinomial logit model
library(car) # ANOVA tests

# data frame
                 cushion             species   LDMC   SLA leafarea area
2 Arenaria_tetraquetra_F Agrostis_nevadensis 321.39  8.90     9.90   18
3 Arenaria_tetraquetra_F  Lotus_corniculatus 250.00 14.89    17.57   18
6 Arenaria_tetraquetra_F Agrostis_nevadensis 303.49  8.13     8.49   19
7 Arenaria_tetraquetra_F Agrostis_nevadensis 267.15 10.75    11.81   19
8 Arenaria_tetraquetra_F  Lotus_corniculatus 226.42 16.29     7.82   19
9 Arenaria_tetraquetra_F  Lotus_corniculatus 214.81 18.01     9.40   19
  individual sample    LA
2          1      1  9.90
3          1      1 17.57
6          1      1  8.49
7          2      1 11.81
8          1      1  7.82
9          2      1  9.40

# build unipartite graph
g <- graph.data.frame(data11, directed=FALSE)

# compute guimera modularity
guimod <- netcarto(g.bip1)

# trait distribution parameters (mean and sd) for each species
ldmc.mean 	<- tapply(data112$ldmc, data112$species, mean, na.rm = TRUE)
sla.mean	<- tapply(data112$sla, data112$species, mean, na.rm = TRUE)
la.mean 	<- tapply(data112$la, data112$species, mean, na.rm = TRUE)
ldmc.sd 	<- tapply(data112$ldmc, data112$species, sd, na.rm = TRUE)
sla.sd 	<- tapply(data112$sla, data112$species, sd, na.rm = TRUE)
la.sd 	<- tapply(data112$la, data112$species, sd, na.rm = TRUE)

# create vector of trait values for removal
# with 25 fixed interval value for each trait
q=1:25

# ldmc
delta.ldmc = (max(ldmc.mean)-(min(ldmc.mean)))/(length(q)-1)
q.ldmc = q
for (i in 0:length(q.ldmc)){q.ldmc[i] <- min(ldmc.mean)+((i-1)*delta.ldmc)}
q.ldmc <- round(q.ldmc,3)
# sla
delta.sla = (max(sla.mean)-(min(sla.mean)))/(length(q)-1)
q.sla = q
for (i in 0:length(q.sla)){q.sla[i] <- min(sla.mean)+((i-1)*delta.sla)}
q.sla<- round(q.sla,3)
# la
delta.la = (max(la.mean)-(min(la.mean)))/(length(q)-1)
q.la = q
for (i in 0:length(q.la)){q.la[i] <- min(la.mean)+((i-1)*delta.la)}
q.la<- round(q.la,3)

# create matrix of trait values randomly sampled [function: rnorm]
# n. replicates (random samples) = 200
repl = 200

ldmc <- matrix(NA, repl, nrow(ft))
colnames(ldmc) <- rownames(ft)
for (i in 1:nrow(ft)){
	ldmc[,i] <- rnorm(repl, ft$ldmc.mean[i], ft$ldmc.sd[i])}

sla <- matrix(NA, repl, nrow(ft))
colnames(sla) <- rownames(ft)
for (i in 1:nrow(ft)){
	sla[,i] <- rnorm(repl, ft$sla.mean[i], ft$sla.sd[i])}

laf <- matrix(NA, repl, nrow(ft))
colnames(laf) <- rownames(ft)
for (i in 1:nrow(ft)){
	laf[,i] <- rnorm(repl, ft$la.mean[i], ft$la.sd[i])}

### pca for ldmc & sla
# create vector of removal
pca02 <- PCA(ft[,c("ldmc.mean","sla.mean")],graph = F) #axis 1 = high SLA e low LDMC
pca02$var$cor
# scores per specie
pca2 <- -round(pca02$ind$coord[,1],3)
delta.pca2 = (max(pca2)-(min(pca2)))/(length(q)-1)
q.pca2 = q
for (i in 0:length(q.pca2)){q.pca2[i] <- min(pca2)+((i-1)*delta.pca2)}
q.pca2 <- round(q.pca2,3)

# create matrix for pca
pca.02 <- matrix(NA, repl, nrow(ft))
colnames(pca.02) <- rownames(ft)

# order: 
for(i in 1:repl){
	pca02 <- PCA(cbind(ldmc[i,],sla[i,]), graph = F)
	ax <- ifelse(sign(pca02$var$cor[1,1])==sign(pca02$var$cor[2,1]), 2, 1) # select right component (i.e. axis)
	pca.02[i,] <- pca02$ind$coord[,ax]
	if(pca02$var$cor[1,ax]<0){pca.02[i,] <- -pca.02[i,]} # higher ldmc e lower sla
	pve.2[i] <- pca02$eig[ax,2]
}

### function for random extinction model

rm_random_nodes <- function(network, frac, no.rm){ #no.rm= id vertex to not remove (i.e."Open")
  if (!is.igraph(network)) {
    stop("Not a graph object")
  }
  
a <- network-V(network)[no.rm]
  
rem = frac*vcount(network)
if (rem>=vcount(a)){rem=vcount(a)}

network <- delete.vertices(network, sample(V(a)$name, rem))
}

### extinction sequences

for(i in 1:length(q)){ # number of steps
for(z in 1:repl){ # number of replicates
# e.g. scenario 1: extinction of species from low LA
	rem <- names(which(laf[z,]<=q.la[i])) # species to remove
	if (length(rem)==0){g1<-g} # nothing to remove
	else{g1 <- delete.vertices(g,rem)} # new network
# for now open is assumed always in the network # later we check if not
# survival
	task$network.size[z+(i-1)*repl] <- max(clusters(g1)$csize)-1   # -1 bc "open"
# absolute primary ext
	task$prim.ext[z+(i-1)*repl] <- vcount(g)-vcount(g1)
# absolute secondary ext
	task$sec.ext[z+(i-1)*repl] <- sum(clusters(g1)$csize==1)
# eigenvector centrality
	task$centr[z+(i-1)*repl] <- centralization.evcent(g1)$centralization
# check if open node remain isolated in the network
	if (length(neighbors(g1, "Open"))==0){
	# if other sp present and networked
	if(clusters(g1)$no!=1){
	task$prim.ext[z+(i-1)*repl] <- task$prim.ext[z+(i-1)*repl]+1	
	task$sec.ext[z+(i-1)*repl] = task$sec.ext[z+(i-1)*repl] -1} # 	
	# or if the last one ("desert")
	else{ 
	task$network.size[z+(i-1)*repl] = 0
	task$prim.ext[z+(i-1)*repl] = vcount(g)-1
	task$sec.ext[z+(i-1)*repl] = 0}}
# presence of species
for(k in 1:length(spft)){
	if(!is.na(match(spft[k], V(g1)$name))){
	task[z+(i-1)*repl,nc+k] <- length(neighbors(g1, spft[k])) # total abundance of sp k
	task[z+(i-1)*repl,nc+k] <- ifelse(task[z+(i-1)*repl,nc+k]==0,-1,1)} # presence or sec ext
	else{task[z+(i-1)*repl,nc+k] = 0}}
# check if foundation species remain isolated in the network
if (task[z+(i-1)*repl, "A.tetraquetra_F"]==-1){
	task$network.size[z+(i-1)*repl] <- task$network.size[z+(i-1)*repl]+1
	task$sec.ext[z+(i-1)*repl] <- task$sec.ext[z+(i-1)*repl] -1}
if (task[z+(i-1)*repl, "F.indigesta_F"]==-1){
	task$network.size[z+(i-1)*repl] <- task$network.size[z+(i-1)*repl]+1
	task$sec.ext[z+(i-1)*repl] <- task$sec.ext[z+(i-1)*repl] -1}
if (task[z+(i-1)*repl, "P.holosteum_F"]==-1){
	task$network.size[z+(i-1)*repl] <- task$network.size[z+(i-1)*repl]+1
	task$sec.ext[z+(i-1)*repl] <- task$sec.ext[z+(i-1)*repl] -1}


#scenario 5: RANDOM
	g1 <- rm_random_nodes(g, q1[i], 3) # new network
# for now open is assumed always in the network # later we check if not
# survival
	task$network.size[z+(i-1)*repl+4*m] <- max(clusters(g1)$csize)-1   # -1 bc "open"
# absolute primary ext
	task$prim.ext[z+(i-1)*repl+4*m] <- vcount(g)-vcount(g1)
# absolute secondary ext
	task$sec.ext[z+(i-1)*repl+4*m] <- sum(clusters(g1)$csize==1)
# eigenvector centrality
	task$centr[z+(i-1)*repl+4*m] <- centralization.evcent(g1)$centralization
# check if open node remain isolated in the network
	if (length(neighbors(g1, "Open"))==0){
	# if other sp present and networked
	if(clusters(g1)$no!=1){
	task$prim.ext[z+(i-1)*repl+4*m] <- task$prim.ext[z+(i-1)*repl+4*m]+1	
	task$sec.ext[z+(i-1)*repl+4*m] = task$sec.ext[z+(i-1)*repl+4*m] -1} # 	
	# or if the last one ("desert")
	else{ 
	task$network.size[z+(i-1)*repl+4*m] = 0
	task$prim.ext[z+(i-1)*repl+4*m] = vcount(g)-1
	task$sec.ext[z+(i-1)*repl+4*m] = 0}}
# presence of species
for(k in 1:length(spft)){
	if(!is.na(match(spft[k], V(g1)$name))){
	task[z+(i-1)*repl+4*m,nc+k] <- length(neighbors(g1, spft[k])) # total abundance of sp k
	task[z+(i-1)*repl+4*m,nc+k] <- ifelse(task[z+(i-1)*repl+4*m,nc+k]==0,-1,1)} # presence or sec ext
	else{task[z+(i-1)*repl+4*m,nc+k] = 0}}
# check if foundation species remain isolated in the network
if (task[z+(i-1)*repl+4*m, "A.tetraquetra_F"]==-1){
	task$network.size[z+(i-1)*repl+4*m] <- task$network.size[z+(i-1)*repl+4*m]+1
	task$sec.ext[z+(i-1)*repl+4*m] <- task$sec.ext[z+(i-1)*repl+4*m] -1}
if (task[z+(i-1)*repl+4*m, "F.indigesta_F"]==-1){
	task$network.size[z+(i-1)*repl+4*m] <- task$network.size[z+(i-1)*repl+4*m]+1
	task$sec.ext[z+(i-1)*repl+4*m] <- task$sec.ext[z+(i-1)*repl+4*m] -1}
if (task[z+(i-1)*repl+4*m, "P.holosteum_F"]==-1){
	task$network.size[z+(i-1)*repl+4*m] <- task$network.size[z+(i-1)*repl+4*m]+1
	task$sec.ext[z+(i-1)*repl+4*m] <- task$sec.ext[z+(i-1)*repl+4*m] -1}

# control for species associations, i.e. to which module subordinate species belong to
task$A.nevadensis[which(task$A.tetraquetra_F==0 & task$A.nevadensis==1)] <- -1
task$A.vulneraria[which(task$P.holosteum_F==0 & task$A.vulneraria==1)] <- -1
task$A.tetraquetra[which(task$A.tetraquetra_F==0 & task$A.tetraquetra==1)] <- -1
task$C.oporinoides[which(task$P.holosteum_F==0 & task$C.oporinoides==1)] <- -1
task$D.brachyanthus[which(task$P.holosteum_F==0 & task$D.brachyanthus==1)] <- -1
task$E.glaciale[which(task$A.tetraquetra_F==0 & task$E.glaciale==1)] <- -1
task$E.nevadense[which(task$A.tetraquetra_F==0 & task$E.nevadense==1)] <- -1
task$E.willkommii[which(task$A.tetraquetra_F==0 & task$E.willkommii==1)] <- -1
task$F.indigesta[which(task$A.tetraquetra_F==0 & task$F.indigesta==1)] <- -1
task$G.pyrenaicum[which(task$P.holosteum_F==0 & task$G.pyrenaicum==1)] <- -1
task$J.amethystina[which(task$F.indigesta_F==0 & task$J.amethystina==1)] <- -1
task$J.humilis[which(task$P.holosteum_F==0 & task$J.humilis==1)] <- -1
task$K.vallesiana[which(task$A.tetraquetra_F==0 & task$K.vallesiana==1)] <- -1
task$L.boryi[which(task$A.tetraquetra_F==0 & task$L.boryi==1)] <- -1
task$L.pectinata[which(task$A.tetraquetra_F==0 & task$L.pectinata==1)] <- -1
task$L.corniculatus[which(task$P.holosteum_F==0 & task$L.corniculatus==1)] <- -1
task$N.purpurea[which(task$P.holosteum_F==0 & task$N.purpurea==1)] <- -1
task$P.holosteum[which(task$A.tetraquetra_F==0 & task$P.holosteum==1)] <- -1
task$R.angiocarpus[which(task$P.holosteum_F==0 & task$R.angiocarpus==1)] <- -1
task$S.boissieri[which(task$P.holosteum_F==0 & task$S.boissieri==1)] <- -1
task$S.boryi[which(task$F.indigesta_F==0 & task$S.boryi==1)] <- -1
task$T.serpylloides[which(task$A.tetraquetra_F==0 & task$T.serpylloides)] <- -1

### correct for sec.ext and network.size
task$network.size[i] <- sum(task[i,c((nc+1):(nc+1+length(spft)-1))]==1)
task$sec.ext[i] <- sum(task[i,c((nc+1):(nc+1+length(spft)-1))]==-1)

### calculating:
# relative netework size = relative survival
task$rel.net.size = task$network.size/(vcount(g)-1)
# relative secondary extinctions to total extinctions
task$rel.sec.ext2 = task$sec.ext/(task$prim.ext+task$sec.ext)

### statystical analysis
# mixed effects
SURV <- task$rel.net.size
SECEXT <- task$rel.sec.ext2
EXTYPE <- factor(task$scenario, labels=c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Random"))
NCUSH <- factor(task$cush, ordered=TRUE)
Q <- factor(task$q, ordered=TRUE)
REPL <- factor(task$repl)

reml=asreml(fixed=SURV~EXTYPE*NCUSH, #treatment model
             random=~Q/REPL,
		 rcov=~id(EXTYPE):ar1(Q):REPL,keep.order=TRUE, #error model
             data=task,
             control=asreml.control(maxiter=50))
test.asreml(reml)

reml=asreml(fixed=SECEXT~EXTYPE*NCUSH, #treatment model
             random=~Q/REPL,
		 rcov=~id(EXTYPE):ar1(Q):REPL,keep.order=TRUE, #error model
             data=task,
             control=asreml.control(maxiter=50))
test.asreml(reml)

### for sec ext species, differences among species
# multivariate logit model

sp2.nul <- multinom(pers~1, data=task2)
sp2.1 <- multinom(pers~species-1, data=task2)
anova(sp2.nul,sp2.1)
Anova(sp2.1)
coef <- coefficients(sp2.1);coef
ci   <- confint(sp2.1)

# testing for species abudance and community membership
task.sp2$spab <- c(98,13,28,4,2,3,34,14,33,18,260,2,48,62,8,249,16,26,249,23,106,3,7,27,7)
task.sp2$secext <- tab.coef2[,4]
task.sp2$surv <- tab.coef2[,1]

spmod2 <- glm(secext~comm*spab, data=task.sp2)
summary(spmod2)
Anova(spmod2)



