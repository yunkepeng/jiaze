pbrblue <- colorRampPalette(brewer.pal(11, "Blues"), interpolate = "spline")
B_palette<-c(pbrblue(11),"black")

# Set the boundary 
boundary_0 <- read.csv("/Users/yunpeng/Desktop/Examples/Data/S_Quercus_robur.csv")
boundary_0<-as.data.frame(na.omit(boundary_0))


# Input the derived GAM results 
mtrx.melt <- read.csv("/Users/yunpeng/Desktop/Examples/Data/P_Quercus_robur.csv")


for (i in 1:5){ 
  boundary1<-boundary_0[boundary_0$gdd>0 & boundary_0$gdd<1000,] 
  x1<-as.matrix(boundary1[,c(1,2)])
  y1=base::unique(x1)
  ahull.t1<- ahull(y1,alpha=100)
  #colnames(mtrx.melt1)=colnames(mtrx.melt2)=colnames(mtrx.melt3)=colnames(mtrx.melt4)=colnames(mtrx.melt5)=c("rtmi","mtco","Taxon")
  detect1<-as.vector(inahull(ahull.t1,p=c(mtrx.melt1$rtmi,mtrx.melt1$mtco)))
  mtrx.melt.a1<-cbind(mtrx.melt1,detect1)
  mtrx.melt.a1<-mtrx.melt.a1[!mtrx.melt.a1$detect1 == FALSE,] 
  #max_pct=max(mtrx.melt.a1$Taxon,mtrx.melt.a2$Taxon,mtrx.melt.a3$Taxon,mtrx.melt.a4$Taxon,mtrx.melt.a5$Taxon)
  #and select something....
  mtrx.melt.a1$pct <- cut(mtrx.melt.a1$Taxon, cuttime)
  
  mtrx_to_dcast1<-mtrx.melt.a1[,1:3]
  mtrx.dcast1<-dcast(mtrx_to_dcast1,rtmi~mtco,value.var="Taxon") 
  rownames(mtrx.dcast1)=mtrx.dcast1$rtmi 
  mtrx.dcast1<-mtrx.dcast1[,-1] 
  p1<-filled.contour(
    x = sort(unique(mtrx_to_dcast1$rtmi)), 
    y = sort(unique(mtrx_to_dcast1$mtco)),
    z = as.matrix(mtrx.dcast1), 
    levels = pretty(c(0,m_pct), 10), 
    xlim= c(0,3), ylim= c(-46.10,26.90), 
    col= c_palette, 
    plot.title = title(main = "GDD5 = 500",cex.main=1, 
                       xlab = "sqrt(Moisture Index)", ylab = "MTCO"), 
    key.title = title(main="Presence 
  Probability",cex.main=0.8))
  
  }



boundary1<-boundary_0[boundary_0$gdd>0 & boundary_0$gdd<1000,] 
boundary2<-boundary_0[boundary_0$gdd>1000 & boundary_0$gdd<2000,] 
boundary3<-boundary_0[boundary_0$gdd>2000 & boundary_0$gdd<3000,] 
boundary4<-boundary_0[boundary_0$gdd>3000 & boundary_0$gdd<4000,] 
boundary5<-boundary_0[boundary_0$gdd>4000,] 

# Set the alphahull
x1<-as.matrix(boundary1[,c(1,2)])
x2<-as.matrix(boundary2[,c(1,2)])
x3<-as.matrix(boundary3[,c(1,2)])
x4<-as.matrix(boundary4[,c(1,2)])
x5<-as.matrix(boundary5[,c(1,2)])

y1=base::unique(x1)
y2=base::unique(x2) 
y3=base::unique(x3)
y4=base::unique(x4)
y5=base::unique(x5)

ahull.t1<- ahull(y1,alpha=100)
ahull.t2<- ahull(y2,alpha=100) 
ahull.t3<- ahull(y3,alpha=100)
ahull.t4<- ahull(y4,alpha=100)
ahull.t5<- ahull(y5,alpha=100)

mtrx.melt1 <- mtrx.melt[mtrx.melt$gdd==500,c(1,2,4)] # extract values when gdd=500 
mtrx.melt2 <- mtrx.melt[mtrx.melt$gdd==1500,c(1,2,4)] # extract values when gdd=1500 
mtrx.melt3 <- mtrx.melt[mtrx.melt$gdd==2500,c(1,2,4)] # extract values when gdd=2500
mtrx.melt4 <- mtrx.melt[mtrx.melt$gdd==3500,c(1,2,4)] # extract values when gdd=3500
mtrx.melt5 <- mtrx.melt[mtrx.melt$gdd==4500,c(1,2,4)] # extract values when gdd=4500

colnames(mtrx.melt1)=colnames(mtrx.melt2)=colnames(mtrx.melt3)=colnames(mtrx.melt4)=colnames(mtrx.melt5)=c("rtmi","mtco","Taxon")

# Use Alphahull to remove outside points
detect1<-as.vector(inahull(ahull.t1,p=c(mtrx.melt1$rtmi,mtrx.melt1$mtco)))
detect2<-as.vector(inahull(ahull.t2,p=c(mtrx.melt2$rtmi,mtrx.melt2$mtco)))
detect3<-as.vector(inahull(ahull.t3,p=c(mtrx.melt3$rtmi,mtrx.melt3$mtco)))
detect4<-as.vector(inahull(ahull.t4,p=c(mtrx.melt4$rtmi,mtrx.melt4$mtco)))
detect5<-as.vector(inahull(ahull.t5,p=c(mtrx.melt5$rtmi,mtrx.melt5$mtco)))

mtrx.melt.a1<-cbind(mtrx.melt1,detect1)
mtrx.melt.a2<-cbind(mtrx.melt2,detect2)
mtrx.melt.a3<-cbind(mtrx.melt3,detect3)
mtrx.melt.a4<-cbind(mtrx.melt4,detect4)
mtrx.melt.a5<-cbind(mtrx.melt5,detect5)

mtrx.melt.a1<-mtrx.melt.a1[!mtrx.melt.a1$detect1 == FALSE,] 
mtrx.melt.a2<-mtrx.melt.a2[!mtrx.melt.a2$detect2 == FALSE,] 
mtrx.melt.a3<-mtrx.melt.a3[!mtrx.melt.a3$detect3 == FALSE,]
mtrx.melt.a4<-mtrx.melt.a4[!mtrx.melt.a4$detect4 == FALSE,]
mtrx.melt.a5<-mtrx.melt.a5[!mtrx.melt.a5$detect5 == FALSE,]

max_pct=max(mtrx.melt.a1$Taxon,mtrx.melt.a2$Taxon,mtrx.melt.a3$Taxon,mtrx.melt.a4$Taxon,mtrx.melt.a5$Taxon)

# Set the cuttime for different color palette
if (max_pct>=0.9){ 
  m_pct=as.integer((max_pct+0.1)*10)/10 
  c_palette=B_palette 
} else if (max_pct>=0.1){
  m_pct=as.integer((max_pct+0.1)*10)/10 
  c_palette=B_palette 
} else if (max_pct>=0.01){
  m_pct=as.integer((max_pct+0.01)*100)/100 
  c_palette=B_palette 
}else { 
  m_pct=as.integer((max_pct+0.001)*1000)/1000 
  c_palette=B_palette } 
cuttime=c(seq(0,m_pct,m_pct/10))

# Use the cuttime 
mtrx.melt.a1$pct <- cut(mtrx.melt.a1$Taxon, cuttime)
mtrx.melt.a2$pct <- cut(mtrx.melt.a2$Taxon, cuttime)
mtrx.melt.a3$pct <- cut(mtrx.melt.a3$Taxon, cuttime)
mtrx.melt.a4$pct <- cut(mtrx.melt.a4$Taxon, cuttime)
mtrx.melt.a5$pct <- cut(mtrx.melt.a5$Taxon, cuttime)

# 2D Plot - Without Point
#pdf(file = "~/Desktop/Examples/Plots/2DP_Quercus_robur.pdf", width = 6,height = 5) 
## Plot when gdd = 500
mtrx_to_dcast1<-mtrx.melt.a1[,1:3]
mtrx.dcast1<-dcast(mtrx_to_dcast1,rtmi~mtco,value.var="Taxon") 
rownames(mtrx.dcast1)=mtrx.dcast1$rtmi 
mtrx.dcast1<-mtrx.dcast1[,-1] 
p1<-filled.contour(
  x = sort(unique(mtrx_to_dcast1$rtmi)), 
  y = sort(unique(mtrx_to_dcast1$mtco)),
  z = as.matrix(mtrx.dcast1), 
  levels = pretty(c(0,m_pct), 10), 
  xlim= c(0,3), ylim= c(-46.10,26.90), 
  col= c_palette, 
  plot.title = title(main = "GDD5 = 500",cex.main=1, 
                     xlab = "sqrt(Moisture Index)", ylab = "MTCO"), 
  key.title = title(main="Presence 
  Probability",cex.main=0.8))