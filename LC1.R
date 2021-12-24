pbrblue <- colorRampPalette(brewer.pal(11, "Blues"), interpolate = "spline")
B_palette<-c(pbrblue(11),"black")

# Set the boundary 
boundary_0 <- read.csv("/Users/yunpeng/Desktop/Examples/Data/S_Quercus_robur.csv")
boundary_0<-as.data.frame(na.omit(boundary_0))

# Input the derived GAM results 
mtrx.melt <- read.csv("/Users/yunpeng/Desktop/Examples/Data/P_Quercus_robur.csv")

#assign(substr("boundary",i,sep="/"), df2) 

for (i in 1:5){ 
  if (i<=4) { 
    ahull.t1<- ahull(base::unique(as.matrix((boundary_0[boundary_0$gdd>(1000*i-1000) & boundary_0$gdd<(1000*i),])[,c(1,2)])),alpha=100)
  } else {
    ahull.t1<- ahull(base::unique(as.matrix((boundary_0[boundary_0$gdd>(1000*i-1000),])[,c(1,2)])),alpha=100)
  }
  
  mtrx.melt1 <- mtrx.melt[mtrx.melt$gdd==1000*i-500,c(1,2,4)] # extract values when gdd=500 
  detect1<-as.vector(inahull(ahull.t1,p=c(mtrx.melt1$rtmi,mtrx.melt1$mtco)))
  
  mtrx.melt.a1<-cbind(mtrx.melt1,detect1)
  mtrx.melt.a1<-mtrx.melt.a1[!mtrx.melt.a1$detect1 == FALSE,] 
  assign(paste("mtrx.melt.a",i,sep=""), mtrx.melt.a1) 
  #max_pct=max(mtrx.melt.a1$Taxon,mtrx.melt.a2$Taxon,mtrx.melt.a3$Taxon,mtrx.melt.a4$Taxon,mtrx.melt.a5$Taxon)
  #and select something....
}

#select manually
max_pct=max(mtrx.melt.a1$Species,mtrx.melt.a2$Species,mtrx.melt.a3$Species,mtrx.melt.a4$Species,mtrx.melt.a5$Species)
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

pdf(file = paste("/Users/yunpeng/Desktop/Examples/Plots/2DP_Quercus_robur.pdf",sep=""), width = 6, height = 5) 

for (i in 1:5){ 
  mtrx.melt.a1 <- eval(parse(text=paste("mtrx.melt.a",i,sep="")))
  mtrx.melt.a1$pct <- cut(mtrx.melt.a1$Species, cuttime)
  # Plot when gdd = 500
  mtrx_to_dcast1<-mtrx.melt.a1[,1:3]
  mtrx.dcast1<-dcast(mtrx_to_dcast1,rtmi~mtco,value.var="Species") 
  rownames(mtrx.dcast1)=mtrx.dcast1$rtmi 
  mtrx.dcast1<-mtrx.dcast1[,-1]

  filled.contour(
    x = sort(unique(mtrx_to_dcast1$rtmi)), 
    y = sort(unique(mtrx_to_dcast1$mtco)),
    z = as.matrix(mtrx.dcast1), 
    levels = pretty(c(0,m_pct), 10), 
    xlim= c(0,3), ylim= c(-46.10,26.90), 
    col= c_palette, 
    plot.title = title(main = paste("GDD5 =",i*1000-500,sep=" "),cex.main=1, 
                       xlab = "sqrt(Moisture Index)", ylab = "MTCO"), 
    key.title = title(main="Presence 
  Probability",cex.main=0.8))
}
dev.off()

