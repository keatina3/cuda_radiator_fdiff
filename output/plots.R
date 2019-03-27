shared <- read.table("~/Documents/5615/Assignments/Assignment_2/output/GPU_sharedpitch_times.csv",header=T,sep=",")
glob <- read.table("~/Documents/5615/Assignments/Assignment_2/output/GPU_glob_times.csv",header=T,sep=",")
cpu <- read.table("~/Documents/5615/Assignments/Assignment_2/output/CPU_times.csv",header=T,sep=",")
double <- read.table("~/Documents/5615/Assignments/Assignment_2/output/GPU_doubles_times.csv", header=T, sep=",")

plotvals <- function(sharedVals, globVals, cpuVals, leg1, leg2){
    par(mfrow=c(1,1))
    shared = log10(rep(cpuVals, each = 6)/sharedVals)
    glob = log10(rep(cpuVals, each = 6)/globVals)
    
    blocks = log10(c(16,32,64,128,256,1024))
    logblocks2 = log10(c(16,32,64,128,256,512,1024))
    xlabs = c(16,32,64,128,256,512, 1024)
    
    ylims = c(floor(min(c(shared))),ceiling(max(c(shared))))
    ylabs = seq(min(ylims),max(ylims),by=1)
    ylabs2 = 10^ylabs
    ylabs2 = paste(ylabs2,"x")
    plot(blocks, shared[1:6], ylim = ylims, pch = 17,col='red',xaxt='n',yaxt='n', ylab="Speedup", xlab = "Block-size")
    axis(1, at = logblocks2, labels = xlabs)
    axis(2, at = ylabs, labels = ylabs2)
    abline(h=0,col="grey",lty=2, lwd=2)
    points(blocks, shared[7:12], pch = 17, col='blue')
    points(blocks, shared[13:18], pch = 17, col='green')
    points(blocks, shared[19:24], pch = 17, col='purple')
    
    lines(blocks, shared[1:6], lwd = 2, col='red')
    lines(blocks, shared[7:12], lwd = 2, col='blue')
    lines(blocks, shared[13:18], lwd = 2, col=3)
    lines(blocks, shared[19:24], lwd = 2, col='purple')
    
    cols=c('red','blue',3,'purple')
    legend("topleft",legend=c("100x100","1000x1000","5000x5000","15360x15360"),col=cols, lwd=2, pch=17)

    ylims = c(floor(min(c(glob))),ceiling(max(c(glob))))
    ylabs = seq(min(ylims),max(ylims),by=1)
    ylabs2 = 10^ylabs
    ylabs2 = paste(ylabs2,"x")
    plot(blocks, glob[1:6], ylim = ylims, pch = 17,col='red',xaxt='n', yaxt ='n',ylab="Speedup", xlab = "Block-size")
    axis(1, at = logblocks2, labels = xlabs)
    axis(2, at = ylabs, labels = ylabs2)
    abline(h=0,col="grey",lty=2, lwd=2)
    points(blocks, glob[7:12], pch = 17, col='blue')
    points(blocks, glob[13:18], pch = 17, col=3)
    points(blocks, glob[19:24], pch = 17, col='purple')
    
    lines(blocks, glob[1:6], lwd = 2, col='red')
    lines(blocks, glob[7:12], lwd = 2, col='blue')
    lines(blocks, glob[13:18], lwd = 2, col=3)
    lines(blocks, glob[19:24], lwd = 2, col='purple')
    legend("topleft",legend=c("100x100","1000x1000","5000x5000","15360x15360"),col=cols, lwd=2, pch=17)

    par(mfrow=c(2,2))
    for (i in c(1:4)){
      ylims = c(floor(min(c(glob[((i-1)*6):(i*6)],shared[((i-1)*6):(i*6)]))),
                        ceiling(max(c(glob[((i-1)*6):(i*6)],shared[((i-1)*6):(i*6)]))))
      plot(blocks, shared[(1+(i-1)*6):(i*6)], ylim = ylims, pch = 17,col=cols[i],xaxt='n', 
                    yaxt='n',ylab="Speedup", xlab = "Block-size",main = paste(cpu$NxM[i]))
      abline(h=0,col="grey",lty=2, lwd=2)
      ylabs = seq(min(ylims),max(ylims),by=1)
      ylabs2 = 10^ylabs
      ylabs2 = paste(ylabs2,"x")
      axis(1, at = logblocks2, labels = xlabs)
      axis(2, at = ylabs, labels = ylabs2)
      lines(blocks, shared[(1+(i-1)*6):(i*6)], lwd = 2, col=cols[i])
      points(blocks, glob[(1+(i-1)*6):(i*6)], pch = 17, col=cols[i])
      lines(blocks, glob[(1+(i-1)*6):(i*6)], lwd = 2, lty=2, col=cols[i])
      legend("topleft",legend=c(paste(leg1),paste(leg2)),col=cols[i], lwd=2,lty=c(1,2),cex=0.75)
    }
}

shareVglob <- function(sharedVals, globVals, log=FALSE){
  par(mfrow=c(1,1))
  if(log)
    shared = log10(globVals/sharedVals)
  else
    shared = globVals/sharedVals
  
  blocks = log10(c(16,32,64,128,256,1024))
  logblocks2 = log10(c(16,32,64,128,256,512,1024))
  xlabs = c(16,32,64,128,256,512, 1024)
  
  if(!log){
    ylims = c(floor(min(c(shared))),ceiling(max(c(shared))))
    ylabs = round(c(seq(min(ylims),max(ylims),length.out = 5)),4)
    ylabs2 = paste(ylabs,"x")
  } else {
    ylims = c(floor(min(c(shared))),ceiling(max(c(shared))))
    ylabs = seq(min(ylims),max(ylims),by=1)
    ylabs2 = 10^ylabs
    ylabs2 = paste(ylabs2,"x")
  }
  plot(blocks, shared[1:6], ylim = ylims, pch = 17,col='red',xaxt='n', yaxt = 'n', ylab="Speedup", xlab = "Block-size")
  axis(1, at = logblocks2, labels = xlabs)
  axis(2, at = ylabs, labels = ylabs2)
  abline(h=0,col="grey",lty=2, lwd=2)
  if(!log) abline(h=1,col="grey",lty=2, lwd=2)
  points(blocks, shared[7:12], pch = 17, col='blue')
  points(blocks, shared[13:18], pch = 17, col='green')
  points(blocks, shared[19:24], pch = 17, col='purple')
  
  lines(blocks, shared[1:6], lwd = 2, col='red')
  lines(blocks, shared[7:12], lwd = 2, col='blue')
  lines(blocks, shared[13:18], lwd = 2, col=3)
  lines(blocks, shared[19:24], lwd = 2, col='purple')
  cols=c('red','blue',3,'purple')
  legend("topleft",legend=c("100x100","1000x1000","5000x5000","15360x15360"),col=cols, lwd=2, pch=17)
}

plotvals(shared$FD.Calc,glob$FD.Calc, cpu$FD.Calc/1000,"shared","glob")
plotvals(shared$Avg.Calc,glob$Avg.Calc, cpu$Avg.Calc/1000,"shared","glob")
plotvals(shared$Total,glob$Total, cpu$Total,"shared","glob")

shareVglob(shared$Alloc,glob$Alloc,1)
shareVglob(shared$Transfer.RAM,glob$Transfer.RAM,1)

plotvals(shared$Total,double$Total, cpu$Total,"float","double")
shareVglob(shared$Total,double$Total,1)

max(rep(cpu$Total,each=6)/shared$Total)