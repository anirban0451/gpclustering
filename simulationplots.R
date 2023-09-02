load("origvsvecchia05050503.Rdata")

replicates = 25
misclasfull = array(dim = c(2, replicates))
mvecchia = c(floor(p/20), floor(p/12), floor(p/10), floor(p/5))
for(i in 1:nrow(l_f_matrix)){
  for(j in 1:replicates){
    
    misclasfull[i, j] = NMI(clusterdata_full[i, j, ],
                            rep(1:2, each = nsamp))
  }
}

misclasvecchia = array(dim = c(nrow(l_f_matrix), replicates, length(mvecchia)))

for(i in 1:nrow(l_f_matrix)){
  
  for(j in 1:replicates){
    
    for(k in 1:length(mvecchia)){
      
      misclasvecchia[i, j, length(mvecchia) - k + 1] = 
        NMI(clusterdata_veccs[i, j, k, ],
            rep(1:2, each = nsamp))
    }
  }
}

groups = as.factor(rep(1:5, each = replicates))

y = c(c(misclasfull[1, ]), c( misclasvecchia[1, , ]))

label = c("Full EM",
          "VEM (m = 60)",
          "VEM (m = 30)",
          "VEM (m = 25)",
          "VEM (m = 15)")
pdf("boxplot2.pdf", height = 4, width = 8)
par(mar = c(1, 3, 1, 1), oma = c(1, 0, 0, 0))
boxplot(y ~ groups, ylab = "", names = label, 
        col = RColorBrewer::brewer.pal(length(label), "Set2"))
mtext("y", side = 2, line = 2)
dev.off()

clustertimefullavg = median(clustertime[1, ])
clustertimeveccavg = apply(clustertime_veccs[1 , , ], 2, median)








load("origvsvecchia.Rdata")

replicates = 25
misclasfull = array(dim = c(2, replicates))
mvecchia = c(floor(p/20), floor(p/12), floor(p/10))
for(i in 1:nrow(l_f_matrix)){
  for(j in 1:replicates){
    
    misclasfull[i, j] = NMI(clusterdata_full[i, j, ],
                            rep(1:2, each = nsamp))
  }
}

misclasvecchia = array(dim = c(nrow(l_f_matrix), replicates, length(mvecchia)))

for(i in 1:nrow(l_f_matrix)){
  
  for(j in 1:replicates){
    
    for(k in 1:length(mvecchia)){
      
      misclasvecchia[i, j, length(mvecchia) - k + 1] = 
        NMI(clusterdata_veccs[i, j, k, ],
            rep(1:2, each = nsamp))
    }
  }
}

groups2 = as.factor(rep(1:4, each = replicates))

y2 = c(c(misclasfull[2, ]), c( misclasvecchia[2, , ]))

label2 = c("Full EM",
          
          "VEM (m = 30)",
          "VEM (m = 25)",
          "VEM (m = 15)")
pdf("boxplot2.pdf", height = 4, width = 6)
par(mar = c(1, 3, 1, 1), oma = c(1, 0, 0, 0))
boxplot(y ~ groups2, ylab = "", names = label2, 
        col = RColorBrewer::brewer.pal(length(label2), "Set2"))
mtext("y", side = 2, line = 2)
dev.off()

clustertimfullavg = median(clustertime[2, ])
clustertimeveccavg = apply(clustertime_veccs[2 , , ], 2, median)