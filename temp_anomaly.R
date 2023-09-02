m_avg = function(x, t = 3){
  
  n = length(x)
  y = rep(0, n - t + 1)
  for( i in 1: length(y)){
    y[i] = mean(x[i: (i + t - 1)])
  }
  return(y)
}

anomalydat = read.csv("tempanomaly.csv")
library(RColorBrewer)
cut_pos = which(anomalydat[ , 1] == 202212)
cut_pos2 = (which(anomalydat[ , 1] >= 190101))[1]
data_modif = anomalydat[cut_pos2:cut_pos, ]

numofyears = nrow(data_modif) / 12 

data_matrix = matrix(nrow = numofyears, ncol = 12)

for(i in 1:12){
  
  positions = which((as.double(data_modif[ ,1]) %% 100) == i)
  data_matrix[ ,i] = data_modif[positions, 2]
}

data_mean = apply(data_matrix, 2, mean)

data_matrix_centered = apply(data_matrix_cent, 1, function(x) x - data_mean)
datamat_ma = apply(data_matrix_centered, 1, function(x) m_avg(x, t = 5))
X = 1:nrow(datamat_ma)/nrow(datamat_ma)
datamat_ma2 = datamat_ma/sd(datamat_ma)
cl_check = gmm.fromscratch.v1.Matern(X = X, Y = t(datamat_ma2), k = 3, itern_em = 50)

p = nrow(datamat_ma2); m = 20
orgordering = 1:p; locs = 1:p/p 
vecchialocspecs = getLocsVecchia(locs = matrix(locs), 
                                 ordering = "maxmin", m = m, seed = 1234)
X = X[vecchialocspecs$MatchOrigtoPerms]; Y = (t(datamat_ma2))[ , vecchialocspecs$MatchOrigtoPerms]
cl_check_vecchia = gmm..v2.Matern.vecchia(X = X, Y = Y, k = 3, vecchia_obj = vecchia_obj, itern_em = 50)

save(cl_check_vecchia, cl_check, data_matrix, file = "tempanomalyclustering.Rdata")

#plots
plot_X = 1900 + (1:((cut_pos - cut_pos2 + 1)/12))

cl_l_f = cl_check_vecchia$l_f
cl_sig_f = cl_check_vecchia$sig_f
cl_sig_n = cl_check_vecchia$sig_n
cl_labels = cl_check_vecchia$cluster
k = 3
cl_means = matrix(0, nrow = length(plot_X), ncol = k)
for(i in 1:k) cl_means[ , i] = 
  apply(data_matrix[, which(cl_labels == i)], 1, mean)

plot_df = data.frame(Year = plot_X,
                     anomaly = c(cl_means),
                     cl = as.factor(rep(1:k, each = length(plot_X))))
library(ggplot2); library(ggsci); ibrary(ggpubr)


pdf("tempanomalycluster.pdf", width = 6, height = 3)
ggplot(data = plot_df, aes(x = Year, y = anomaly, color = cl)) + 
  geom_line() + theme_bw() + 
  scale_color_aaas(name = "", labels = c("Summer", "Winter", "Transition months")) + 
  xlim(1901,2022) +
  theme(legend.position = c(0.29, 0.916),
        legend.background = element_rect(colour = "black"),
        legend.direction = "horizontal") +
  ylab("Anomaly (in celcius)") +
  expand_limits(x = c(1905, 2020))
dev.off()

pdf("tempanomalyoriginal.pdf", width = 6, height = 3)
ggplot(data = data.frame(Year = rep(plot_X, 12), 
                         y = c(data_matrix), 
                         month = as.factor(rep(1:12, each = length(plot_X)))),
       aes(x = Year, y = y, color = month)) + 
  geom_line() + theme_bw() + theme(legend.position = "None") +
  scale_color_manual(values = get_palette("lancet", 12)) +
  ylab("Anomaly (in celcius)")
dev.off()