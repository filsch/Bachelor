source('~/Documents/Bachelor/auxiliary_function_procedure.R')
source('~/Documents/Bachelor/auxiliary_function_stats.R')

n = 60;
tempdata = read.table('~/Documents/Bachelor/tempdata.txt', header = TRUE)
sail_lines = read.table('~/Documents/Bachelor/sail_lines.txt')

names(tempdata)[1] <- 'x'; names(tempdata)[2] <- 'y'; names(tempdata)[3] <- 'z'
sail_lines = list(run1 = data.frame(y = tempdata$x[sail_lines$V1], x = tempdata$y[sail_lines$V1]),
                  run2 = data.frame(y = tempdata$x[sail_lines$V2], x = tempdata$y[sail_lines$V2]),
                  run3 = data.frame(y = tempdata$x[sail_lines$V3], x = tempdata$y[sail_lines$V3]),
                  run4 = data.frame(y = tempdata$x[sail_lines$V4], x = tempdata$y[sail_lines$V4]),
                  run5 = data.frame(y = tempdata$x[sail_lines$V5], x = tempdata$y[sail_lines$V5]),
                  run6 = data.frame(y = tempdata$x[sail_lines$V6], x = tempdata$y[sail_lines$V6]),
                  run7 = data.frame(y = tempdata$x[sail_lines$V7], x = tempdata$y[sail_lines$V7]),
                  run8 = data.frame(y = tempdata$x[sail_lines$V8], x = tempdata$y[sail_lines$V8]))
tempdata = data.frame(x = c(1, tempdata$x), y = c(1, tempdata$y), z = c(5.513, tempdata$z))

par(mfrow=c(2,2))
for (i in 1:8){
i = attr(sail_lines,"names")[i]
runs = matrix(numeric(max(tempdata$x)*max(tempdata$y)),ncol=max(tempdata$x),nrow=max(tempdata$y))
runs[cbind(sail_lines[[i]]$x,sail_lines[[i]]$y)] = runs[cbind(sail_lines[[i]]$x,sail_lines[[i]]$y)] + 1
image.plot(runs, col=c('black','lightgreen'))
}
map = reshapeMap(map = tempdata, grid_size = n, type='reduced')

new_tempdata = xyz.to.grid(tempdata)
par(mfrow=c(1,2))
image.plot(new_tempdata)
image.plot(map)



