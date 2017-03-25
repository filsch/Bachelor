source('~/Documents/Bachelor/auxiliary_function_procedure.R')
source('~/Documents/Bachelor/auxiliary_function_stats.R')

tempdata = read.table('~/Documents/Bachelor/tempdata.txt', header = TRUE)
sail_lines = read.table('~/Documents/Bachelor/sail_lines.txt')

names(tempdata)[1] <- 'x'; names(tempdata)[2] <- 'y'; names(tempdata)[3] <- 'z'
tempdata = data.frame(x = c(1, tempdata$x), y = c(1, tempdata$y), z = c(5.513, tempdata$z))

map = reshapeMap(map = tempdata, grid_size=100, type='reduced')

new_tempdata = xyz.to.grid(tempdata)
par(mfrow=c(1,2))
image.plot(new_tempdata)
image.plot(map)


