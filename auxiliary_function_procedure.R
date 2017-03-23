#Auxilliary functions for bachelor-thesis, NOT relevant for stats

#----------------------------------------------------------------------------------------

#Constructing square map with adapted axes as a matrix
constructSquareMap <- function(map, grid_size){
  square_map = interp(x = map$x, y = map$y, z = map$z,
                      xo = seq(1,max(map$x),length.out = grid_size + 2)[2:(grid_size + 1)],
                      yo = seq(1,max(map$y),length.out = grid_size + 2)[2:(grid_size + 1)])$z
  #square_map = apply(square_map,1,rev)
  return(square_map)
}

#----------------------------------------------------------------------------------------
#Opposite of grid.to.xyz, means loss of "true" coordinates
xyz.to.grid <- function(map){
  nrows = max(map$x)
  ncolumns = max(map$y)
  grid = matrix(numeric(nrows*ncolumns), nrow = nrows, ncol = ncolumns)
  sorted_mapx = sort(unique(map$x))
  sorted_mapy = sort(unique(map$y))
  for (i in 1:nrows){
    for (j in 1:ncolumns){
      grid[i,j] = map$z[intersect(which(map$x == sorted_mapx[i]), which(map$y == sorted_mapy[j]))]
    }
  }
  return(grid)
}
#----------------------------------------------------------------------------------------

#Constructing distance matrix in dim_samplepoints x dim_prediction_points
constructDistanceMatrix <- function(sample_points, prediction_points){
  
  nrows = dim(prediction_points)[1]; ncolumns = dim(sample_points)[1]
  dist_matrix = matrix(numeric(nrows*ncolumns),nrow=nrows, ncol=ncolumns)
  zero_matrix = matrix(numeric(2*ncolumns),nrow=2)
  
  #Calculating distances, one column of the matrix at the time
  for (i in 1:nrows){
    dist_matrix[i,] = sqrt( colSums((t(sample_points) - (zero_matrix + as.numeric(prediction_points[i,])) )^2))
  }
  return(dist_matrix)
}

#----------------------------------------------------------------------------------------

#Samples n samples with position from the grid
#(To check the values: z[ intersect(which(square_map$y==a),which(square_map$x==b))) ]
gridSampler <- function(n, map, design, noise = 0){
  xmax = dim(map)[1]; xmin = 1; ymax = dim(map)[2]; ymin = 1;
  #xmax = max(map$x);   ymax = max(map$y);   xmin = min(map$x);   ymin = min(map$y)
  if (xmin == 0){
    xmin = 1
  }
  if (ymin == 0){
    ymin = 1
  }
  
  if (design == 'regular'){
    #Designed to have a perfect regular grid, disregarding the actual n
    xincrement = round( (xmax - xmin)/(round(sqrt(n)) - 1) )
    yincrement = round( (ymax - ymin)/(round(sqrt(n)) - 1) )
    
    samples = expand.grid(x = seq(xmin, xmax, xincrement), y = seq(ymin, ymax, yincrement))
    points = cbind(x = c(samples$x), y = c(samples$y))
    samples$z = map[points]
  } 
  else if (design == 'random'){
    x = sample(xmax, size=n, replace=TRUE)
    y = sample(ymax, size=n, replace=TRUE)
    points = cbind(x = c(samples$x), y = c(samples$y))
    samples = list(x = x, y = y, z = map[points])
  }
  
  samples$z = samples$z + rnorm(length(samples$z),0, noise)    
  return(samples)
}

#----------------------------------------------------------------------------------------
#Reshapes the points if needed for plotting
reshapePoints <- function(list_points, size_to_x, size_to_y){
}