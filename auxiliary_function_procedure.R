#Auxilliary functions for bachelor-thesis, NOT relevant for stats


#Constructing map with adapted axes as a matrix
reshapeMap <- function(map, grid_size = 0, type){
  if (type == 'square'){
    grid_x = grid_size
    grid_y = grid_size
    
  } else if (type == 'reduced'){
    dimx = max(map$x); dimy = max(map$y)
    proportion = dimx/(dimy + dimx)
    grid_x = round(grid_size*proportion)
    grid_y = grid_size - grid_x
  }
  map = interp(x = map$x, y = map$y, z = map$z,
               xo = seq(1,max(map$x),length.out = grid_x + 2)[2:(grid_x + 1)],
               yo = seq(1,max(map$y),length.out = grid_y + 2)[2:(grid_y + 1)])$z
  #square_map = apply(square_map,1,rev)
  return(map)
}

#Opposite of grid.to.xyz, means loss of "true" coordinates
xyz.to.grid <- function(map){
  nrows = max(map$x)
  ncolumns = max(map$y)
  grid = matrix(numeric(nrows*ncolumns), ncol = ncolumns, nrow = nrows)
  sorted_mapx = sort(unique(map$x))
  sorted_mapy = sort(unique(map$y))
  for (i in 1:nrows){
    for (j in 1:ncolumns){
      #cat('(',i,',',j,')\n')
      grid[i,j] = map$z[intersect(which(map$x == sorted_mapx[i]), which(map$y == sorted_mapy[j]))]
    }
  }
  return(grid)
}

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

#Samples n samples with position from the grid
#(To check the values: z[ intersect(which(square_map$y==a),which(square_map$x==b))) ]
gridSampler <- function(nx = 0, ny = 0, map, design, noise = 0){
  xmax = dim(map)[1]; xmin = 1; ymax = dim(map)[2]; ymin = 1;
  #xmax = max(map$x);   ymax = max(map$y);   xmin = min(map$x);   ymin = min(map$y)
  if (xmin == 0){
    xmin = 1
  }
  if (ymin == 0){
    ymin = 1
  }
  
  if (design == 'regular'){

    samples = expand.grid(x = seq(xmin, xmax - xmax/nx, length.out=nx),
                          y = seq(ymin, ymax - ymax/ny, length.out=ny))
    shift = round((min(samples$x) - xmin + xmax - max(samples$x))/2)
    samples$x = round(samples$x + shift)
    samples$y = round(samples$y + shift)
    points = cbind(samples$x, samples$y)
    samples$z = map[points]
  } 
  else if (design == 'random'){
    x = sample(xmax, size=n, replace=TRUE)
    y = sample(ymax, size=n, replace=TRUE)
    points = cbind(x = c(samples$x), y = c(samples$y))
    samples = list(x = x, y = y, z = map[points])
  }
  
  samples$z = samples$z + rnorm(length(samples$z),0, noise)
  samples$noise = noise
  return(samples)
}

#Adapts the sailing lines to the map you'ld like to predict on.
#Input: 
# lines             - sailing lines in the form of xy-coordinates
# prediction_map    - the map that one is to predict on, in the form of matrix grid
# original_map      - the original mapping that the saillines are from. 
#Output:
# adapted_lines     - sailing lines in the form of xy-coordinates adapted to the axes of prediction_map
adaptLines <- function(lines, prediction_map, original_map){
  dims_p = dim(prediction_map)
  dims_o = dim(original_map)
  
  ratiox = dims_p[1]/dims_o[1]
  ratioy = dims_p[2]/dims_o[2]
  
  lines$x = round(lines$x*ratiox)
  lines$y = round(lines$y*ratioy)
  return(unique(lines))
}