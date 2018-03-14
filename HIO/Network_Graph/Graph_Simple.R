# Libraries
rm(list=ls())
library(igraph)


# Define system 
dim_states <- 6

names_states <- c()
for (i in 1:dim_states) {
  names_states <- c(names_states,paste('x',i, sep = ''))
}
  

A = matrix( data = c(      0 ,       0 ,       0 ,      0 ,      0 ,      0 ,
                           0 ,       0 ,       0 ,      0 ,      0 ,      0 ,
                           0 ,     'f' ,       0 ,      0 ,      0 ,      0 ,
                           0 ,       0 ,     'f' ,      0 ,      0 ,      0 ,
                         'f' ,       0 ,     'f' ,      0 ,      0 ,      0 ,
                           0 ,       0 ,       0 ,    'f' ,      0 ,      0 )
            , byrow = TRUE , nrow = dim_states )



# Define Edges and Knots

directed <- TRUE

edges <- c()
for (i in 1:dim_states) {
  for (j in 1:dim_states) {
    if (A[i,j] != '0') {
      edges <- c(edges, c( names_states[j]  , names_states[i] )  )
    }
  }
}


# Test

G1 <- graph( edges =  edges, directed = directed)

plot(G1, vertex.label.cex=1, vertex.label.dist=-4, edge.curved=0.13)