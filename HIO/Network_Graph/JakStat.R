# Libraries
rm(list=ls())
library(igraph)


# Define system 
dim_states <- 4

names_states <- c()
for (i in 1:dim_states) {
  names_states <- c(names_states,paste('x',i, sep = ''))
}
  

A = matrix( data = c(    'a' ,       0 ,       0 ,      0 ,  
                         'b' ,     'a' ,       0 ,      0 ,   
                           0 ,     'b' ,     'a' ,      0 ,   
                           0 ,       0 ,     'b' ,      0 )
            , byrow = TRUE , nrow = dim_states )



# Define nominal Edges and Knots

directed <- TRUE

edges <- c()
for (i in 1:dim_states) {
  for (j in 1:dim_states) {
    if (A[i,j] != '0') {
      edges <- c(edges, c( names_states[j]  , names_states[i] )  )
    }
  }
}

# make the IGraph Object
G <- graph( edges =  edges, directed = directed)
E(G)$weight <- 1 # nominal model
V(G)$weight <- 1

# Define hidden edge
edge_w <- c('x4','x1')
# add to the IGraph Object
G <- add_edges(G, edge_w, weight = 0)
#E(G)[length(E(G))]$weight <- 0 # not in the nominal model

# define known input vertex and edge
G <- add_vertices(G, nv = 1,attr = c(name = 'u') , weight = 2)
G <- add_edges(G,c('u','x1'), weight = 2)
#E(G)[length(E(G))]$weight <- 2 # known input
#V(E)['u']$weight <- 2 # known input

# Plot
# Set nominal edges blue and hidden edges red
E(G)$color[E(G)$weight == 0] <- 'red'    # hidden
E(G)$color[E(G)$weight == 1] <- 'blue'   # nominal
E(G)$color[E(G)$weight == 2] <- 'orange' # known input

V(G)$color[V(G)$weight == 1] <- 'black'   # nominal
V(G)$color[V(G)$weight == 2] <- 'white' # known input

elength <- c(1,1,1,1,1,1,5,2)



plot(G, weight.node.edge.dist = 19  ,vertex.label.cex=1, vertex.label.dist=4, edge.curved=0.19)