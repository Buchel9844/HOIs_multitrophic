#################################################-
## Code from William Petry ----
# https://github.com/wpetry/StructuralCoexistence/blob/master/server.R
# copied the 19 sept 2020
#################################################-

# mk_graph_3sp ====
# make interaction network plot

mk_graph_3sp <- function(alphamat, rs){
    rs = rs/max(rs)
    alphamat.pos <- unname(abs(alphamat))
    g <- igraph::graph_from_adjacency_matrix(alphamat.pos  > 0)
    E(g)$weight <- as.numeric(alphamat.pos )
    widths <- E(g)$weight *4
    color.mat <-  -as.numeric(alphamat)
    color.mat[which(color.mat > 0)] <- 1
    color.mat[which(color.mat < 0)] <- 3
    E(g)$lty <- color.mat
    #widths[widths > 1] <- sqrt(widths)
    as.ggplot(function() plot(g,
         #main = "Species interaction network",
         #main = title,
         margin = c(0, -0.15, 0, -0.15),
         xlim = c(-1.25, 1.25), 
         ylim = c(-1.25, 1.25),
         vertex.label = c("Radish","Field   \nbean   ", "Tomato"),
         vertex.label.family="Helvetica",   
         vertex.label.cex = 1,
         vertex.label.dist= 5,
         vertex.label.degree = c(pi/2,pi/2,-pi/2),
         vertex.label.color = "black",
         vertex.size = 10 * rs,
         vertex.color = "grey80",
         vertex.frame.color = "transparent",
         edge.curved = TRUE,
         edge.width = widths,
         edge.arrow.size = 2,
         edge.arrow.mode = c(0, 2, 2,
                             2, 0, 2,
                             2, 2, 0),
         edge.color = "black",
         edge.loop.angle = 0.75,
         layout = matrix(c(4, 0, 0, 0, 2, sqrt(3)/2), ncol = 2,
                         byrow = TRUE)))
   # int.network <- list()
    #int.network[["1"]] <- plot.int.network
    #return(int.network[["1"]])
}

