
identify_neigbours = function( nb, index, n_neighbours=1 ) {
    out = index
    for (k in 1:n_neighbours) {
        aoi = out
        for (j in aoi) {
            out = sort( unique( c(aoi, nb$nbs[[j]]) ) )
        }
    }
    return(out)
}