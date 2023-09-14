
identify_neigbours = function( nb, index, num_levels=1 ) {
    out = index
    for (k in 1:num_levels) {
        aoi = out
        for (j in aoi) {
            out = sort( unique( c(aoi, nb$nbs[[j]]) ) )
        }
    }
    return(out)
}