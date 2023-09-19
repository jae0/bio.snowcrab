size_structure_extract = function(x) {
    out = list() 
    for (i in 1:length(x)) {
        kk = names( x[[i]] )
        out[[i]] = list(
            space = x[[i]] [["space"]],
            time =  x[[i]] [["time"]],
            N = x[[i]] [["N"]],
            peaks = ifelse( "peaks" %in% kk, x[[i]] [["peaks"]], NA ),
            peak_values = ifelse( "pkvalue" %in% kk, x[[i]] [["pkvalue"]], NA ),
            troughs = ifelse( "troughs" %in% kk, x[[i]] [["troughs"]], NA )
        )
    }
    return (out)
}
