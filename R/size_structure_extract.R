size_structure_extract = function(x) {
    out = list() 
    for (i in 1:length(x)) {
        kk = names( x[[i]] )
        out[[i]] = list(
            space = x[[i]] [["space"]],
            time =  x[[i]] [["time"]],
            peaks = ifelse( "peaks" %in% kk, x[[i]] [["peaks"]], NA ),
            troughs = ifelse( "troughs" %in% kk, x[[i]] [["troughs"]], NA ),
            peak_values = ifelse( "peak_values" %in% kk, x[[i]] [["peak_values"]], NA ),
            trough_values = ifelse( "trough_values" %in% kk, x[[i]] [["trough_values"]], NA ),
            N = x[[i]] [["N"]]
        )
    }
    return (out)
}
