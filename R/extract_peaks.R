extract_peaks = function(x, nm="peaks") {
    out = NULL
    for (i in 1:length(x)) {
        kk = names( x[[i]] )
        if (length(kk) > 0 ) {
            if ( nm %in% kk ) {
                out = c(out,  x[[i]] [[nm]] )
            }
        } 
    }
    return (out)
}
