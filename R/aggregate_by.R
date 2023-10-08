
aggregate_by = function(M, 
    agg_by = c("year", "sex", "mat", "region", "zi", "ti" ), 
    xvals = NULL, 
    agg_function = function(x)  { mean( x, na.rm=TRUE )  },
    recale_density_to_numerical_density=TRUE, 
    add_offset=NULL ) {

    # names of variables with KD estimates "V###"
    cols = colnames(M)
    cols = cols[grep("^V[[:digit:]]+",cols)]
 
    
    agg_formula = paste( "variable ~", paste( agg_by,  collapse="+"  ) )
    
    if (recale_density_to_numerical_density) {
        # simple transformations (-> no/km^2)
        M = cbind(M[, ..agg_by], M[ , .SD * Neffective , .SDcols =cols ] )
    }

    # aggregation
    if (!is.null(add_offset)) {
        if (is.logical(add_offset)) {
            if (add_offset) {
                # mimic observation error (lowest nonzero value)
                add_offset = zapsmall( as.matrix(M[, ..cols ] ) )
                add_offset = min( add_offset[add_offset>0], na.rm=TRUE)
            } 
        }
        M = cbind(M[, ..agg_by], M[,  zapsmall(.SD + add_offset), .SDcols =cols  ] )
    }
 

    MA = M[ ,  lapply(.SD, agg_function ), .SDcols =cols, by=agg_by  ]

    # reshape    
    kdb = dcast( melt(MA, id.vars=agg_by ), agg_formula )
    
    if (is.null(xvals)) xvals=1:nrow(kdb)
    kdb$logcw = xvals
    kdb$cw = exp(kdb$logcw)
    
    return(kdb)
}
