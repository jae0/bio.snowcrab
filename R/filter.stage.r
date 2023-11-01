
filter.stage = function(M, mds=NULL) {

    M$logcw = log(M$cw)

    M$sex= as.character(M$sex)
    M$sex[ which(M$sex=="0") ] = "m"
    M$sex[ which(M$sex=="1") ] = "f"
 
    M$mat= as.character(M$mat)
    M$mat[ which(M$mat=="0") ] = "i"
    M$mat[ which(M$mat=="1") ] = "m"
 
    out = rep( NA, nrow(M))
    for ( j in 1:nrow(mds) ) {
        i = M[ sex==mds$sex[j] & mat==mds$mat[j] & 
            logcw >= mds$predicted_lb[j] & logcw < mds$predicted_ub[j], which=TRUE]
        if ( length(i) > 0 ) out[i] = mds$stage[j]
    }
    return(out)
}
