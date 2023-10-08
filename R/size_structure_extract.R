size_structure_extract = function(KD,
    sexes=c("0", "1"), 
    mats=c("0", "1"), 
    regions=c("cfanorth", "cfasouth", "cfa4x"),
    aus=NULL,
    years=as.character(1996:2023),
    zlevels=c(0, 100),
    tlevels= c(-2, 6),
    strata="yasm",
    ...
    ) {

    peaks = data.table()
    troughs = data.table()
    peak_values = data.table()
    trough_values = data.table()

    if (strata=="smryzt") {

        for ( s in sexes ) {
        for ( m in mats ) {
        for ( r in regions) {
        for ( y in years ) {
        for ( z in zlevels ) {
        for ( t in tlevels ) {
            vn = paste(y,s,m,r,z,t, sep="_" )
            if (!exists(vn, KD)) next()

            mds = identify_modes( Z=as.vector(t(KD[, ..vn])), ... )   
            if (is.null(mds)) next()
            if (inherits(mds, "try-error")) next()
            # out[[s]][[m]][[r]][[y]][[z]][[t]] = mds
            peaks = rbind(peaks, cbind(s, m, r, y, z, t, t(t(mds[["peaks"]])) ))
            troughs = rbind(peaks, cbind(s, m, r, y, z, t, t(t(mds[["troughs"]])) ))
            peak_values = rbind(peak_values, cbind(s, m, r, y, z, t, t(t(mds[["peak_values"]])) ))
            trough_values = rbind(trough_values, cbind(s, m, r, y, z, t, t(t(mds[["trough_values"]])) ))
        }}} }}}
        
        setnames(peaks, "V7", "peaks")
        setnames(troughs, "V7", "troughs")
        setnames(peak_values, "V7", "peak_values")
        setnames(trough_values, "V7", "trough_values")

    } else if (strata=="yasm") {    
  
        for ( s in sexes ) {
        for ( m in mats ) {
        for ( a in aus) {
        for ( y in years ) {
            vn = paste(y,a,s,m, sep="_" )
            if (!exists(vn, KD)) next()
            mds = identify_modes( Z=as.vector(t(KD[, ..vn])), ... )   
            if (is.null(mds)) next()
            if (inherits(mds, "try-error")) next()
            peaks = rbind(peaks, cbind(s, m, a, y, t(t(mds[["peaks"]])) ))
            troughs = rbind(peaks, cbind(s, m, a, y, t(t(mds[["troughs"]])) ))
            peak_values = rbind(peak_values, cbind(s, m, a, y, t(t(mds[["peak_values"]])) ))
            trough_values = rbind(trough_values, cbind(s, m, a, y, t(t(mds[["trough_values"]])) ))
        }}} }
 
        setnames(peaks, "V5", "peaks")
        setnames(troughs, "V5", "troughs")
        setnames(peak_values, "V5", "peak_values")
        setnames(trough_values, "V5", "trough_values")
         
    }

    peaks$peaks = as.numeric(peaks$peaks)
    peak_values$peak_values = as.numeric(peak_values$peak_values)
    troughs$troughs = as.numeric(troughs$troughs)
    trough_values$trough_values = as.numeric(trough_values$trough_values)

    out = list()
    out$peaks=peaks
    out$peak_values=peak_values
    out$troughs=troughs
    out$trough_values=trough_values

    return(out)
 
}
