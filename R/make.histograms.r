

  make.histograms = function(set, det, hvar="cw", breaks=NULL) {

    setDT(det)
    setDT(set)

    det$breaks = cut(det$cw, breaks=breaks, labels=breaks[-1] )
    o = det[ , .(counts=.N), by=.(sid, breaks) ]
    o$breaks = factor2number( o$breaks)
    out0 = CJ( sid=sort(unique(set$sid)), breaks = breaks[-1] )
    o = o[ out0, on=.(sid, breaks) ][is.na(counts), counts := 0][]

    out = dcast( o, sid~breaks, value.var="counts", fun.aggregate=function(x) sum(x, na.rm=T) )
    out = set[,.(sid, sa)][ out, on=.(sid)] 
    
    out = out[, .SD / sa, .SDcols=3:ncol(out)]  
    setDF(out)

    rownames(out) = set$sid
    set$sid = NULL
    set$sa = NULL

    colnames(out) = breaks[-1]
    return(out)
  }


