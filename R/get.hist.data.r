    get.hist.data = function(set, det, save=F) {
      dd = det[is.finite(det$cw),]
      hdata = merge(dd, set, by=c("trip", "set"), all.x=T, all.y=F, sort=F)
      if (save) read_write_fast(data=hdata, fn="hist.rdz" )
      return (hdata)
    }


