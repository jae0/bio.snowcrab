
  errorbars = function(b, v, br="Sturges", nfilter=0, lowess=0.4, trim=0, xlab, ylab, p, ...)  {
    # ... params to send to errbar (labels, lty, etc)
    lookup.table = snowcrab.db( p=p, DS="data.transforms" )

    b[,v[2]] = bio.snowcrab::variable.recode( b[,v[2]], v[2], direction="forward",  lookup.table=lookup.table )

    good = which(is.finite(rowSums(b[,v], na.rm=T)))
    x = b[good, v[1]]
    y = b[good, v[2]]
    h = hist(x, breaks=br, plot=F)
    mids = h$mids
    delta = h$breaks[2] - h$breaks[1]
    d = cut(x, breaks=h$breaks)
    e =  tapply(X=y, INDEX=d, FUN=mean, simplify=T, na.rm=T, trim=trim)
    e =  data.frame(mean=e)
    f =  tapply(X=y, INDEX=d, FUN=var, simplify=T, na.rm=T)
    f =  data.frame(var=f)
    g =  tapply(X=y, INDEX=d, FUN=length, simplify=T)
    g =  data.frame(n=g)
    out = cbind(e,f,g)
    out$se = sqrt(out$var/(out$n-1))
    out$mids = h$mids
    out = out[ out$n > nfilter ,]
    r = which (is.finite(out$mids+out$mean) )

    yplus = bio.snowcrab::variable.recode( x=(out$mean+out$se), v[2],  direction="backward", lookup.table=lookup.table)
    yminus = bio.snowcrab::variable.recode( x=(out$mean-out$se), v[2],  direction="backward", lookup.table=lookup.table )
    out$mean = bio.snowcrab::variable.recode( x=(out$mean), v[2],  direction="backward", lookup.table=lookup.table )

    errbar( x=out$mids, y=out$mean, yplus=yplus, yminus=yminus,
            xlab=xlab, ylab=ylab, xlim=c(min(out$mids, na.rm=T)-delta, max(out$mids, na.rm=T)+delta),  axes=F )  # from Hmisc

   lines(lowess(out$mids[r], out$mean[r], f=lowess), lty="solid", col="orange", lwd=4 )
#    lines(
#        loess(out$mean[r] ~ out$mids[r], weights=out$n[r], na.action="na.omit", span=lowess),
#        lty="solid", col="orange", lwd=4
#        )
      axis( 1 )
      axis( 2 )
      abline(h=0, lty=3)

    text(  x=out$mids, y=out$mean, labels=out$n, pos=4 )
    points(  x=out$mids, y=out$mean )

    return (out)
  }


