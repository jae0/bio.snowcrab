
  # ------------------------
  # mean body size

  bodysize = function(x, factors=c("trip", "set"), variable, logtransform=TRUE) {

    if (logtransform) x[[variable]] = log10(x[[variable]] )

    x = x[is.finite(get(variable)),]

    m = x[, .(
        mn=mean(get(variable), na.rm=T),
        vr=var(get(variable), na.rm=T),
        n =.N
      ),
      by = factors 
    ]

    names(m) = c(factors, paste(variable, c("mean", "var", "n"), sep="."))
  
    return(m)
  }


