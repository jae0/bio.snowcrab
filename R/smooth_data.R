
smooth_data = function(x, y, sk = kernel("modified.daniell", c(3, 1)) ){
    nz = length(y)
    z = data.frame(x = 1:nz, y = y)
    nd = length(sk$coef) - 1
    vv = unique(diff(z$x))
    xr = min(vv[vv > 0])
    uu = unique(round(vv, abs(log10(xr))))
    if (length(uu) > 1) {
        ifunc = approxfun(x = z$x, y = z$y, method = "linear", 
            rule = 2)
        r0 = range(z$x, na.rm = TRUE)
        fx = seq(r0[1], r0[2], by = min(uu))
        fy = ifunc(fx)
    } else if (length(uu) == 1) {
        fx = z$x
        fy = z$y
    }
    fn = length(fx)
    si = (nd + 1):(fn - nd)
    fp = rep(NA, fn)
    fp[si] = kernapply(as.vector(z$y), sk)
    pfunc = approxfun(x = fx, y = fp, method = "linear", rule = 2)
    return( pfunc(z$x) )
}

