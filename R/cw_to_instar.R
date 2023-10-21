
cw_to_instar = function(logcw, sex) {
    if (sex=="f") {
        # cw = exp(2.198848 + 0.315026 * (instar - 4) )
        instar = floor( 4 + ( logcw - 2.198848 ) /  0.315026 )
    }
    if (sex=="m") {
        # cw = exp(1.917564 + 0.298914 * (instar - 3) )
        instar = floor( 3 + ( logcw - 1.917564 ) /  0.298914 )
    }
    return(instar)
}
