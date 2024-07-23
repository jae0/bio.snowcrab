# -----------------------
# obtain a point estimate of maturity from the data provided
  est.size.at.maturity = function(d,sex) {

    stop ("LD50 is deprecated? use direct computation from coef: -slope/int ... check")
    # recode maturity: 2 is immature -> 0
    imm = which(d$mat==2)
    d$mat[imm] = 0

    i = which(is.finite(d$mat))

    outputvector = NULL

    if (length(i) > 20) {
      d = d[i,]

    #  d$surfacearea = d$surfacearea / min(d$surfacearea[d$surfacearea>0] )
      dummy = d[1,]
    #  dummy$surfacearea = floor(length(i)/10)
      dmat = dimm = dummy
      dmat$mat = 1
      dimm$mat = 0
      if (sex==male) {
        dimm$cw = 30
        dmat$cw = 150
      } else {
        dimm$cw = 30
        dmat$cw = 80
      }

      d = rbind(d, dmat, dmat, dmat, dimm, dimm, dimm, dimm)

      r =  glm (mat ~ cw, data=d, family=binomial(link="logit") )
      res = coef(summary(r))
      res2 = confint(r)
      fin = cbind(res[,c("Estimate","Std. Error")],res2)
      colnames(fin) = c("Estimate","Std. Error", colnames(res2))
      cw50 = doBy::dose.LD50(r,lambda=c(1,NA)) # lambda is a vector of model coef, with NA for inverse prediction variable.
      names(cw50) = c("cw50", "cw50lower", "cw50upper")
      olist = c("Estimate","Std. Error")

      r1 = res[1,olist] ; names(r1) = c("a0", "a0.se" )
      r2 = res[2,olist] ; names(r2) = c("a1", "a1.se" )
      r3 = c(length(i), r$deviance, r$aic); names(r3)=c("n", "deviance", "aic")
      r4 = c(floor(r$converged)); names(r4)="converged"

      outputvector = c( r1, r2, r3, r4, cw50)

    }
    return(outputvector)
  }



