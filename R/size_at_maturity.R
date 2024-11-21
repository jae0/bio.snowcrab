
size_at_maturity = function(p) {
   
    M = size_distributions(p=p, toget="base_data" )
    setDT(M)
    M$logcw = log(M$cw)
    
    M = M[ mat %in% c("0", "1") & sex %in% c("0", "1") & is.finite(logcw) , ]
    M$mat= as.numeric( as.character(M$mat ))

    Mf = M[sex=="1",]
    Mm = M[sex=="0",]
    
    of = glm( mat ~ logcw + t + log(z) + factor(year) + factor(region), data=Mf, family=binomial(link="logit"), na.action=na.omit  )
    AIC(of)

    om = glm( mat ~ logcw, data=Mm, family=binomial(link="logit")  )
  AIC(om)

    Mf$predicted = NA
    Mf$predicted[ -of$na.action ] = predict(of, data=Mf, type="link", na.action=na.pass)
    
    Mm$predicted = NA
    Mm$predicted[ -of$na.action ] = predict(of, data=Mm, type="link", na.action=na.pass)


    summary(of)
    cor(Mf$predicted,  Mf$mat, use="pairwise.complete.obs")  # 0.687;  0.8186
    plot(Mf$predicted,  Mf$mat)  

    summary(om)
    cor(Mm$predicted,  Mm$mat)  # 0.55
    plot(Mm$predicted,  Mm$mat)  

}
