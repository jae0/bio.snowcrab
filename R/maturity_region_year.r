
 

maturity_region_year = function( p  ) {
    

    regions = c("cfanorth", "cfasouth", "cfa4x")
    nregions = length(regions)
  

    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    # maturity codes
    immature = 0
    mature = 1
    mat.unknown = 2

    det = snowcrab.db( p=p, DS="det.georeferenced" )
    setDT(det)
    det = det[det$sex %in% c(0,1), ]
    det = det[det$mat %in% c(0,1), ]
    det = det[det$yr >= 1999, ]
 
    det$region = NA
    for ( reg in regions) {
        r = polygon_inside(x = det, region = aegis.polygons::polygon_internal_code(reg), planar=FALSE)
        det$region[r] = reg
    }
    
    det$region = factor(det$region, levels=regions )
    det$sex = as.factor(det$sex)
 
    det$yr = as.factor(det$yr)

    output = NULL
 #   r =  glm(mat ~ cw + region + sex, data=det, family=binomial(link="logit") )
    rM =  glm(mat ~ cw + region + yr+ region*cw + yr*cw, data=det[sex=="0",], family=binomial(link="logit") )
    rF =  glm(mat ~ cw + region + yr+ region*cw + yr*cw, data=det[sex=="1",], family=binomial(link="logit") )

    outM = expand.grid(cw=40:130, region=unique(det$region), yr=unique(det$yr))
    outF = expand.grid(cw=30:80, region=unique(det$region), yr=unique(det$yr))
    
    outM$mat = predict(rM, outM, type="response")
    outF$mat = predict(rF, outF, type="response")
 

    require(ggplot2)

    color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
    
    
    outMM = ggplot(outM, aes(x=cw, y=mat, fill=yr, colour=yr ) )+
      geom_line( alpha=0.9, linewidth=1.2 ) +
      # geom_point(aes(shape=region), size=3, alpha=0.7 ) +
      # geom_errorbar(aes(ymin=cw50lower,ymax=cw50upper), linewidth=0.8, alpha=0.8, width=0.3)  +
      labs(x=NULL, y=NULL) +
      # labs(x="Year", y="", size = rel(1.5)) +
      # scale_colour_manual(values=color_map) +
      # scale_fill_manual(values=color_map) +
      # scale_shape_manual(values = c(15, 17, 19)) +
      scale_colour_hue() +
      theme_light( base_size = 22) + 
      theme( 
        legend.text = element_text(size=16),
        legend.title=element_blank(),
        legend.key.size = unit(0.1, "cm"),
        legend.position = "bottom",
        legend.justification="left"
      ) + 
      geom_hline(yintercept=0.5)  + facet_grid(region ~ .)
      fn = file.path(  p$datadir, "output", "size_at_maturity_male.png" )
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=outMM, device="png", width=12, height = 8)
      print( fn )
  
    outFF = ggplot(outF, aes(x=cw, y=mat, fill=yr, colour=yr ) )+
      geom_line( alpha=0.9, linewidth=1.2 ) +
      # geom_point(aes(shape=region), size=3, alpha=0.7 ) +
      # geom_errorbar(aes(ymin=cw50lower,ymax=cw50upper), linewidth=0.8, alpha=0.8, width=0.3)  +
      labs(x=NULL, y=NULL) +
      # labs(x="Year", y="", size = rel(1.5)) +
      # scale_colour_manual(values=color_map) +
      # scale_fill_manual(values=color_map) +
      # scale_shape_manual(values = c(15, 17, 19)) +
      scale_colour_hue() +
      theme_light( base_size = 22) + 
      theme( 
        legend.text = element_text(size=16),
        legend.title=element_blank(),
        legend.key.size = unit(0.1, "cm"),
        legend.position = "bottom",
        legend.justification="left"
      ) + 
      geom_hline(yintercept=0.5)  + facet_grid(region ~ .)
      fn = file.path(  p$datadir, "output", "size_at_maturity_female.png" )
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=outFF, device="png", width=12, height = 8)
      print( fn )
  
    return( list( predictions=out, male_plot=outMM, female_plot=outFF ) )
}
 