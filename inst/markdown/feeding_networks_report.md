---
title: "Snow crab predators"

keywords: 
  - basic feeding networks
abstract: |
  Details of stomach contents data base. 

metadata-files:
  - _metadata.yml

params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
  todo: [fishery_results,ecosystem]

--- 


<!-- to render a presentation using revealjs 

make quarto FN=feeding_networks_report.md DOCTYPE=revealjs  PARAMS="-P year_assessment:2024 -P todo:[fishery_results,ecosystem,redo_data]"  --directory=~/bio/bio.snowcrab/inst/markdown

-->
 
 

<!-- NOTES: 
Data are in ./data/FH.*.csv
 
Here is the requested data. There is some variability in the resolution of records, hence you’ll find that some predators may have full weight information while others do not, and the same will apply to prey species.  This will limit how informative this data will be and how far you can use these analytically.

The new database schema includes unique identifiers to ensure records are as intended (prior to this there was a lot of duplication of records).  These are as described:
SET_SEQ – relates to the set information details
PRED_SEQ – relates to the predator fish the stomach was taken from
PREY_SEQ – relates to the identified prey within a stomach

Some column headers may need further description:
DEPTH – average set depth in meters
GEAR – set fishing gear code (list attached)
SPEC – predator id code (species list attached)
FWT – predator weight in grams
FLEN – predator length in cm
STOWGT/EMPTYWGT – stomach weights before and after processing in grams
FULLNESS – level of fullness associated with stomachs, where codes are described from 0 through 6:
0             empty - no food contents
1             less than ¼ full
2             ¼ to ½ full
3             ½ to ¾ full
4             ¾ full to full
5             everted 
6             regurgitated
FGEN – Gender code of predator:
0             Unknown
1             Male
2             Female
PREYSPECCD – prey id code (species list attached)
PWT – prey weight in grams
PLEN – prey length in cm
DIGESTION – level of digestion of prey, where:
1             Good Condition
2             Partly Digested
3             Well Digested
4             Unidentifiable
9             Unidentified digestion state
 
-->


## First set up environment and get the data. 


```{r}
#| eval: true
#| output: false
#| warning: false
#| error: false 
# startup libraries and directories

require(aegis)

year_assessment = params$year_assessment

p = bio.snowcrab::load.environment( year.assessment=year_assessment )

project_directory = p$data_root

data_dir = file.path( project_directory, "data", "diets" )

require(data.table) # for speed
require(lubridate)
require(stringr) 
require(gt)  # table formatting
require(ggplot2)

require(aegis) # map-related 
require(bio.taxonomy)  # handle species codes

# assimilate the CSV data tables:
# diet = get_feeding_data( data_dir, redo=TRUE )  # if there is a data update
diet = get_feeding_data( data_dir, redo=FALSE )

```

<!--
  # misc options for aegis should we need them: 

  bycatch_dir = file.path( p$annual.results, "bycatch")
  years = as.character(1996: year_assessment)

  if (0) {
    loadfunctions( "aegis")
    loadfunctions( "bio.snowcrab")  # in case of local edits
  }

-->



## What are the predators of snow crab?

```{r}
#| eval: true
#| output: true
#| label: predator_list
#| tbl-cap: "All predators of snow crab on Scotian Shelf of Atlantic Canada. Of 58287 finfish stomach samples, 159 had snow crab (0.28%). There is no indormation on snow crab diet in the database."
tx = taxa_to_code("snow crab")  
# matching codes are 
#  spec    tsn                  tx                   vern tx_index
#1  528 172379        BENTHODESMUS           BENTHODESMUS     1659
#2 2522  98427        CHIONOECETES SPIDER QUEEN SNOW UNID      728
#3 2526  98428 CHIONOECETES OPILIO        SNOW CRAB QUEEN      729
# 2 and 3 are correct

snowcrab_predators = diet[ preyspeccd %in% c(2522, 2526), ]  # n=159 oservations out of a total of 58287 observations in db (=0.28% of all data)
snowcrab_predators$Species = code_to_taxa(snowcrab_predators$spec)$vern
snowcrab_predators$Predator = factor(snowcrab_predators$Species)
 
counts = snowcrab_predators[ , .(Frequency=.N), by=.(Species)]
setorderv(counts, "Frequency", order=-1)

gt::gt(counts) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))

```

## Where are the predators with snow crab in their stomachs?


```{r}
#| eval: false
#| output: false
#| label: temp_depth
#| fig-cap: "Temperature and depths of snow crab predation on the Scotian Shelf of Atlantic Canada. Grey is all species observations in diet data base. Red is snow crab as prey."
#| fig-dpi: 144
#| fig-height: 4
#plot( depth ~ bottom_temperature, diet, pch=18, col="grey")
#points( depth ~ bottom_temperature, snowcrab_predators, col="red", pch=19 )
```

```{r}
#| eval: true
#| output: true
#| label: map_observations_pre2010
#| fig-cap: "Locations where snow crab predation was observed on the Scotian Shelf of Atlantic Canada, 2000-2010."
#| fig-dpi: 144
#| fig-height: 6
ggplot() + 
    borders("world", fill = "lightgray", colour = "grey80") + 
    xlim(c( -65.5, -57.1)) + ylim(c(42.1, 47.1)) +
    geom_point(data=snowcrab_predators[year(timestamp) %in% c(2000:2010),], aes(x=slongdd, y=slatdd, colour=Predator ), size=2.5 ) +
    labs(x="Longitude", y="Latitude", caption="2000-2010") +
    theme(legend.position="inside", legend.position.inside =c(0.16, 0.65), legend.title=element_blank(), legend.text=element_text(size=6.0) 
 ) 
```


```{r}
#| eval: true
#| output: true
#| label: map_observations_post2010
#| fig-cap: "Locations where snow crab predation was observed on the Scotian Shelf of Atlantic Canada, 2011-2020."
#| fig-dpi: 144
#| fig-height: 6
ggplot() + 
    borders("world", fill = "lightgray", colour = "grey80") + 
    xlim(c( -65.5, -57.1)) + ylim(c(42.1, 47.1)) +
    geom_point(data=snowcrab_predators[year(timestamp) %in% c(2011:2020),], aes(x=slongdd, y=slatdd, colour=Predator ), size=2.5 ) +
    labs(x="Longitude", y="Latitude", caption="2011-2020") +
    theme(legend.position="inside", legend.position.inside =c(0.16, 0.725), legend.title=element_blank(), legend.text=element_text(size=6.0) )
    
```


## Project location of predators upon species composition information

Basic idea is that we have a reasonable undertanding of demersal communities. Based upon this, we can project where these predation events are most likely. 

Note: 

```
#| eval: false
#| output: true
# Loadings are eigenvectors scaled by the square roots of the respective eigenvalues
loadings = eigenvectors %*% sqrt(diag( eigenvalues, nrow = length( eigenvalues ))) 

# scores are sum of dot products
scores = indat %*% t(pracma::pinv(loadings ))  
```


```{r} 
#| eval: true
#| output: true
#| warning: false
#| error: false 

yrs = 1999:year_assessment
carstm_model_label="default"
require(aegis)
require(aegis.speciescomposition)

p = speciescomposition_parameters( yrs=yrs, carstm_model_label=carstm_model_label )

# speciescomposition_db( DS="speciescomposition.ordination.redo", p=p )  # analsysis
# speciescomposition_db( DS="speciescomposition.redo", p=p  ) # compute planar coords and remove dups
 
pca = speciescomposition_db( DS="pca", p=p )  # analsysis
names(pca)
# [1] "cm"               "indat"            "rotate"           "eigenvalues"     
# [5] "eigenvectors"     "loadings"         "scores"           "variance"        
# [9] "variance_percent" "total_variance"  

predator_names = unique(snowcrab_predators$Species)
 
pcares = as.data.frame( pca$loadings )
pcares$vern = taxonomy.recode( from="spec", to="taxa", tolookup=rownames( pcares ) )$vern
pcares$sc_predator = 0
pcares$sc_predator[ which(pcares$vern %in% predator_names) ] = 1   

toblank = which( pcares$sc_predator==0)
toblank = setdiff( toblank,  which(pcares$vern=="SNOW CRAB QUEEN" ) )
pcares$label = pcares$vern
pcares$label[toblank] = "+" 


plot( PC2 ~ PC1, pcares, type="n")
text( PC2 ~ PC1, labels=label, data=pcares, cex=0.8 )

pcares$uv = pcares$sc_predator
pcares$uv[ which( pcares$vern=="SNOW CRAB QUEEN")] = 1
setDT(pcares)


# re-weight species by frequency
setDT(counts)
counts[, weight := Frequency / sum(Frequency) ]

pcares$Species = pcares$vern
pcares = counts[pcares, on="Species"]
pcares$weight[ which(is.na(pcares$weight)) ] = 0
pcares$weight = pcares$weight * 3
pcares$weight[ pcares$Species =="SNOW CRAB QUEEN" ] = 3


# uvs =pcares$uv %*% t(pracma::pinv(pca$loadings ))    
# points( uvs[2] ~ uvs[1], col="red" )
# segments(0, 0, x1 = uvs[1], y1 = uvs[2], col = "red", lwd = 2)

# uvsw = pcares$uv*pcares$weight %*% t(pracma::pinv(pca$loadings ))    
# points( uvsw[2] ~ uvsw[1], col="green" )
# segments(0, 0, x1 = uvsw[1], y1 = uvsw[2], col = "green", lwd = 2)

```
  

## Question: Where do these predators and snow crab co-occur?

Next use eigenvectors from community field and apply to individual predators's consumption patterns ... this identifies where along the axes they are focussing in terms of diet.

Plot distances along axes of key species.

 
## Issues[^1]
 
## Bibliography  


~~~
~~~


[^1]: test

