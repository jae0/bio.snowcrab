## Archved snippets / figures /tables, etc.




## To do next:

Estimate via carstm: requires a Poisson model of each species of catch (number) with offset of landings and covariates ... soon, ever?


##### Entanglements of megafauna

```{r}
#| label: tbl-bycatch_entanglements_all
#| warning: false
#| error: false 
#| tbl-cap: "Entanglements of megafauna in Maritimes"
#| eval: true
#| output: true

region = "cfaall"

# observer data .. extra care needed as there are duplicated records, etc
oss = observer.db( DS="bycatch_clean_data", p=p,  yrs=yrs_observer )  # Prepare at sea observed data
i = polygon_inside( oss[,  c("lon", "lat")], region=region )
oss =  oss[i,]
 
whales = oss[grep("whale", common, ignore.case=TRUE), .(whales=.N), by=.(yr)] 
leatherback = oss[grep("leatherback",  common, ignore.case=TRUE), .(leatherback=.N), by=.(yr)] 
basking_shark = oss[ grep("BASKING SHARK",  common, ignore.case=TRUE), .(basking_shark=.N), by=.(yr)] 

out = data.table( yr=yrs_observer )
out = whales[out, on="yr"]
out = leatherback[out, on="yr"]
out = basking_shark[out, on="yr"]
out[ is.na(out) ] = 0

colnames(out) = c("Year", "Whale", "Leatherback turtle", "Basking shark")

gt::gt(out) |> gt::tab_options(table.font.size = 12, data_row.padding = gt::px(1), 
  summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
  footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
  row_group.padding = gt::px(1))
  
```


##### Map of locations of megafauna entanglements:

#```{r}
#| label: fig-bycatch_entanglements_map_all
#| warning: false
#| error: false 
#| fig-cap: "Entanglement locations in Maritimes of megafauna since 2000. Whales (red), Leatherback turtles (green), Basking shark (blue)."
#| fig-dpi: 144
#| fig-height: 5

loc = file.path( data_loc, "output", "maps", "observer.entanglements" )
fn = file.path( loc, "observed_bycatch_entanglements.png" )

include_graphics( fn ) 
   
```
