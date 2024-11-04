
  observer.historical.data = function() {


      loc = file.path( project.datadirectory("bio.snowcrab"), "data", "observer", "archive", "2003", "Zone 24", "Port Sample", "Semaine 25", "24P03053.txt" )
      write.table( scan( url( loc),  what="", sep="\n", strip.white=F ), file="tmp.out",  row.names=F, col.names=F, quote=F, sep=";" )
      dataformat = observer.data.formats()  # see auxillary function, below
      rawdata = read.fwf( file="tmp.out", widths=dataformat$w, as.is=T, colClasses="character", header=F, strip.white=F)
      colnames(rawdata) = dataformat$v
      rawdata = rawdata[, dataformat$keep]
      str(rawdata)

  rawdata = concatenate.historical.data( type="seasamples",   outfile="odb.historical.datadump.seasamples.rdata" )
  rawdata =  concatenate.historical.data( type="portsamples", outfile="odb.historical.datadump.portsamples.rdata" )


### TODO incomplete


  }


