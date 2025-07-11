
  get.marked = function(DS="file" ) {

    tags.datadir= file.path( project.datadirectory("bio.snowcrab"), "data", "tagging" )
    marked.file="tags.1996_2001.csv"
    marked = NULL

    if (DS=="redo") {
      marked = read.table( file.path( tags.datadir, marked.file), sep=";", header=T, as.is=T)
      marked$Ncrabs = NULL
      f = which(marked$lon>-55)
      marked[f,] = NA
      marked$timestamp = lubridate::mdy( marked$date)
      read_write_fast(data=marked, fn=file.path(tags.datadir, "marked.rdz") )
    }
    if (DS=="file") marked = read_write_fast( file.path(tags.datadir, "marked.rdz"))

    return (marked)

  }


