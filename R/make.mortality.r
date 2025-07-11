
make.mortality = function(p, redo=F) {
  # est total mortality
  outfilename = file.path( p$annual.results, "mortality.rates.rdz" )
  if (!redo) {
    
    return ( read_write_fast( outfilename ) )
  }

  read_write_fast( data=xx, fn=outfilename )
  return(xx)
} 


