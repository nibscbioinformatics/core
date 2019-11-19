ena_search_accession <- function(query){

  require(xml2)
  require(tidyverse)
  # https://www.ebi.ac.uk/ena/data/warehouse/search?query=%22AY585242%22&result=sequence_release&display=report
  base_url <- "https://www.ebi.ac.uk/ena/data/warehouse/search?query=%22"
  
    url <- paste0(base_url, query, "%22&result=sequence_release&display=report")
    url2 <- URLencode(url)
    
    x <- read_delim(url2, delim = "\t")

  return(x)
}


# https://www.ebi.ac.uk/ena/data/view/AY585242&display=text