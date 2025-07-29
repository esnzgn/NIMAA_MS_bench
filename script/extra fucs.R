# some extra fucs
library(httr)
library(jsonlite)

library(httr)
library(jsonlite)

# Function to fetch data for a batch of UniProt IDs
fetch_uniprot_batch <- function(ids) {
  results <- lapply(ids, function(id) {
    url <- paste0("https://rest.uniprot.org/uniprotkb/", id, ".json")
    res <- GET(url)
    if (status_code(res) == 200) {
      data <- fromJSON(content(res, "text", encoding = "UTF-8"))
      return(data.frame(
        accession = id,
        organism  = data$organism$scientificName,
        taxon_id  = data$organism$taxonId,
        stringsAsFactors = FALSE
      ))
    } else {
      return(data.frame(accession = id, organism = NA, taxon_id = NA))
    }
  })
  do.call(rbind, results)
}

# Wrapper to handle large lists of IDs in smaller chunks
get_uniprot_organisms <- function(uniprot_ids, batch_size = 50) {
  all_results <- list()
  n <- length(uniprot_ids)

  for (i in seq(1, n, by = batch_size)) {
    batch_ids <- uniprot_ids[i:min(i + batch_size - 1, n)]
    cat("Fetching batch:", i, "to", i + length(batch_ids) - 1, "\n")
    batch_result <- fetch_uniprot_batch(batch_ids)
    all_results[[length(all_results) + 1]] <- batch_result
    Sys.sleep(0.5) # avoid rate-limiting
  }

  return(do.call(rbind, all_results))
}

