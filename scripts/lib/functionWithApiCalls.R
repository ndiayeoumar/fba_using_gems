library(httr)
library(jsonlite)


#' @param reaction_ids vector of reaction ids
#' @return descriptions: a vector of string values of the reactions complete names
getApiReactionDesc <- function(reaction_ids) {
  # get all available reactions using VMH API
  api_res <- GET('https://www.vmh.life/_api/reactions/?page_size=19313')
  # retrieve only the real the content of the API response
  all_reactions <- fromJSON(rawToChar(api_res$content))
  # clean reaction ids given by compass output
  reaction_ids <- gsub(x = reaction_ids, pattern = "_pos|_neg", replacement = "")
  # extract complete names of the reaction in the description field
  descriptions <- all_reactions$results$description[match(
    reaction_ids,
    all_reactions$results$abbreviation)]
  return(descriptions)
}


#' @param reaction_id vector of reaction ids
#' @return genes data: a vector of string values of the reactions complete names
getApiReactionGene <- function(reaction_id) {
  # clean reaction ids given by compass output
  reaction_id <- gsub(x = reaction_id, pattern = "_pos|_neg", replacement = "")
  ## default uri
  uri <- paste0("https://www.vmh.life/_api/reactiongenes/", reaction_id)
  # get all available reactions using VMH API
  api_res <- GET(uri)
  api_res <- fromJSON(rawToChar(api_res$content))
  return(api_res)
}
