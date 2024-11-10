# Check that sample labels match metadata rows
samples_match_metadata <- function(samples.vector, metadata) {
  if (all(samples.vector==metadata$SampleName)) { return(TRUE) }
  else { return(FALSE) }
}