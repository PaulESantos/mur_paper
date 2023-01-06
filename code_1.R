# Read the R script as a string
get_accum_data <- function(comm, method = "rarefaction"){
  METHODS <- c("collector", "random", "exact", "rarefaction",
               "coleman")
  method <- match.arg(method, METHODS)

  accurve <- specaccum(comm, method = method)
  tibble(sites = accurve$sites,
                richness = accurve$richness,
                sd = accurve$sd,
                method = accurve$method) |>
  mutate(min = richness - (sd * 2),
                  max = richness + (sd * 2))
}
