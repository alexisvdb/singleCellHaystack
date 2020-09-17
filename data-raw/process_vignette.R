library(here)

unlink(here("vignettes/examples/a07_moca_100k.Rmd"))
unlink(here("vignettes/examples/figure"), recursive = TRUE)
knitr::opts_knit$set(base.dir = here("vignettes/examples"))
knitr::knit(here("data-raw/source_a07_moca_100k.Rmd"), output = here("vignettes/examples/a07_moca_100k.Rmd"))
