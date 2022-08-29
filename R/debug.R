plot_compare_ranks <- function(res1, res2, sort_by="log.p.vals") {
  sum1 <- res1$results
  sum1 <- sum1[order(sum1[[sort_by]]), ]

  sum2 <- res2$results
  sum2 <- sum2[order(sum2[[sort_by]]), ]

  d1 <- data.frame(rank1=seq_len(nrow(sum1)), gene1=rownames(sum1))
  d2 <- data.frame(rank2=seq_len(nrow(sum2)), gene2=rownames(sum2))

  d1 <- d1[order(d1[["gene1"]]), ]
  d2 <- d2[order(d2[["gene2"]]), ]

  if (!identical(d1[["gene1"]], d2[["gene2"]]))
    stop("Features in res1 and res2 results not identical.")

  d <- cbind(d1, d2)

  ggplot(d, aes(.data[["rank1"]], .data[["rank2"]])) +
    geom_point(size=.1) +
    geom_abline(slope=1, intercept=0, color="limegreen") +
    theme(panel.grid=element_line(size=.1))
}
