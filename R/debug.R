#' plot_compare_ranks
#'
#' @param res1 haystack result.
#' @param res2 haystack result.
#' @param sort_by column to sort results (default: log.p.vals).
#'
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

#' plot_rand_KLD
#'
#' Plots the distribution of randomized KLD for each of the genes, together with
#' the mean and standard deviation, the 0.95 quantile and the 0.95 quantile from
#' a normal distribution with mean and standard deviations from the distribution
#' of KLDs. The logCV is indicated in the subtitle of each plot.
#'
#' @param x haystack result.
#' @param n number of genes from randomization set to plot.
#' @param log whether to use log of KLD.
#' @param tail whether the genes are chosen from the tail of randomized genes.
#'
plot_rand_KLD <- function(x, n=12, log=TRUE, tail=FALSE) {
  kld_rand <- x$info$randomization$KLD_rand
  if (log) kld_rand <- log(kld_rand)
  logcv <- x$info$randomization$mean$observed$x

  N <- seq_len(n)
  if (tail) N <- rev(N)

  ngenes <- nrow(kld_rand)
  select_k <- function(i, tail=FALSE) {
    if (tail)
      i <- ngenes - i
    i
  }
  p <- lapply(N, function(i) {
    k <- select_k(i, tail)
    m <- mean(kld_rand[k, ])
    s <- sd(kld_rand[k, ])
    q95 <- quantile(kld_rand[k, ], 0.95)
    n95 <- qnorm(0.95, mean=m, sd=s)
    qplot(kld_rand[k, ], bins=30) +
      geom_vline(xintercept=m, color="limegreen") +
      geom_vline(xintercept=c(m-s, m+s), color="violetred") +
      geom_vline(xintercept=q95, color="steelblue1") +
      geom_vline(xintercept=n95, color="red") +
      theme(axis.text=element_blank()) +
      labs(x=k, subtitle=paste0("log10CV: ", format(logcv[k], digits=2, nsmall=2)))
  })
  
  patchwork::wrap_plots(p)
}
