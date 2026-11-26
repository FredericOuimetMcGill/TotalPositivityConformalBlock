# Strict total positivity of a Gauss hypergeometric kernel
# Appendix D.2: certify the negative conformal-block minor of order five.
#
# Requires R >= 4.3 and the CRAN package flint (install.packages("flint")).
# Tested with R 4.5.0 and flint 0.1.5.
# Run with: Rscript negative_minor_order_five.R

library(flint)

verify_minor <- function() {
    previous_precision <- flintPrec()
    on.exit(flintPrec(previous_precision))
    flintPrec(256L)

    # R fills matrices by columns. Within each column, Delta runs from
    # 21/64 to 25/64; the column points are z = 1 - 2^(-24), ..., 1 - 2^(-28).
    dimension_indices <- 21:25
    exponents <- 24:28
    delta <- arb(rep(dimension_indices, times = length(exponents))) / 64L
    t <- arb(2L)^(-rep(exponents, each = length(dimension_indices)))
    z <- 1L - t

    check_determinant <- function(entries, method) {
        dim(entries) <- c(length(dimension_indices), length(exponents))
        scaled <- arb(10L)^53L * det(entries)
        stopifnot(isTRUE(scaled > -4450L), isTRUE(scaled < -4449L))
        cat(method, ": -4450 < 10^53 * determinant < -4449\n", sep = "")
    }

    # Direct evaluation of G_delta(z) = z^delta * 2F1(delta, delta; 2*delta; z).
    values <- z^delta * arb_hypgeom_2f1(delta, delta, 2L * delta, z)
    check_determinant(values, "Hypergeometric function")

    # Independent check using the zero-balanced series. At the start of
    # iteration k, weight is ((delta)_k / k!)^2 * t^k.
    log_term <- -log(t)
    weight <- arb(1L)
    partial_sum <- arb(0L)
    for (k in 0:9) {
        psi_difference <- arb_hypgeom_polygamma(0L, k + 1L) -
                          arb_hypgeom_polygamma(0L, k + delta)
        partial_sum <- partial_sum + weight * (log_term + 2L * psi_difference)
        weight <- weight * t * ((delta + k) / (k + 1L))^2L
    }

    # The omitted tail belongs to [0, tail_bound]. The ball below is [0, 1].
    tail_bound <- (log_term + arb(1L) / 5L) * t^10L / (1L - t)
    tail_interval <- arb(mid = 0.5, rad = 0.5) * tail_bound
    values <- z^delta * (partial_sum + tail_interval) / arb_hypgeom_beta(delta, delta)
    check_determinant(values, "Ten terms and a tail bound")
}

verify_minor()
