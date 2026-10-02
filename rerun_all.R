# =====================================================================
# Rerun of all experiments for
#   "Counterfactual Bandits for Minimum Necessary Dose Identification"
# using the REVISED Algorithms 1 (MND-B UCB) and 2 (MND-H UCB).
#
# Output (in ./results):
#   Figures  : PDF regret plots (one per figure in the paper, see FIGURES below)
#   Tables   : table4.tex, table6.tex, table7.tex  (LaTeX, ready to \input)
#              table_identification.tex          (chosen dose for every MND run)
#   Text     : ground_truth.txt  (theta values to quote in the paper)
#
# Usage:  Rscript rerun_all.R            # full run (100 reps, T = 10000)
#         Rscript rerun_all.R quick      # smoke test (3 reps, T = 500)
#
# Conventions (match the revised paper):
#   * regret is charged on the dose actually ALLOCATED (a_t), Eq. (3);
#   * regret is counted over rounds t = 1..T, excluding initialization;
#   * the chosen dose is the most frequently allocated dose over t = 1..T.
# =====================================================================

args  <- commandArgs(trailingOnly = TRUE)
QUICK <- length(args) > 0 && args[1] == "quick"
REPS  <- if (QUICK) 3 else 100
TT    <- if (QUICK) 500 else 10000
CORES <- max(1, parallel::detectCores() - 1)
OUT   <- "results"
dir.create(OUT, showWarnings = FALSE)
KAPPAS <- c(1, 0.1, 0.01)

# ---------------------------------------------------------------------
# Settings (dose-response f(k); Y = f(X) + U, U ~ N(0,1))
# ---------------------------------------------------------------------
f_main <- function(k) -k^2 - k + 1
f_quant <- function(p0, beta) function(k) qnorm(p0 - beta * k)   # P(Y_k > 0) = p0 - beta k

# ---------------------------------------------------------------------
# Ground truth
# ---------------------------------------------------------------------
truth <- function(type, K, f, yu, yl) {
  k  <- 0:K
  p1 <- 1 - pnorm(yu - f(k))            # P(Y_k > y^u)
  q  <- pnorm(yu - f(k)) - pnorm(yl - f(k))   # P(y^l < Y_k <= y^u)
  th <- numeric(K + 1)
  if (type == "B") {
    th[2:(K + 1)] <- p1[1:K] - p1[2:(K + 1)]
    th[1] <- 1 + p1[K + 1] - p1[1]
  } else {
    th[2:(K + 1)] <- pmin(p1[1:K] - p1[2:(K + 1)], q[2:(K + 1)])
    th[1] <- 1 - sum(th[2:(K + 1)])
  }
  list(theta = th, benefit = 1 - p1, harmless = q)
}

# ---------------------------------------------------------------------
# Revised Algorithm 1 (type "B") and Algorithm 2 (type "H").
# confirm = FALSE gives the "without confirmation task" baseline (a_t = k_t).
# Returns the allocated doses a_1..a_T (0-based dose labels).
# ---------------------------------------------------------------------
mnd_ucb <- function(type, K, f, yu, yl, TT, kappa, Tinit, confirm = TRUE) {
  N <- numeric(K + 1); s1 <- numeric(K + 1); s2 <- numeric(K + 1)
  pull <- function(d) {                  # d is 1-based index (dose d-1)
    y <- f(d - 1) + rnorm(1)
    N[d]  <<- N[d] + 1
    s1[d] <<- s1[d] + (y > yu)
    s2[d] <<- s2[d] + (y > yl)
  }
  for (d in 1:(K + 1)) for (i in seq_len(Tinit)) pull(d)
  alloc <- integer(TT)
  lo <- 1:K; hi <- 2:(K + 1)             # (dose k-1, dose k) for k = 1..K
  for (t in 1:TT) {
    p1 <- s1 / N; p2 <- s2 / N; lg <- log(t)
    r  <- sqrt(2 * lg * (1 / N[lo] + 1 / N[hi]))
    if (type == "B") {
      U  <- c(1 + p1[K + 1] - p1[1] + sqrt(2 * lg * (1 / N[1] + 1 / N[K + 1])),
              p1[lo] - p1[hi] + r)
    } else {
      m  <- pmin(p1[lo] - p1[hi], p2[hi] - p1[hi])
      U  <- c(1 - sum(m) + sum(r), m + r)
    }
    k <- which.max(U)                     # first maximum = smallest dose on ties
    if (k == 1) {                         # dose 0
      j <- if (type == "B") K + 1 else 1 + which.min(N[2:(K + 1)])
    } else j <- k - 1
    a <- if (confirm && t >= 2 && N[k] >= kappa / lg * (N[j] + N[k]) * N[j]) j else k
    pull(a); alloc[t] <- a - 1
  }
  alloc
}

# ---------------------------------------------------------------------
# Baselines (Algorithm 3 and MED-based Thompson sampling).
# crit = "B": success = (Y <= y^u);  crit = "H": success = (y^l < Y <= y^u)
# ---------------------------------------------------------------------
success <- function(y, crit, yu, yl) if (crit == "B") (y <= yu) else (y <= yu && y > yl)

ucb_baseline <- function(crit, K, f, yu, yl, TT, Tinit = 1) {
  N <- numeric(K + 1); S <- numeric(K + 1)
  for (d in 1:(K + 1)) for (i in seq_len(Tinit)) {
    N[d] <- N[d] + 1; S[d] <- S[d] + success(f(d - 1) + rnorm(1), crit, yu, yl)
  }
  alloc <- integer(TT)
  for (t in 1:TT) {
    a <- which.max(S / N + sqrt(2 * log(t) / N))
    N[a] <- N[a] + 1; S[a] <- S[a] + success(f(a - 1) + rnorm(1), crit, yu, yl)
    alloc[t] <- a - 1
  }
  alloc
}

# MED-TS (Filloon, 1995): smallest dose whose sampled success prob >= phi
#   (if none, the dose with the largest sample).
# TS_A (Aziz et al., 2021): dose whose sampled success prob is closest to phi.
ts_baseline <- function(rule, crit, phi, K, f, yu, yl, TT) {
  S <- numeric(K + 1); F <- numeric(K + 1); alloc <- integer(TT)
  for (t in 1:TT) {
    th <- rbeta(K + 1, 1 + S, 1 + F)
    if (rule == "MED") {
      ok <- which(th >= phi)
      a <- if (length(ok)) min(ok) else which.max(th)
    } else a <- which.min(abs(th - phi))
    s <- success(f(a - 1) + rnorm(1), crit, yu, yl)
    S[a] <- S[a] + s; F[a] <- F[a] + 1 - s; alloc[t] <- a - 1
  }
  alloc
}

# ---------------------------------------------------------------------
# Runner: REPS independent runs in parallel; returns mean cumulative regret
# (w.r.t. theta of the MND problem) and the distribution of chosen doses.
# ---------------------------------------------------------------------
run <- function(id, algo, K, theta) {
  res <- parallel::mclapply(1:REPS, function(r) {
    set.seed(100000 * id + r)
    al  <- algo()
    list(cum = cumsum((max(theta) - theta)[al + 1]),
         chosen = which.max(tabulate(al + 1, nbins = K + 1)) - 1)
  }, mc.cores = CORES)
  cum <- rowMeans(sapply(res, `[[`, "cum"))
  stopifnot(all(cum <= (1:TT) * max(max(theta) - theta) + 1e-9))  # sanity check
  ch  <- tabulate(sapply(res, `[[`, "chosen") + 1, nbins = K + 1)
  list(cum = cum, chosen = ch, mode = which.max(ch) - 1, pct = round(100 * max(ch) / REPS))
}
fmt <- function(x, kstar) {
  s <- sprintf("%d (%d\\%%)", x$mode, x$pct)
  if (x$mode == kstar) paste0("\\textbf{", x$mode, "} (", x$pct, "\\%)") else s
}
JOB <- 0; nextid <- function() { JOB <<- JOB + 1; JOB }

# ---------------------------------------------------------------------
# Plot helper (base R, same style as the paper: red solid / green dashed / blue dotted)
# ---------------------------------------------------------------------
plot_regret <- function(file, runs, ylim = NULL) {
  pdf(file.path(OUT, file), width = 5, height = 5)
  ys <- sapply(runs, `[[`, "cum")
  if (is.null(ylim)) ylim <- c(0, 1.05 * max(ys))
  plot(1:TT, ys[, 1], type = "l", col = "red", lty = 1, ylim = ylim,
       xlab = "Round", ylab = "Mean cumulative regret")
  lines(1:TT, ys[, 2], col = "green3", lty = 2)
  lines(1:TT, ys[, 3], col = "blue", lty = 3)
  legend("topleft", legend = c(expression(kappa == 1), expression(kappa == 0.1),
                               expression(kappa == 0.01)),
         col = c("red", "green3", "blue"), lty = 1:3, bty = "o")
  dev.off()
  # report the vertical order at T for the caption
  ord <- order(ys[TT, ], decreasing = TRUE)
  cat(sprintf("%-28s order at T (top->bottom): %s\n", file,
              paste(c("kappa=1", "kappa=0.1", "kappa=0.01")[ord], collapse = " > ")))
}

mnd_set <- function(type, K, f, yu, yl, Tinit, th) {
  lapply(KAPPAS, function(kap) run(nextid(),
    function() mnd_ucb(type, K, f, yu, yl, TT, kap, Tinit), K, th))
}

gt <- file(file.path(OUT, "ground_truth.txt"), "w")
logtruth <- function(name, x) writeLines(sprintf("%-22s %s", name,
  paste(sprintf("%.3f", x), collapse = ", ")), gt)
ident <- list()   # rows for table_identification.tex
addid <- function(setting, Tinit, runs, kstar) ident[[length(ident) + 1]] <<-
  sprintf("%s & %d & %s & %s & %s \\\\", setting, Tinit,
          fmt(runs[[1]], kstar), fmt(runs[[2]], kstar), fmt(runs[[3]], kstar))

# =====================================================================
# 1. Main experiment (Section 5): K = 2, y^u = 0, y^l = -1, T_init = 1
#    -> Figures 1-2, Table 4
# =====================================================================
cat("== Main experiment\n")
K <- 2; yu <- 0; yl <- -1; f <- f_main
tB <- truth("B", K, f, yu, yl); tH <- truth("H", K, f, yu, yl)
logtruth("main theta^MNB", tB$theta); logtruth("main theta^MNH", tH$theta)
logtruth("main theta^B", tB$benefit); logtruth("main theta^H", tH$harmless)
kB <- which.max(tB$theta) - 1; kH <- which.max(tH$theta) - 1

mB <- mnd_set("B", K, f, yu, yl, 1, tB$theta); plot_regret("Fig1_bene.pdf", mB)
mH <- mnd_set("H", K, f, yu, yl, 1, tH$theta); plot_regret("Fig1_harm.pdf", mH)
addid("Main, MND-B", 1, mB, kB); addid("Main, MND-H", 1, mH, kH)

baselines <- function(K, f, yu, yl, thB, thH, phiB, phiH) {
  b <- list(); h <- list()
  b[["MND-B UCB without confirmation task"]] <- run(nextid(), function()
    mnd_ucb("B", K, f, yu, yl, TT, 1, 1, confirm = FALSE), K, thB)
  h[["MND-H UCB without confirmation task"]] <- run(nextid(), function()
    mnd_ucb("H", K, f, yu, yl, TT, 1, 1, confirm = FALSE), K, thH)
  b[["Benefit UCB"]]  <- run(nextid(), function() ucb_baseline("B", K, f, yu, yl, TT), K, thB)
  h[["Harmless UCB"]] <- run(nextid(), function() ucb_baseline("H", K, f, yu, yl, TT), K, thH)
  for (rule in c("MED", "TSA")) {
    nm <- if (rule == "MED") "MED-TS" else "TS\\_A"
    for (phi in phiB) b[[sprintf("%s algorithm for benefit ($\\phi=%s$)", nm, phi)]] <-
      run(nextid(), function() ts_baseline(rule, "B", phi, K, f, yu, yl, TT), K, thB)
    for (phi in phiH) h[[sprintf("%s algorithm for harmless ($\\phi=%s$)", nm, phi)]] <-
      run(nextid(), function() ts_baseline(rule, "H", phi, K, f, yu, yl, TT), K, thH)
  }
  list(b = b, h = h)
}
base2 <- baselines(K, f, yu, yl, tB$theta, tH$theta, c(0.1, 0.8, 0.95), c(0.05, 0.1, 0.3))

rows <- function(m, base, label, kstar) c(
  sprintf("MND-%s UCB ($\\kappa=%s$) & %s", label, KAPPAS, sapply(m, fmt, kstar = kstar)),
  "\\hline",
  sprintf("%s & %s", names(base), sapply(base, fmt, kstar = kstar)))
rb <- rows(mB, base2$b, "B", kB); rh <- rows(mH, base2$h, "H", kH)
writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  sprintf(paste0("\\caption{The most frequently chosen dose with $T=%d$ by each benefit and harmless ",
    "algorithm across %d repeated runs and the proportion of runs in which it is chosen. ",
    "The MND-B is %d and the MND-H is %d.}"), TT, REPS, kB, kH),
  "\\label{tab:4}", "\\small", "\\begin{tabular}{lc|lc}", "\\hline",
  "Algorithm (Benefit) & Chosen dose & Algorithm (Harmless) & Chosen dose \\\\", "\\hline",
  ifelse(rb == "\\hline", "\\hline", paste0(rb, " & ", rh, " \\\\")),
  "\\hline", "\\end{tabular}", "\\end{table}"), file.path(OUT, "table4.tex"))

# =====================================================================
# 2. Appendix E.1.1: larger T_init (Figures 3-8)
# =====================================================================
cat("== E.1.1\n")
for (Ti in c(10, 100, 1000)) {
  r1 <- mnd_set("B", K, f, yu, yl, Ti, tB$theta)
  plot_regret(sprintf("FigE11_bene_Tinit%d.pdf", Ti), r1); addid("Main, MND-B", Ti, r1, kB)
  r2 <- mnd_set("H", K, f, yu, yl, Ti, tH$theta)
  plot_regret(sprintf("FigE11_harm_Tinit%d.pdf", Ti), r2); addid("Main, MND-H", Ti, r2, kH)
}

# =====================================================================
# 3. Appendix E.2.1: MND-B with y^u = 0.5 (Figures 9-11)
# =====================================================================
cat("== E.2.1\n")
t21 <- truth("B", 2, f_main, 0.5, -1); logtruth("E.2.1 theta^MNB", t21$theta)
for (Ti in c(1, 10, 100)) {
  r <- mnd_set("B", 2, f_main, 0.5, -1, Ti, t21$theta)
  plot_regret(sprintf("FigE21_bene_Tinit%d.pdf", Ti), r)
  addid("E.2.1, MND-B", Ti, r, which.max(t21$theta) - 1)
}

# =====================================================================
# 4. Appendix E.2.2: MND-H with y^l = -2 (Figures 12-14) + no-confirmation baseline
# =====================================================================
cat("== E.2.2\n")
t22 <- truth("H", 2, f_main, 0, -2); logtruth("E.2.2 theta^MNH", t22$theta)
for (Ti in c(1, 10, 100)) {
  r <- mnd_set("H", 2, f_main, 0, -2, Ti, t22$theta)
  plot_regret(sprintf("FigE22_harm_Tinit%d.pdf", Ti), r)
  addid("E.2.2, MND-H", Ti, r, which.max(t22$theta) - 1)
}
nc <- run(nextid(), function() mnd_ucb("H", 2, f_main, 0, -2, TT, 1, 1, confirm = FALSE),
          2, t22$theta)
writeLines(sprintf("E.2.2 MND-H UCB without confirmation (T_init=1): chosen dose %d (%d%%)",
                   nc$mode, nc$pct), gt)

# =====================================================================
# 5. Appendix E.2.3: K = 3 (Tables 6-7, Figures 15-16)
# =====================================================================
cat("== E.2.3 (K = 3)\n")
t3B <- truth("B", 3, f_main, 0, -1); t3H <- truth("H", 3, f_main, 0, -1)
logtruth("K=3 theta^MNB", t3B$theta); logtruth("K=3 theta^MNH", t3H$theta)
k3B <- which.max(t3B$theta) - 1; k3H <- which.max(t3H$theta) - 1
m3B <- mnd_set("B", 3, f_main, 0, -1, 1, t3B$theta); plot_regret("FigE23_bene.pdf", m3B)
m3H <- mnd_set("H", 3, f_main, 0, -1, 1, t3H$theta); plot_regret("FigE23_harm.pdf", m3H)
base3 <- baselines(3, f_main, 0, -1, t3B$theta, t3H$theta, c(0.1, 0.8, 0.95), c(0.05, 0.1, 0.3))
onetable <- function(file, label, m, base, kstar, crit) writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  sprintf(paste0("\\caption{The most frequently chosen dose with $T=%d$ by each %s algorithm across ",
    "%d repeated runs ($K=3$) and the proportion of runs in which it is chosen. The MND-%s is %d.}"),
    TT, crit, REPS, label, kstar),
  sprintf("\\label{tab:%s}", sub(".tex", "", file)),
  "\\begin{tabular}{lc}", "\\hline", "Algorithm & Chosen dose \\\\", "\\hline",
  paste(rows(m, base, label, kstar)[-4], "\\\\")[1:3], "\\hline",
  paste(rows(m, base, label, kstar)[-(1:4)], "\\\\"),
  "\\hline", "\\end{tabular}", "\\end{table}"), file.path(OUT, file))
onetable("table6.tex", "B", m3B, base3$b, k3B, "benefit")
onetable("table7.tex", "H", m3H, base3$h, k3H, "harmless")
addid("$K=3$, MND-B", 1, m3B, k3B); addid("$K=3$, MND-H", 1, m3H, k3H)

# =====================================================================
# 6. Appendix E.2.4 / E.2.5: K = 10 and K = 20, T_init = 10 (Figures 17-20)
# =====================================================================
for (cfg in list(list(K = 10, bB = 0.034, bH = 0.037), list(K = 20, bB = 0.019, bH = 0.018))) {
  cat(sprintf("== K = %d\n", cfg$K))
  fB <- f_quant(0.4, cfg$bB); fH <- f_quant(0.4, cfg$bH)
  tKB <- truth("B", cfg$K, fB, 0, -2); tKH <- truth("H", cfg$K, fH, 0, -2)
  logtruth(sprintf("K=%d theta^MNB", cfg$K), tKB$theta)
  logtruth(sprintf("K=%d theta^MNH", cfg$K), tKH$theta)
  rB <- mnd_set("B", cfg$K, fB, 0, -2, 10, tKB$theta)
  rH <- mnd_set("H", cfg$K, fH, 0, -2, 10, tKH$theta)
  plot_regret(sprintf("FigE2_K%d_bene.pdf", cfg$K), rB)
  plot_regret(sprintf("FigE2_K%d_harm.pdf", cfg$K), rH)
  addid(sprintf("$K=%d$, MND-B", cfg$K), 10, rB, which.max(tKB$theta) - 1)
  addid(sprintf("$K=%d$, MND-H", cfg$K), 10, rH, which.max(tKH$theta) - 1)
}

# =====================================================================
# 7. Identification summary for every MND run (useful for Appendix E text)
# =====================================================================
writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  sprintf(paste0("\\caption{The most frequently chosen dose of MND-B UCB and MND-H UCB with $T=%d$ ",
    "across %d repeated runs and the proportion of runs in which it is chosen. ",
    "Bold: the chosen dose equals the MND.}"), TT, REPS),
  "\\label{tab:ident}", "\\small", "\\begin{tabular}{lcccc}", "\\hline",
  "Setting & $T^{\\text{initial}}$ & $\\kappa=1$ & $\\kappa=0.1$ & $\\kappa=0.01$ \\\\", "\\hline",
  unlist(ident), "\\hline", "\\end{tabular}", "\\end{table}"),
  file.path(OUT, "table_identification.tex"))
close(gt)
cat("Done. Results in ./", OUT, "\n", sep = "")
