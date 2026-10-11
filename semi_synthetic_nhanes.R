# =====================================================================
# Semi-synthetic insulin-initiation experiment (Appendix E.3) for
#   "Counterfactual Bandits for Minimum Necessary Dose Identification"
#
# Real patients, simulated dose response.
#   * Patients are resampled from a real cohort: insulin-naive adults with
#     diagnosed diabetes and HbA1c >= 7% in NHANES 2017-2018
#     (file nhanes_cohort.csv; see prepare_nhanes_cohort.R). Each patient
#     keeps their real baseline HbA1c and real body weight.
#   * The response to a daily basal insulin dose d (U/day) is simulated with
#     an Emax model on the excess HbA1c above 5%:
#         Y_d = 5 + (A0 - 5) * (1 - Emax * S d / (ED50 + S d)) + e,
#     where S = (median weight / weight) * exp(eta) makes the effect depend on
#     the dose per kilogram, eta is unexplained between-patient variation,
#     and e is assay noise shared across doses.
#   * Benefit: HbA1c <= 7.0% (ADA target). Harm: HbA1c <= 6.0% (the
#     intensive target that increased mortality and severe hypoglycemia in
#     the ACCORD trial).
# Because the simulator returns ALL potential outcomes of each patient, we
# can also measure counterfactual quantities (who was over-treated).
#
# Output (./results_nhanes):
#   FigNHANES_bene.pdf, FigNHANES_harm.pdf  regret of MND-B / MND-H UCB
#   table_nhanes_cohort.tex                 real cohort characteristics
#   table_nhanes_truth.tex                  per-dose ground truth
#   table_nhanes_choice.tex                 chosen dose of every algorithm
#   table_nhanes_clinical.tex               in-trial clinical outcomes
#   appendix_E3.tex                         appendix text with numbers filled in
#   nhanes_summary.txt                      all numbers in plain text
#
# Usage:  Rscript semi_synthetic_nhanes.R          # full run (100 reps, T = 10000)
#         Rscript semi_synthetic_nhanes.R quick    # smoke test
# Requires only base R. nhanes_cohort.csv must be in the working directory.
# =====================================================================

args  <- commandArgs(trailingOnly = TRUE)
QUICK <- length(args) > 0 && args[1] == "quick"
REPS  <- if (QUICK) 4 else 100
TT    <- if (QUICK) 1000 else 10000
TINIT <- 10
CORES <- max(1, parallel::detectCores() - 1)
OUT   <- "results_nhanes"
dir.create(OUT, showWarnings = FALSE)
KAPPAS <- c(1, 0.1, 0.01)
KAPPAS_H_SMALL <- c(0.003, 0.001)            # extra values for MND-H (see appendix text)
PHIS   <- c(0.5, 0.7, 0.8)                   # MED thresholds within the attainable range

# ---------------------------------------------------------------------
# Real cohort and simulator
# ---------------------------------------------------------------------
cohort <- read.csv("nhanes_cohort.csv")
stopifnot(all(c("hba1c", "weight_kg") %in% names(cohort)))
A1C <- cohort$hba1c; WT <- cohort$weight_kg; WMED <- median(WT)

DOSES <- c(0, 10, 20, 30, 40)                # basal insulin, U/day; dose index k = 0..4
K     <- length(DOSES) - 1
YU <- 7.0; YL <- 6.0                         # HbA1c (%): ADA target / ACCORD intensive target
PAR <- list(floor = 5, Emax = 0.7, ED50 = 10, eta = 0.3, sd = 0.15)

# n x (K+1) matrix of potential outcomes Y_0..Y_K for n patients resampled from the cohort
patients <- function(n) {
  i  <- sample.int(length(A1C), n, replace = TRUE)
  A0 <- A1C[i]
  S  <- (WMED / WT[i]) * exp(rnorm(n, 0, PAR$eta))   # dose-per-kg effect x residual variation
  e  <- rnorm(n, 0, PAR$sd)                         # assay noise, shared across doses
  sapply(DOSES, function(d)
    PAR$floor + (A0 - PAR$floor) * (1 - PAR$Emax * S * d / (PAR$ED50 + S * d)) + e)
}

# ---------------------------------------------------------------------
# Ground truth from 2,000,000 simulated patients
# ---------------------------------------------------------------------
set.seed(2027)
Ybig <- patients(if (QUICK) 2e5 else 2e6)
p1 <- colMeans(Ybig > YU); p2 <- colMeans(Ybig > YL)
lo <- 1:K; hi <- 2:(K + 1)
thB   <- 1 - p1                                          # P(HbA1c <= 7)
thH   <- p2 - p1                                         # P(6 < HbA1c <= 7)
thMNB <- c(1 + p1[K + 1] - p1[1], p1[lo] - p1[hi])       # Theorem 4
idMNH <- c(NA, pmin(p1[lo] - p1[hi], p2[hi] - p1[hi]))   # Theorem 8 (identified)
idMNH[1] <- 1 - sum(idMNH[-1])
okH <- (Ybig[, lo] > YU) & (Ybig[, hi] <= YU) & (Ybig[, hi] > YL)   # true counterfactual
indH <- ifelse(rowSums(okH) > 0, max.col(okH, ties.method = "first"), 0)
thMNH <- tabulate(indH + 1, nbins = K + 1) / nrow(Ybig)
intens <- colMeans(Ybig <= YL)                           # over-intensive control
over   <- c(0, colMeans(Ybig[, lo] <= YU))               # target also reached at next lower dose
A3err  <- max(abs(idMNH - thMNH))                        # Assumption 3 check
red    <- colMeans(Ybig[, 1] - Ybig[, -1])               # mean HbA1c reduction vs. no insulin
rm(Ybig, okH, indH); invisible(gc())

kB <- which.max(thMNB) - 1; kH <- which.max(thMNH) - 1
gapB <- sort(thMNB, TRUE)[1] - sort(thMNB, TRUE)[2]
gapH <- sort(thMNH, TRUE)[1] - sort(thMNH, TRUE)[2]

# ---------------------------------------------------------------------
# Algorithms. Each run gets a matrix Y of potential outcomes, one row per
# arriving patient: rows 1..(K+1)*TINIT are used for initialization and
# the next TT rows are the trial patients. Returns allocated dose indices.
# ---------------------------------------------------------------------
NINIT <- (K + 1) * TINIT

mnd_ucb <- function(type, Y, kappa, confirm = TRUE) {
  N <- numeric(K + 1); s1 <- numeric(K + 1); s2 <- numeric(K + 1); row <- 0
  pull <- function(a) { row <<- row + 1; y <- Y[row, a]
    N[a] <<- N[a] + 1; s1[a] <<- s1[a] + (y > YU); s2[a] <<- s2[a] + (y > YL) }
  for (a in 1:(K + 1)) for (i in seq_len(TINIT)) pull(a)
  alloc <- integer(TT)
  for (t in 1:TT) {
    q1 <- s1 / N; q2 <- s2 / N; lg <- log(t)
    r <- sqrt(2 * lg * (1 / N[lo] + 1 / N[hi]))
    if (type == "B") {
      U <- c(1 + q1[K + 1] - q1[1] + sqrt(2 * lg * (1 / N[1] + 1 / N[K + 1])), q1[lo] - q1[hi] + r)
    } else {
      m <- pmin(q1[lo] - q1[hi], q2[hi] - q1[hi]); U <- c(1 - sum(m) + sum(r), m + r)
    }
    k <- which.max(U)
    j <- if (k > 1) k - 1 else if (type == "B") K + 1 else 1 + which.min(N[2:(K + 1)])
    a <- if (confirm && t >= 2 && N[k] >= kappa / lg * (N[j] + N[k]) * N[j]) j else k
    pull(a); alloc[t] <- a
  }
  alloc - 1
}

success <- function(y, crit) if (crit == "B") y <= YU else (y <= YU && y > YL)

ucb_baseline <- function(crit, Y) {
  N <- numeric(K + 1); S <- numeric(K + 1); row <- 0
  for (a in 1:(K + 1)) for (i in seq_len(TINIT)) {
    row <- row + 1; N[a] <- N[a] + 1; S[a] <- S[a] + success(Y[row, a], crit) }
  alloc <- integer(TT)
  for (t in 1:TT) {
    a <- which.max(S / N + sqrt(2 * log(t) / N)); row <- row + 1
    N[a] <- N[a] + 1; S[a] <- S[a] + success(Y[row, a], crit); alloc[t] <- a
  }
  alloc - 1
}

# MED-TS (Filloon, 1995): smallest dose whose posterior sample >= phi
# TS_A   (Aziz et al., 2021): dose whose posterior sample is closest to phi
ts_baseline <- function(rule, crit, phi, Y) {
  S <- numeric(K + 1); F <- numeric(K + 1); alloc <- integer(TT); row <- NINIT
  for (t in 1:TT) {
    th <- rbeta(K + 1, 1 + S, 1 + F)
    a <- if (rule == "MED") { ok <- which(th >= phi); if (length(ok)) min(ok) else which.max(th) }
         else which.min(abs(th - phi))
    row <- row + 1; s <- success(Y[row, a], crit)
    S[a] <- S[a] + s; F[a] <- F[a] + 1 - s; alloc[t] <- a
  }
  alloc - 1
}

# ---------------------------------------------------------------------
# Runner: regret curve, chosen dose and in-trial clinical outcomes
# ---------------------------------------------------------------------
JOB <- 0
run <- function(algo, theta) {
  JOB <<- JOB + 1; id <- JOB
  res <- parallel::mclapply(1:REPS, function(r) {
    set.seed(100000 * id + r)
    Y  <- patients(NINIT + TT)
    al <- algo(Y)
    Yt <- Y[NINIT + 1:TT, , drop = FALSE]
    ya <- Yt[cbind(1:TT, al + 1)]
    yl <- ifelse(al >= 1, Yt[cbind(1:TT, pmax(al, 1))], Inf)
    list(cum = cumsum((max(theta) - theta)[al + 1]),
         chosen = which.max(tabulate(al + 1, nbins = K + 1)) - 1,
         clin = c(units = mean(DOSES[al + 1]), cured = mean(ya <= YU),
                  intens = mean(ya <= YL), over = mean(al >= 1 & yl <= YU)))
  }, mc.cores = CORES)
  ch <- tabulate(sapply(res, `[[`, "chosen") + 1, nbins = K + 1)
  list(cum = rowMeans(sapply(res, `[[`, "cum")), chosen = ch,
       mode = which.max(ch) - 1, pct = round(100 * max(ch) / REPS),
       clin = rowMeans(sapply(res, `[[`, "clin")))
}
lab <- function(x, kstar) {
  s <- sprintf("%d U/day (%d\\%%)", DOSES[x$mode + 1], x$pct)
  if (x$mode == kstar) paste0("\\textbf{", DOSES[x$mode + 1], " U/day} (", x$pct, "\\%)") else s
}

cat("Running MND-B / MND-H UCB ...\n")
mB  <- lapply(KAPPAS, function(kap) run(function(Y) mnd_ucb("B", Y, kap), thMNB))
mH  <- lapply(KAPPAS, function(kap) run(function(Y) mnd_ucb("H", Y, kap), thMNH))
mHs <- lapply(KAPPAS_H_SMALL, function(kap) run(function(Y) mnd_ucb("H", Y, kap), thMNH))
cat("Running baselines ...\n")
bB <- list(); bH <- list()
bB[["MND-B UCB without confirmation task"]] <- run(function(Y) mnd_ucb("B", Y, 1, FALSE), thMNB)
bH[["MND-H UCB without confirmation task"]] <- run(function(Y) mnd_ucb("H", Y, 1, FALSE), thMNH)
bB[["Benefit UCB"]]  <- run(function(Y) ucb_baseline("B", Y), thMNB)
bH[["Harmless UCB"]] <- run(function(Y) ucb_baseline("H", Y), thMNH)
for (rule in c("MED", "TSA")) for (phi in PHIS) {
  nm <- if (rule == "MED") "MED-TS" else "TS\\_A"
  bB[[sprintf("%s for benefit ($\\phi=%s$)", nm, phi)]] <-
    run((function(rule, phi) function(Y) ts_baseline(rule, "B", phi, Y))(rule, phi), thMNB)
  bH[[sprintf("%s for harmless ($\\phi=%s$)", nm, phi)]] <-
    run((function(rule, phi) function(Y) ts_baseline(rule, "H", phi, Y))(rule, phi), thMNH)
}

# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------
plot_regret <- function(file, runs) {
  pdf(file.path(OUT, file), width = 5, height = 5)
  ys <- sapply(runs, `[[`, "cum")
  plot(1:TT, ys[, 1], type = "l", col = "red", lty = 1, ylim = c(0, 1.05 * max(ys)),
       xlab = "Round", ylab = "Mean cumulative regret")
  lines(1:TT, ys[, 2], col = "green3", lty = 2); lines(1:TT, ys[, 3], col = "blue", lty = 3)
  legend("topleft", legend = c(expression(kappa == 1), expression(kappa == 0.1),
         expression(kappa == 0.01)), col = c("red", "green3", "blue"), lty = 1:3)
  dev.off()
}
plot_regret("FigNHANES_bene.pdf", mB); plot_regret("FigNHANES_harm.pdf", mH)

# ---------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
boldmax <- function(x) { v <- f3(x); m <- which.max(x); v[m] <- paste0("\\textbf{", v[m], "}"); v }
trow <- function(label, x, bold = FALSE) paste0(label, " & ", paste(if (bold) boldmax(x) else f3(x), collapse = " & "), " \\\\")
qs <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, .25), quantile(x, .75))

writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  sprintf("\\caption{Real cohort used in the semi-synthetic experiment: insulin-naive adults with diagnosed diabetes and HbA1c $\\ge 7\\%%$ in NHANES 2017--2018 ($n=%d$). Median [interquartile range].}", nrow(cohort)),
  "\\label{tab:nhanes_cohort}", "\\small", "\\begin{tabular}{lc}", "\\hline",
  sprintf("Age (years) & %s \\\\", qs(cohort$age)),
  sprintf("Body weight (kg) & %s \\\\", qs(WT)),
  sprintf("Baseline HbA1c (\\%%) & %s \\\\", qs(A1C)),
  "\\hline", "\\end{tabular}", "\\end{table}"), file.path(OUT, "table_nhanes_cohort.tex"))

writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  "\\caption{Ground truth of each basal insulin dose in the semi-synthetic experiment, computed from $2\\times10^6$ patients resampled from the NHANES cohort. Over-treated: the patient would also have reached the target at the next lower dose. Bold: the optimum of each criterion.}",
  "\\label{tab:nhanes_truth}", "\\small", "\\begin{tabular}{lccccc}", "\\hline",
  paste0("Dose (U/day) & ", paste(DOSES, collapse = " & "), " \\\\"), "\\hline",
  trow("$\\theta^{\\text{B}}$ (HbA1c $\\le7\\%$)", thB, TRUE),
  trow("$\\theta^{\\text{H}}$ ($6\\%<$ HbA1c $\\le7\\%$)", thH, TRUE),
  trow("$\\theta^{\\text{MNB}}$", thMNB, TRUE),
  trow("$\\theta^{\\text{MNH}}$", thMNH, TRUE),
  trow("HbA1c $\\le6\\%$", intens),
  trow("Over-treated", over),
  "\\hline", "\\end{tabular}", "\\end{table}"), file.path(OUT, "table_nhanes_truth.tex"))

rowsB <- c(sprintf("MND-B UCB ($\\kappa=%s$) & %s", KAPPAS, sapply(mB, lab, kstar = kB)),
           "\\hline", sprintf("%s & %s", names(bB), sapply(bB, lab, kstar = kB)))
rowsH <- c(sprintf("MND-H UCB ($\\kappa=%s$) & %s", KAPPAS, sapply(mH, lab, kstar = kH)),
           "\\hline", sprintf("%s & %s", names(bH), sapply(bH, lab, kstar = kH)))
writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  sprintf("\\caption{Most frequently chosen dose in the semi-synthetic experiment ($T=%d$, $T^{\\text{initial}}=%d$, %d runs) and the proportion of runs in which it is chosen. The MND-B is %d U/day and the MND-H is %d U/day.}", TT, TINIT, REPS, DOSES[kB + 1], DOSES[kH + 1]),
  "\\label{tab:nhanes_choice}", "\\small", "\\resizebox{\\textwidth}{!}{%", "\\begin{tabular}{lc|lc}", "\\hline",
  "Algorithm (Benefit) & Chosen dose & Algorithm (Harmless) & Chosen dose \\\\", "\\hline",
  ifelse(rowsB == "\\hline", "\\hline", paste0(rowsB, " & ", rowsH, " \\\\")),
  "\\hline", "\\end{tabular}}", "\\end{table}"), file.path(OUT, "table_nhanes_choice.tex"))

key <- c(list("MND-B UCB ($\\kappa=0.01$)" = mB[[3]], "MND-H UCB ($\\kappa=0.01$)" = mH[[3]]),
         setNames(mHs, sprintf("MND-H UCB ($\\kappa=%s$)", KAPPAS_H_SMALL)),
         list("Benefit UCB" = bB[["Benefit UCB"]], "Harmless UCB" = bH[["Harmless UCB"]],
              "MED-TS for benefit ($\\phi=0.7$)" = bB[["MED-TS for benefit ($\\phi=0.7$)"]],
              "TS\\_A for benefit ($\\phi=0.7$)" = bB[["TS\\_A for benefit ($\\phi=0.7$)"]]))
pc <- function(x) sprintf("%.1f\\%%", 100 * x)
writeLines(c(
  "\\begin{table}[tb]", "\\centering",
  sprintf("\\caption{In-trial outcomes over the $T=%d$ trial patients, averaged over %d runs: mean basal insulin dose, and the proportions of patients who reached the target (HbA1c $\\le7\\%%$), were driven to over-intensive control (HbA1c $\\le6\\%%$), and were over-treated (would also have reached the target at the next lower dose).}", TT, REPS),
  "\\label{tab:nhanes_clinical}", "\\small", "\\begin{tabular}{lcccc}", "\\hline",
  "Algorithm & Mean dose (U/day) & HbA1c $\\le7\\%$ & HbA1c $\\le6\\%$ & Over-treated \\\\", "\\hline",
  sapply(names(key), function(n) { c <- key[[n]]$clin
    sprintf("%s & %.1f & %s & %s & %s \\\\", n, c["units"], pc(c["cured"]), pc(c["intens"]), pc(c["over"])) }),
  "\\hline", "\\end{tabular}", "\\end{table}"), file.path(OUT, "table_nhanes_clinical.tex"))

# ---------------------------------------------------------------------
# Plain-text summary
# ---------------------------------------------------------------------
sink(file.path(OUT, "nhanes_summary.txt"))
cat(sprintf("Cohort: n = %d, median HbA1c %.1f, median weight %.1f\n", nrow(cohort), median(A1C), WMED))
cat("Doses (U/day):", DOSES, "\n")
for (nm in c("thB", "thH", "thMNB", "thMNH", "idMNH", "intens", "over"))
  cat(sprintf("%-6s %s\n", nm, paste(f3(get(nm)), collapse = " ")))
cat("Mean HbA1c reduction by dose:", sprintf("%.2f", red), "\n")
cat(sprintf("Assumption 3 check: max |identified - true theta^MNH| = %.2e\n", A3err))
cat(sprintf("MND-B = %d U/day (gap %.3f), MND-H = %d U/day (gap %.3f)\n", DOSES[kB + 1], gapB, DOSES[kH + 1], gapH))
cat(sprintf("Theorem 9 kappa bound: %.5f; Theorem 10 kappa bound: %.6f\n", gapB^2 / 32, gapH^2 / (32 * K^2)))
show <- function(nm, x) cat(sprintf("%-45s chosen %s  final regret %7.1f  dose %.1f target %.3f intens %.3f over %.3f\n",
  nm, paste(x$chosen, collapse = "/"), x$cum[TT], x$clin["units"], x$clin["cured"], x$clin["intens"], x$clin["over"]))
for (i in seq_along(KAPPAS)) show(sprintf("MND-B kappa=%s", KAPPAS[i]), mB[[i]])
for (i in seq_along(KAPPAS)) show(sprintf("MND-H kappa=%s", KAPPAS[i]), mH[[i]])
for (i in seq_along(KAPPAS_H_SMALL)) show(sprintf("MND-H kappa=%s", KAPPAS_H_SMALL[i]), mHs[[i]])
for (n in names(bB)) show(n, bB[[n]]); for (n in names(bH)) show(n, bH[[n]])
sink()

# ---------------------------------------------------------------------
# Appendix text with the numbers filled in
# ---------------------------------------------------------------------
pct  <- function(x) sprintf("%.1f\\%%", 100 * x)
idr  <- function(x, k) round(100 * x$chosen[k + 1] / REPS)
du   <- function(x) DOSES[x$mode + 1]
sciTex <- function(x) { if (x < 1e-12) return("0"); e <- floor(log10(x)); sprintf("%.1f\\times10^{%d}", x / 10^e, e) }
mb <- mB[[3]]$clin; bu <- bB[["Benefit UCB"]]$clin
D0 <- max(thMNH) - thMNH[1]; kapStar <- (D0 / (2 * K))^2; kapStarB <- (max(thMNB) - thMNB[1])^2 / 2
iB  <- which.max(thB); iH <- which.max(thH)
smallTxt <- paste(sprintf("%d\\%% with $\\kappa=%s$", sapply(mHs, idr, k = kH), KAPPAS_H_SMALL), collapse = " and ")
mndHokSmall <- any(sapply(mHs, idr, k = kH) >= 90)

writeLines(c(
  "\\subsection{Semi-synthetic experiment with NHANES patients}\\label{appB3}",
  "\\noindent{\\bf Data and simulator.}",
  sprintf("We build a semi-synthetic benchmark in which the patients are real and only their responses to insulin are simulated. We take the %d insulin-naive adults with diagnosed diabetes and HbA1c $\\ge7\\%%$ in the National Health and Nutrition Examination Survey (NHANES) 2017--2018 \\citep{nhanes2017}, i.e., the patients for whom basal insulin initiation is typically considered (Table~\\ref{tab:nhanes_cohort}). Each simulated patient is drawn with replacement from this cohort and keeps their measured baseline HbA1c $A_0$ and body weight $W$.", nrow(cohort)),
  sprintf("For a daily basal insulin dose $d_k\\in\\{%s\\}$ U/day ($K=%d$), the HbA1c after treatment follows an $E_{\\max}$ model on the excess HbA1c above $%g\\%%$,", paste(DOSES, collapse = ", "), K, PAR$floor),
  sprintf("\\[Y_k=%g+(A_0-%g)\\Big(1-E_{\\max}\\frac{S\\,d_k}{ED_{50}+S\\,d_k}\\Big)+\\varepsilon,\\qquad S=\\frac{\\bar W}{W}\\,e^{\\eta},\\]", PAR$floor, PAR$floor),
  sprintf("with $E_{\\max}=%.1f$, $ED_{50}=%d$ U/day, $\\bar W=%.1f$\\,kg the median cohort weight, so that the effect depends on the dose per kilogram, $\\eta\\sim N(0,%.1f^2)$ unexplained between-patient variation, and $\\varepsilon\\sim N(0,%.2f^2)$ assay noise shared across doses. With these values, the mean HbA1c reduction is %s percentage points at %s U/day, respectively (%.2f--%.2f U/kg for a patient of median weight), of the same order as the reductions reported when basal insulin is added to oral therapy \\citep{riddle2003}. The parameters of the dose response are chosen to be clinically plausible but are not fitted to trial data. We use the ADA target HbA1c $\\le7\\%%$ \\citep{ada2025glycemic} as the benefit ($y^u=7.0$) and the intensive target HbA1c $\\le6\\%%$, which increased mortality and severe hypoglycemia in the ACCORD trial \\citep{accord2008}, as the harm ($y^l=6.0$).",
    PAR$Emax, PAR$ED50, WMED, PAR$eta, PAR$sd, paste(sprintf("%.1f", red), collapse = ", "),
    paste(DOSES[-1], collapse = ", "), DOSES[2] / WMED, DOSES[K + 1] / WMED),
  sprintf("Since $Y_k$ is decreasing in $k$ for every patient, Assumptions~\\ref{MON} and~\\ref{MON2} hold. Assumption~\\ref{MIN} is not guaranteed because the latent variable $(A_0,W,\\eta,\\varepsilon)$ is multi-dimensional; because the simulator returns all potential outcomes, we compute the true $\\theta^{\\text{MNH}}$ from the joint counterfactual outcomes and compare it with the identified expression in Theorem~\\ref{theo5}, finding a maximum discrepancy of $%s$.", sciTex(A3err)),
  "",
  sprintf("Table~\\ref{tab:nhanes_truth} shows the ground truth. The Benefit criterion is maximized at %d U/day and the Harmless criterion at %d U/day, whereas the MND-B is %d U/day and the MND-H is %d U/day: %d U/day is exactly necessary for the largest share (%.1f\\%%) of patients. At %d U/day, %s of patients are driven to HbA1c $\\le6\\%%$ and %s are over-treated, compared with %s and %s at %d U/day.",
    DOSES[iB], DOSES[iH], DOSES[kB + 1], DOSES[kH + 1], DOSES[kB + 1], 100 * max(thMNB),
    DOSES[iB], pct(intens[iB]), pct(over[iB]), pct(intens[kB + 1]), pct(over[kB + 1]), DOSES[kB + 1]),
  "",
  "\\noindent{\\bf Results.}",
  sprintf("We set $T^{\\text{initial}}=%d$ and $T=%d$, and repeat each experiment %d times. Table~\\ref{tab:nhanes_choice} shows the chosen doses, Figures~\\ref{fig:nhanesB} and~\\ref{fig:nhanesH} the cumulative regret, and Table~\\ref{tab:nhanes_clinical} the outcomes of the trial patients.", TINIT, TT, REPS),
  sprintf("MND-B UCB identifies the MND-B in %d\\%%, %d\\%% and %d\\%% of runs with $\\kappa=1$, $0.1$ and $0.01$, respectively. Benefit UCB converges to %d U/day and Harmless UCB to %d U/day, and the doses chosen by the MED-based methods depend on $\\phi$. Compared with Benefit UCB, MND-B UCB ($\\kappa=0.01$) administers %.1f instead of %.1f U/day on average and reduces over-intensive control (HbA1c $\\le6\\%%$) from %s to %s and over-treatment from %s to %s, while %s instead of %s of patients reach the target. This reflects the trade-off that the MND makes explicit: it targets the dose that is necessary for the largest share of patients rather than the dose that brings the most patients to target.",
    idr(mB[[1]], kB), idr(mB[[2]], kB), idr(mB[[3]], kB), du(bB[["Benefit UCB"]]), du(bH[["Harmless UCB"]]),
    mb["units"], bu["units"], pct(bu["intens"]), pct(mb["intens"]), pct(bu["over"]), pct(mb["over"]),
    pct(mb["cured"]), pct(bu["cured"])),
  "",
  sprintf("MND-H UCB identifies the MND-H in %d\\%%, %d\\%% and %d\\%% of runs with $\\kappa=1$, $0.1$ and $0.01$, and in %s. The difficulty is the confidence radius of dose~$0$, $r_0=\\sum_{k=1}^K r_k$. When dose~$0$ leads, its confirmation partners are sampled only until their counts reach about $\\log(t)/\\kappa$ (Step~1 of the proof of Theorem~\\ref{theo6}), so each $r_k$ stays near $2\\sqrt{\\kappa}$ and $r_0$ is of order $2K\\sqrt{\\kappa}$; dose~$0$ is therefore abandoned only if roughly $\\kappa\\lesssim(\\Delta_0/2K)^2\\approx%.4f$ in this setting ($\\Delta_0=%.3f$). For MND-B UCB, the radius of dose~$0$ involves only doses $0$ and $K$ and stays near $\\sqrt{2\\kappa}$, giving the much weaker requirement $\\kappa\\lesssim\\Delta_0^2/2\\approx%.2f$. %s This limitation grows with $K$, consistent with the factor $K^2$ in Theorem~\\ref{theo6}; a tighter confidence bound for $\\theta^{\\text{MNH}}(0)$ that exploits the telescoping structure of $\\sum_k\\theta^{\\text{MNH}}(k)$ is an interesting direction for future work.",
    idr(mH[[1]], kH), idr(mH[[2]], kH), idr(mH[[3]], kH), smallTxt, kapStar, D0, kapStarB,
    if (mndHokSmall) "With $\\kappa$ below this threshold, MND-H UCB identifies the MND-H, but a smaller $\\kappa$ also lengthens confirmation (the term $\\log(T)/\\kappa$ in Theorem~\\ref{theo6}), so in this setting MND-H UCB requires $\\kappa$ in a narrow range that depends on the unknown gaps."
    else "Values of $\\kappa$ below this threshold make confirmation last for most of the horizon (the term $\\log(T)/\\kappa$ in Theorem~\\ref{theo6} exceeds the number of rounds available per dose), so within $T=10000$ rounds MND-H UCB does not identify the MND-H in this setting."),
  "",
  "\\begin{figure}[tb]",
  "\\begin{minipage}[t]{.48\\textwidth}\\centering",
  "\\includegraphics[width=\\linewidth]{FigNHANES_bene.pdf}",
  sprintf("\\caption{Mean cumulative regret of MND-B UCB in the semi-synthetic NHANES experiment over %d runs (red solid: $\\kappa=1$; green dashed: $\\kappa=0.1$; blue dotted: $\\kappa=0.01$).}\\label{fig:nhanesB}", REPS),
  "\\end{minipage}\\hfill",
  "\\begin{minipage}[t]{.48\\textwidth}\\centering",
  "\\includegraphics[width=\\linewidth]{FigNHANES_harm.pdf}",
  sprintf("\\caption{Mean cumulative regret of MND-H UCB in the semi-synthetic NHANES experiment over %d runs (red solid: $\\kappa=1$; green dashed: $\\kappa=0.1$; blue dotted: $\\kappa=0.01$).}\\label{fig:nhanesH}", REPS),
  "\\end{minipage}",
  "\\end{figure}",
  "",
  "\\input{table_nhanes_cohort}", "\\input{table_nhanes_truth}", "\\input{table_nhanes_choice}", "\\input{table_nhanes_clinical}"),
  file.path(OUT, "appendix_E3.tex"))
cat("Done. Results in ./", OUT, "\n", sep = "")
