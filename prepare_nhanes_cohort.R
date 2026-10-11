# =====================================================================
# Build nhanes_cohort.csv from the official NHANES 2017-2018 files.
#
# Download these four SAS transport files from the CDC NHANES website
# (https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017)
# into the working directory:
#   DEMO_J.xpt  demographics        RIDAGEYR (age in years)
#   DIQ_J.xpt   diabetes question.  DIQ010 (doctor told you have diabetes: 1 = yes)
#                                   DIQ050 (taking insulin now: 1 = yes)
#   GHB_J.xpt   glycohemoglobin     LBXGH  (HbA1c, %)
#   BMX_J.xpt   body measures       BMXWT  (weight, kg)
#
# Cohort: adults (>= 18 years) with diagnosed diabetes (DIQ010 == 1), not
# currently taking insulin (DIQ050 != 1), with measured HbA1c >= 7% and weight.
#
# Usage: Rscript prepare_nhanes_cohort.R
# Requires the 'foreign' package, which ships with R.
# =====================================================================
library(foreign)
rd <- function(f) read.xport(f)
demo <- rd("DEMO_J.xpt")[, c("SEQN", "RIDAGEYR")]
diq  <- rd("DIQ_J.xpt")[, c("SEQN", "DIQ010", "DIQ050")]
ghb  <- rd("GHB_J.xpt")[, c("SEQN", "LBXGH")]
bmx  <- rd("BMX_J.xpt")[, c("SEQN", "BMXWT")]
d <- Reduce(function(a, b) merge(a, b, by = "SEQN", all.x = TRUE), list(demo, diq, ghb, bmx))

cohort <- subset(d, RIDAGEYR >= 18 & DIQ010 %in% 1 & !(DIQ050 %in% 1) &
                    !is.na(LBXGH) & !is.na(BMXWT) & LBXGH >= 7)
out <- data.frame(SEQN = cohort$SEQN, age = cohort$RIDAGEYR,
                  weight_kg = cohort$BMXWT, hba1c = cohort$LBXGH)
write.csv(out[order(out$SEQN), ], "nhanes_cohort.csv", row.names = FALSE)
cat(sprintf("Wrote nhanes_cohort.csv with %d patients (median HbA1c %.1f, median weight %.1f kg)\n",
            nrow(out), median(out$hba1c), median(out$weight_kg)))
