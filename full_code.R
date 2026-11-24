library(reticulate)
library(ggplot2)
library(refund)
library(tidyverse)
library(dplyr)
library(face)
py_install("lightkurve")
lk <- import("lightkurve")

mcquillan <- read.csv("data/Table_1_Periodic.txt")

g_stars <- subset(mcquillan, Teff >= 5200 & Teff <= 6000)
set.seed(1)
sample_ids <- sample(g_stars$KID, 100)
targets <- paste("KIC", sample_ids)

get_lc <- function(target, outdir = "data/lightcurves", quarter = NULL) {
  dir.create(outdir, showWarnings = FALSE)
  file <- file.path(outdir, paste0(gsub(" ", "_", target), ".csv"))
  if (file.exists(file)) return(read.csv(file))
  
  sr <- lk$search_lightcurve(target, mission = "Kepler", author = "Kepler",
                             exptime = 1800L, quarter = quarter)
  if (py_len(sr) == 0) return(NULL)
  
  lc <- sr$download_all()$stitch()$remove_nans()
  df <- lc$to_pandas()
  df$time <- as.numeric(rownames(df))
  rownames(df) <- NULL
  df <- data.frame(target = target, df[, c("time", "flux", "flux_err", "quality")])
  
  write.csv(df, file, row.names = FALSE)
  df
}

# 3. Download in a loop, run once
all_lc <- vector("list", length(targets))
for (i in seq_along(targets)) {
  message(i, "/", length(targets), ": ", targets[i])
  all_lc[[i]] <- tryCatch(get_lc(targets[i]),
                          error = function(e) { message("  failed"); NULL })
}

# 4. Combine into one data frame
files <- list.files("data/lightcurves", pattern = "\\.csv$", full.names = TRUE)
df_all <- do.call(rbind, lapply(files, read.csv))
df_all <- subset(df_all, quality == 0)
df_all$KID <- as.integer(sub("KIC ", "", df_all$target))
df_all <- merge(df_all,
                mcquillan[, c("KID", "Prot", "Prot_err")],
                by = "KID")

df_fold <- df_all |>
  group_by(KID) |>
  arrange(time, .by_group = TRUE) |>
  mutate(
    flux_smooth = runmed(flux, k = 25),       # ~12 hours, reduces noise
    t_min = time[which.min(flux_smooth)]          # time of the minimum
  ) |>
  filter(time >= t_min, time < t_min + 90) |>     # 90-day window
  mutate(
    t0    = time - t_min,                         # minimum now at t0 = 0
    phase = (t0 %% Prot) / Prot
  ) |>
  ungroup()

ggplot(df_fold, aes(phase, flux, group=KID, color=KID)) + geom_point(size = 0.1, alpha = 0.3)

ydata <- df_fold %>% 
  select(c(KID, phase, flux)) %>% 
  rename(subj=KID, argvals=phase, y=flux)
fit_smoothed <- face.sparse(ydata, pve = 0.95)
