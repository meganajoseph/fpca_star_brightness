library(reticulate)
library(ggplot2)
library(refund)
library(tidyverse)
library(dplyr)
library(face)
library(fields)
py_install("lightkurve")
lk <- import("lightkurve")
gc_py <- import("gc")

mcquillan <- read.csv("data/Table_1_Periodic.txt")

g_stars <- mcquillan %>% filter(Teff >= 5300 & Teff <= 6000 & !is.na(Prot))

get_lc <- function(target, outdir = "data/lightcurves_fold",
                   window = 90, thin = 10, min_obs = 300) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  file <- file.path(outdir, paste0(gsub(" ", "_", target), ".csv"))
  skip <- paste0(file, ".skip")
  if (file.exists(file)) return(read.csv(file))
  if (file.exists(skip)) return(NULL)   # already known to fail
  
  sr <- lk$search_lightcurve(target, mission = "Kepler", author = "Kepler",
                             exptime = 1800L)
  if (py_len(sr) == 0) { file.create(skip); return(NULL) }
  
  lc <- sr$download_all()$stitch()$remove_nans()
  pdf <- lc$to_pandas()
  df <- data.frame(target = target,
                   time = as.numeric(rownames(pdf)),
                   flux = pdf$flux, flux_err = pdf$flux_err,
                   quality = pdf$quality)
  rm(lc, pdf); gc_py$collect()                    # free the Python memory
  
  df <- df[df$quality == 0, ]
  df <- df[order(df$time), ]
  
  # find the minimum at full resolution, before thinning
  df$flux_smooth <- runmed(df$flux, k = 25)
  t_min <- df$time[which.min(df$flux_smooth)]
  
  df <- df[df$time >= t_min & df$time < t_min + window, ]
  df$t0 <- df$time - t_min
  df <- df[seq(1, nrow(df), by = thin), ]          # every 10th observation
  
  if (nrow(df) < min_obs) { file.create(skip); return(NULL) }
  
  df$rel_flux <- df$flux / median(df$flux) - 1
  write.csv(df, file, row.names = FALSE)
  df
}

set.seed(1)
candidates <- sample(g_stars$KID)    # shuffled list of all G stars
n_target <- 100
results <- list()

# 3. Download in a loop, run once
for (kid in candidates) {
  if (length(results) >= n_target) break
  t <- paste("KIC", kid)
  message(length(results), " collected; trying ", t)
  res <- tryCatch(get_lc(t), error = function(e) { message("  failed"); NULL })
  if (!is.null(res)) results[[t]] <- res
}

# Combine and add phase
df_fold <- bind_rows(results) %>%
  mutate(KID = as.integer(sub("KIC ", "", target))) %>%
  left_join(mcquillan %>% select(KID, Prot, Prot_err), by = "KID") %>%
  mutate(phase = (t0 %% Prot) / Prot)

df_fold %>% count(KID) %>% summarise(min(n), max(n), stars = n())

# 4. Combine into one data frame
files <- list.files("data/lightcurves_fold", pattern = "\\.csv$", full.names = TRUE)

df_fold <- map_dfr(files, read.csv) %>%
  mutate(KID = as.integer(sub("KIC ", "", target))) %>%
  left_join(mcquillan %>% select(KID, Prot, Prot_err), by = "KID") %>%
  mutate(phase = (t0 %% Prot) / Prot)

ggplot(df_fold, aes(phase, flux, group=KID, color=KID)) + geom_point(size = 0.1, alpha = 0.3)

# FPCA
ydata <- df_fold %>% 
  select(c(KID, phase, flux)) %>% 
  rename(subj=KID, argvals=phase, y=flux)
fit_smoothed_001 <- face.sparse(ydata, pve = 0.95, lambda=0.01)
fit_smoothed_1 <- face.sparse(ydata, pve = 0.95, lambda=1)
fit_smoothed_1000 <- face.sparse(ydata, pve = 0.95, lambda=1000)
fit_smoothed_100 <- face.sparse(ydata, pve = 0.95, lambda=100)

# data vis
image(seq(min(ydata$argvals), max(ydata$argvals), length.out=100), seq(min(ydata$argvals), max(ydata$argvals), length.out=100), fit_smoothed$Chat.new)

png("images/smoothed_covariances.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), oma = c(2, 2, 0, 0))
image.plot(fit_smoothed_001$Chat.new, zlim = c(-4e-5, 5e-5), 
           main = bquote(bold("Smoothed (" ~ tau ~ "= 0.01)" )), legend.mar = 10)

image.plot(fit_smoothed_1$Chat.new, zlim = c(-4e-5, 5e-5),
           main = bquote(bold("Smoothed (" ~ tau ~ "= 1)" )), legend.mar = 10)

image.plot(fit_smoothed_100$Chat.new, zlim = c(-4e-5, 5e-5),
           main = bquote(bold("Smoothed (" ~ tau ~ "= 100)" )), legend.mar = 10)

image.plot(fit_smoothed_1000$Chat.new, zlim = c(-4e-5, 5e-5),
           main = bquote(bold("Smoothed (" ~ tau ~ "= 1000)" )), legend.mar = 1)

mtext("Phase of Period", side = 1, line = -1, outer = TRUE)
mtext("Phase of Period", side = 2, line = -1, outer = TRUE)
dev.off()

image.plot(fit_smoothed$Chat.new)
