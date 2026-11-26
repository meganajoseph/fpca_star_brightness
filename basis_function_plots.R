library(reticulate)
library(ggplot2)
library(refund)
library(tidyverse)
library(dplyr)
library(face)
library(ggpubr)

# ---- Data ----
mcquillan <- read.csv("fpca_star_brightness/data/Table_1_Periodic.txt")
g_stars <- mcquillan %>% filter(Teff >= 5300 & Teff <= 6000 & !is.na(Prot))
set.seed(1)

# Combine all files into one data frame
files <- list.files("fpca_star_brightness/data/lightcurves_fold", 
                    pattern = "\\.csv$", full.names = TRUE)

# Calculate and add phase
df_fold <- map_dfr(files, read.csv) |>
  mutate(KID = as.integer(sub("KIC ", "", target))) |>
  left_join(mcquillan |>
              select(KID, Prot, Prot_err), by = "KID") |>
  mutate(phase = (t0 %% Prot) / Prot)

ggplot(df_fold, aes(phase, flux, group=KID, color=KID)) +
  geom_point(size = 0.1, alpha = 0.3)

ydata <- df_fold |>
  select(c(KID, phase, flux)) |>
  rename(subj=KID, argvals=phase, y=flux)

# ---- Fit various levels of lambda ----
fit_smoothed <- face.sparse(ydata, pve = 0.95)
fit_unsm <- face.sparse(ydata, lambda = exp(-3), pve = 0.95)
f <- face.sparse(ydata, lambda = 0, pve = 0.95)
fk10 <- face.sparse(ydata, knots = 10, lambda = exp(-3), pve = 0.95)
f001 <- face.sparse(ydata, lambda = 0.01, pve = 0.95)
f1 <- face.sparse(ydata, lambda = 1, pve = 0.95)
f100 <- face.sparse(ydata, lambda = 100, pve = 0.95)

fits <- list(
  "lambda = 8103" = fit_smoothed,
  "lambda = 100" = f100,
  "lambda = 1" = f1,
  "lambda = 0.05" = fit_unsm,
  "lambda = 0.05, knots = 10" = fk10,
  "lambda = 0.01" = f001,
  "lambda = 0" = f
)

# Quick summary of each fit
data.frame(fit    = names(fits),
           npc    = sapply(fits, function(f) f$npc),
           sigma2 = sapply(fits, function(f) f$sigma2))

# ---- Basis functions in long format ----
K <- max(sapply(fits, function(f) ncol(f$eigenfunctions)))
get_s <- function(f) if (!is.null(f$argvals.new)) f$argvals.new else f$argvals
ref <- fits[[which.max(sapply(fits, function(f) ncol(f$eigenfunctions)))]]

# long data frame for one method, sign-aligned to reference
make_df <- function(label, s, phi) {
  do.call(rbind, lapply(seq_len(min(K, ncol(phi))), function(k) {
    ref_k <- approx(get_s(ref), ref$eigenfunctions[, k], xout = s)$y
    data.frame(s = s,
               value = phi[, k] * sign(sum(phi[, k] * ref_k, na.rm = TRUE)),
               method = label,
               k = k)
  }))
}

df_phi <- do.call(rbind, lapply(names(fits), function(nm)
  make_df(nm, get_s(fits[[nm]]), fits[[nm]]$eigenfunctions)))

table(df_phi$k, useNA = "ifany")     # should list components 1..K, no NA

lev <- c(names(fits), "PCA")
df_phi$method <- factor(df_phi$method, levels = lev[lev %in% df_phi$method])
df_phi$panel  <- factor(df_phi$k, levels = 1:K,
                        labels = paste0("phi[", 1:K, "](s)"))

lev <- levels(df_phi$method)
pal <- setNames(scales::viridis_pal(end = 0.9)(length(lev)), lev)
if ("PCA" %in% lev) pal["PCA"] <- "red"

# ---- Plot ----
p_basis <- ggplot(df_phi, aes(s, value, color = method, group = method)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, nrow = 2, labeller = label_parsed) +
  scale_color_manual(values = pal, name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.2))) +
  labs(x = "Phase", y = "Basis function") +
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"))
p_basis
