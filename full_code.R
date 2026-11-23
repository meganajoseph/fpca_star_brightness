library(reticulate)
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

# 3. Download in a loop
all_lc <- vector("list", length(targets))
for (i in seq_along(targets)) {
  message(i, "/", length(targets), ": ", targets[i])
  all_lc[[i]] <- tryCatch(get_lc(targets[i]),
                          error = function(e) { message("  failed"); NULL })
}

# 4. Combine into one data frame
df_all <- do.call(rbind, Filter(Negate(is.null), all_lc))
df_all <- subset(df_all, quality == 0)
df_all$KID <- as.integer(sub("KIC ", "", df_all$target))

# 5. Attach catalog properties
df_all <- merge(df_all,
                g_stars[, c("KID", "Teff", "logg", "Mass", "Prot", "Prot_err", "Rper", "LPH")],
                by = "KID")

