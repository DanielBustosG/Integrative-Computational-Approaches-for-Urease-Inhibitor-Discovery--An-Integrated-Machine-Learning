# ======================================================
# WT-MetaD (2 CVs) – Full analysis per replica, per ligand, and across ligands
# Outputs organized under each ligand folder and a global AllLigands folder.
# ======================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

theme_set(
  theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2)
    )
)

#------------------ Readers (robust) -------------------

# .cvseq with multiple cv_k captured
read_cvseq_multi <- function(path){
  if(!file.exists(path)) { warning("cvseq not found: ", path); return(tibble()) }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*$|^\\s*#", lines)]
  if(!length(lines)) return(tibble())
  
  parse_line <- function(s){
    s <- gsub(",", ".", s, fixed = TRUE); s <- gsub("\\[|\\]", " ", s)
    toks <- strsplit(s, "\\s+")[[1]]; toks <- toks[nzchar(toks)]
    if(length(toks) < 4) return(NULL)
    i_time <- which(toks == "time")[1]
    i_pot  <- which(toks %in% c("potential","bias"))[1]
    t  <- if(!is.na(i_time) && i_time < length(toks)) suppressWarnings(as.numeric(toks[i_time+1])) else NA_real_
    po <- if(!is.na(i_pot)  && i_pot  < length(toks)) suppressWarnings(as.numeric(toks[i_pot+1]))  else NA_real_
    i_cvs <- grep("^cv_\\d+$", toks)
    vals <- c(time_ps = t, bias_at_cv = po)
    if(length(i_cvs)){
      for(k in seq_along(i_cvs)){
        idx <- i_cvs[k]
        vals <- c(vals, setNames(suppressWarnings(as.numeric(toks[idx+1])), paste0("cv", k)))
      }
    }
    vals
  }
  
  rows <- lapply(lines, parse_line)
  rows <- rows[ sapply(rows, Negate(is.null)) ]
  if(!length(rows)) return(tibble())
  df <- as_tibble(do.call(rbind, rows))
  df <- df |>
    mutate(across(everything(), as.numeric), time_ns = time_ps*0.001) |>
    relocate(time_ps, time_ns)
  df
}

# .kerseq N-D: time_ps, height, (center1,width1, center2,width2, ...)
read_kers_nd <- function(path){
  if(!file.exists(path)) { warning("kerseq not found: ", path); return(tibble()) }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*$|^\\s*#", lines)]
  if(!length(lines)) return(tibble())
  df <- fread(text = lines, header = FALSE)
  ncv <- as.integer((ncol(df)-2)/2)
  if(2 + 2*ncv != ncol(df)) stop("Unexpected kerseq columns: ", ncol(df))
  nm <- c("time_ps","height")
  for(i in 1:ncv){ nm <- c(nm, paste0("center",i), paste0("width",i)) }
  setnames(df, nm)
  as_tibble(df) |>
    mutate(across(everything(), as.numeric), time_ns = time_ps*0.001) |>
    relocate(time_ps, time_ns)
}

# .fes (1D or 2D)
read_fes_flex <- function(path){
  if(!file.exists(path)) { warning("fes not found: ", path); return(tibble()) }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*$|^\\s*#", lines)]
  if(!length(lines)) return(tibble())
  df <- fread(text = lines, header = FALSE) |> as_tibble()
  if(ncol(df) == 2){
    setNames(df, c("cv","F"))
  } else if(ncol(df) >= 3){
    setNames(df[,1:3], c("cv1","cv2","F"))
  } else tibble()
}

#------------------ FES helpers ------------------------

offset_global_min <- function(df){ df |> mutate(F_off = F - min(F, na.rm = TRUE)) }

offset_window_min <- function(df, window = 0.2){
  if("cv" %in% names(df)){
    x0 <- df$cv[which.min(df$F)]
    near <- df$cv >= x0-window & df$cv <= x0+window
  } else {
    i0 <- which.min(df$F); x0 <- df$cv1[i0]; y0 <- df$cv2[i0]
    near <- (df$cv1 - x0)^2 + (df$cv2 - y0)^2 <= window^2
  }
  df |> mutate(F_off = F - min(F[near], na.rm = TRUE))
}

offset_ref_cv <- function(df, ref){
  if("cv" %in% names(df)){
    i <- which.min(abs(df$cv - ref))
  } else {
    i <- which.min((df$cv1 - ref[1])^2 + (df$cv2 - ref[2])^2)
  }
  df |> mutate(F_off = F - F[i])
}

# Regular-grid checks & bilinear resampling (no external deps)
is_regular_grid <- function(df){
  all(c("cv1","cv2","F_off") %in% names(df)) &&
    (length(unique(df$cv1)) * length(unique(df$cv2)) == nrow(df))
}

reshape_to_matrix <- function(df){
  X <- sort(unique(df$cv1)); Y <- sort(unique(df$cv2))
  Z <- matrix(NA_real_, nrow=length(X), ncol=length(Y))
  ix <- match(df$cv1, X); iy <- match(df$cv2, Y)
  Z[cbind(ix, iy)] <- df$F_off
  list(X=X, Y=Y, Z=Z)
}

bilinear_point <- function(x, y, X, Y, Z){
  if(x < min(X) || x > max(X) || y < min(Y) || y > max(Y)) return(NA_real_)
  i <- findInterval(x, X, all.inside = TRUE)
  j <- findInterval(y, Y, all.inside = TRUE)
  i2 <- min(i+1, length(X)); j2 <- min(j+1, length(Y))
  x1 <- X[i]; x2 <- X[i2]; y1 <- Y[j]; y2 <- Y[j2]
  Q11 <- Z[i , j ]; Q21 <- Z[i2, j ]; Q12 <- Z[i , j2]; Q22 <- Z[i2, j2]
  if(anyNA(c(Q11,Q21,Q12,Q22))){
    di <- which.min(c((x-x1)^2+(y-y1)^2,(x-x2)^2+(y-y1)^2,(x-x1)^2+(y-y2)^2,(x-x2)^2+(y-y2)^2))
    return(c(Q11,Q21,Q12,Q22)[di])
  }
  tx <- if(x2==x1) 0 else (x - x1)/(x2 - x1)
  ty <- if(y2==y1) 0 else (y - y1)/(y2 - y1)
  (1-tx)*(1-ty)*Q11 + tx*(1-ty)*Q21 + (1-tx)*ty*Q12 + tx*ty*Q22
}

resample_to_grid <- function(df, grid){
  stopifnot(is_regular_grid(df))
  M <- reshape_to_matrix(df)
  XX <- rep(grid$x, times=length(grid$y))
  YY <- rep(grid$y, each = length(grid$x))
  vals <- mapply(bilinear_point, XX, YY, MoreArgs = list(X=M$X, Y=M$Y, Z=M$Z))
  tibble(cv1 = XX, cv2 = YY, F_off = vals)
}

# Common grid as intersection of replica ranges (avoids extrapolation)
common_grid_2d_intersection <- function(fes_list, step_x = NULL, step_y = NULL){
  rx <- sapply(fes_list, \(d) range(d$cv1, na.rm=TRUE))
  ry <- sapply(fes_list, \(d) range(d$cv2, na.rm=TRUE))
  xr <- c(max(rx[1,]), min(rx[2,])); yr <- c(max(ry[1,]), min(ry[2,]))
  if(xr[1] >= xr[2] || yr[1] >= yr[2]) stop("No common overlap across replicas.")
  if(is.null(step_x)) step_x <- median(sapply(fes_list, \(d) median(diff(sort(unique(d$cv1))), na.rm=TRUE)), na.rm=TRUE)
  if(is.null(step_y)) step_y <- median(sapply(fes_list, \(d) median(diff(sort(unique(d$cv2))), na.rm=TRUE)), na.rm=TRUE)
  list(x = seq(xr[1], xr[2], by = step_x), y = seq(yr[1], yr[2], by = step_y))
}

#------------------ Plots -------------------------------

plot_cv_timeseries_multi <- function(cvseq){
  if(!nrow(cvseq)) return(ggplot() + labs(title="CV time series — no data"))
  df <- cvseq |> select(time_ns, starts_with("cv")) |>
    pivot_longer(-time_ns, names_to = "cv", values_to = "val")
  ggplot(df, aes(time_ns, val, color=cv)) +
    geom_line(linewidth=0.7) +
    labs(x="Time (ns)", y="CV (Å)", title="CV time series (all CVs)")
}

plot_cv_hists_multi <- function(cvseq, bins=40){
  if(!nrow(cvseq)) return(ggplot() + labs(title="CV histogram — no data"))
  df <- cvseq |> select(starts_with("cv")) |>
    pivot_longer(everything(), names_to="cv", values_to="val")
  ggplot(df, aes(val, fill=cv)) +
    geom_histogram(aes(y=after_stat(density)), bins=bins, alpha=0.4, position="identity") +
    labs(x="CV (Å)", y="Probability density", title="CV histograms")
}

plot_fes2d <- function(fes2, title="FES 2D"){
  if(!nrow(fes2) || !all(c("cv1","cv2","F_off") %in% names(fes2)))
    return(ggplot() + labs(title=paste(title, "— no data")))
  ggplot(fes2, aes(cv1, cv2, fill = F_off)) +
    geom_raster(interpolate=TRUE) +
    stat_contour(aes(z=F_off), color="white", linewidth=0.35) +
    scale_fill_viridis_c() +
    labs(x="CV1 (Å)", y="CV2 (Å)", fill="Free energy (kcal/mol)", title=title)
}

plot_fes2d_mean <- function(df, title="FES 2D (mean ± SEM)"){
  if(!nrow(df)) return(ggplot() + labs(title=paste(title, "— no data")))
  ggplot(df, aes(cv1, cv2, fill = F_mean)) +
    geom_raster(interpolate = TRUE) +
    stat_contour(aes(z=F_mean), color="white", linewidth=0.35) +
    scale_fill_viridis_c() +
    labs(x="CV1 (Å)", y="CV2 (Å)", fill="F_mean (kcal/mol)", title=title)
}

plot_projection_1d <- function(df, title, xlab){
  if(!nrow(df)) return(ggplot() + labs(title=paste(title,"— no data")))
  ggplot(df, aes(cv, F_mean)) + geom_line(linewidth=1) +
    labs(x = xlab, y = "Free energy (offset to 0)", title = title)
}

#------------------ Reconstructions --------------------

# 2D reconstruction from kernels up to t_ns
reconstruct_bias_from_kers_2d <- function(kers, grid, t_ns){
  sel <- kers |> dplyr::filter(time_ns <= t_ns)
  if(nrow(sel)==0) return(NULL)
  X <- grid$x; Y <- grid$y
  V <- matrix(0, nrow=length(X), ncol=length(Y))
  for(r in seq_len(nrow(sel))){
    v1 <- exp(-0.5 * ((X - sel$center1[r]) / sel$width1[r])^2)
    v2 <- exp(-0.5 * ((Y - sel$center2[r]) / sel$width2[r])^2)
    V <- V + sel$height[r] * (outer(v1, v2))
  }
  F <- -V
  expand_grid(cv1 = X, cv2 = Y) |> mutate(F = as.vector(F))
}

#------------------ Projection (Boltzmann) --------------

# F_proj(x) = -kT * log( mean_y exp(-beta F(x,y)) )
project_1d_from_2d <- function(fes2, which = c("cv1","cv2"), T = 300){
  which <- match.arg(which)
  if(!nrow(fes2)) return(tibble())
  E <- if("F_off" %in% names(fes2)) "F_off" else if("F_mean" %in% names(fes2)) "F_mean" else "F"
  kT <- 0.0019872041 * T
  if(which=="cv1"){
    fes2 %>% filter(is.finite(.data[[E]])) %>%
      group_by(cv1) %>%
      summarise(F_mean = -kT * log(mean(exp(-.data[[E]]/kT))), .groups="drop") %>%
      rename(cv = cv1)
  } else {
    fes2 %>% filter(is.finite(.data[[E]])) %>%
      group_by(cv2) %>%
      summarise(F_mean = -kT * log(mean(exp(-.data[[E]]/kT))), .groups="drop") %>%
      rename(cv = cv2)
  }
}

#==================== Aggregations across replicas (CV hist & time series)

common_bins_from_reps <- function(cvseq_list, cv_col = "cv1", nbins = 60){
  vals <- unlist(lapply(cvseq_list, \(d) d[[cv_col]]), use.names = FALSE)
  rng  <- range(vals, na.rm = TRUE)
  seq(rng[1], rng[2], length.out = nbins + 1)
}

mean_sem_hist <- function(cvseq_list, cv_col = "cv1", nbins = 60){
  breaks <- common_bins_from_reps(cvseq_list, cv_col, nbins)
  mids   <- 0.5*(breaks[-1] + breaks[-length(breaks)])
  dens_mat <- sapply(cvseq_list, function(df){
    if(!nrow(df) || !cv_col %in% names(df)) return(rep(NA_real_, length(mids)))
    h <- hist(df[[cv_col]], breaks = breaks, plot = FALSE)
    h$density
  })
  m <- rowMeans(dens_mat, na.rm = TRUE)
  s <- apply(dens_mat, 1, sd, na.rm = TRUE)
  n_eff <- rowSums(is.finite(dens_mat))
  sem <- s / sqrt(pmax(1, n_eff))
  tibble(mid = mids, dens_mean = m, dens_sem = sem)
}

common_time_grid <- function(cvseq_list, dt = 0.1){
  tmins <- sapply(cvseq_list, \(d) if(nrow(d)) min(d$time_ns, na.rm=TRUE) else Inf)
  tmaxs <- sapply(cvseq_list, \(d) if(nrow(d)) max(d$time_ns, na.rm=TRUE) else -Inf)
  t0 <- max(tmins[is.finite(tmins)])
  t1 <- min(tmaxs[is.finite(tmaxs)])
  if(!is.finite(t0) || !is.finite(t1) || t1 <= t0) return(numeric())
  seq(t0, t1, by = dt)
}

interp_rep_to_timegrid <- function(df, cv_col, grid_t){
  if(!nrow(df) || !cv_col %in% names(df) || !length(grid_t)) return(rep(NA_real_, length(grid_t)))
  a <- approx(x = df$time_ns, y = df[[cv_col]], xout = grid_t, rule = 2, ties = "ordered")
  a$y
}

mean_sem_timeseries <- function(cvseq_list, cv_col = "cv1", dt = 0.1){
  grid_t <- common_time_grid(cvseq_list, dt = dt)
  if(!length(grid_t)) return(tibble())
  mat <- sapply(cvseq_list, interp_rep_to_timegrid, cv_col = cv_col, grid_t = grid_t)
  m <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  n_eff <- rowSums(is.finite(mat))
  sem <- s / sqrt(pmax(1, n_eff))
  tibble(time_ns = grid_t, mean = m, sem = sem)
}

plot_cv_hist_avg_two <- function(h1, h2, title = "CV histograms (mean ± SEM)"){
  ggplot() +
    geom_ribbon(data=h1, aes(x=mid, ymin=dens_mean-dens_sem, ymax=dens_mean+dens_sem, fill="cv1"), alpha=0.20) +
    geom_line  (data=h1, aes(x=mid, y=dens_mean, color="cv1"), linewidth=1.0) +
    geom_ribbon(data=h2, aes(x=mid, ymin=dens_mean-dens_sem, ymax=dens_mean+dens_sem, fill="cv2"), alpha=0.20) +
    geom_line  (data=h2, aes(x=mid, y=dens_mean, color="cv2"), linewidth=1.0) +
    scale_color_manual(values=c("cv1"="#F8766D","cv2"="#00BFC4")) +
    scale_fill_manual(values =c("cv1"="#F8766D","cv2"="#00BFC4")) +
    labs(x="CV (Å)", y="Probability density", title=title, color=NULL, fill=NULL)
}

plot_cv_timeseries_avg_two <- function(ts1, ts2, title = "CV time series (mean ± SEM)"){
  ggplot() +
    geom_ribbon(data=ts1, aes(x=time_ns, ymin=mean-sem, ymax=mean+sem, fill="cv1"), alpha=0.20) +
    geom_line  (data=ts1, aes(x=time_ns, y=mean, color="cv1"), linewidth=0.9) +
    geom_ribbon(data=ts2, aes(x=time_ns, ymin=mean-sem, ymax=mean+sem, fill="cv2"), alpha=0.20) +
    geom_line  (data=ts2, aes(x=time_ns, y=mean, color="cv2"), linewidth=0.9) +
    scale_color_manual(values=c("cv1"="#F8766D","cv2"="#00BFC4")) +
    scale_fill_manual(values =c("cv1"="#F8766D","cv2"="#00BFC4")) +
    labs(x="Time (ns)", y="CV (Å)", title=title, color=NULL, fill=NULL)
}

#------------------ Per-replica analysis ----------------

analyze_replica_2cv <- function(prefix_without_ext,
                                times_ns = c(25,50,75,100),
                                out_dir = dirname(prefix_without_ext),
                                offset_method = c("global_min","window_min","ref_cv"),
                                offset_window = 0.2,
                                offset_ref = c(NA,NA)){
  offset_method <- match.arg(offset_method)
  
  cvseq <- read_cvseq_multi(paste0(prefix_without_ext, ".cvseq"))
  kers  <- read_kers_nd   (paste0(prefix_without_ext, ".kerseq"))
  fes   <- read_fes_flex  (paste0(prefix_without_ext, ".fes"))
  
  message("Rows: cvseq=", nrow(cvseq), "  kers=", nrow(kers), "  fes=", nrow(fes))
  
  p_ts   <- plot_cv_timeseries_multi(cvseq)
  p_hsts <- plot_cv_hists_multi(cvseq)
  
  # FES 2D: use engine FES if present; else reconstruct from kernels
  fes2 <- tibble()
  if(all(c("cv1","cv2","F") %in% names(fes))){
    fes2 <- offset_global_min(fes)
  } else if(nrow(kers) && all(c("center1","center2") %in% names(kers))){
    rx <- range(cvseq$cv1, na.rm = TRUE); ry <- range(cvseq$cv2, na.rm = TRUE)
    grid <- list(x = seq(rx[1], rx[2], length.out = 60),
                 y = seq(ry[1], ry[2], length.out = 60))
    recon <- reconstruct_bias_from_kers_2d(kers, grid, max(kers$time_ns, na.rm = TRUE))
    if(!is.null(recon)) fes2 <- offset_global_min(recon)
  }
  
  p_fes2 <- if(nrow(fes2)) plot_fes2d(fes2, "FES 2D (offset to 0)") else ggplot() + labs(title="FES 2D — no data")
  
  # Convergence maps (if kernels exist)
  conv2_plots <- list()
  if(nrow(kers)){
    rx <- range(cvseq$cv1, na.rm = TRUE); ry <- range(cvseq$cv2, na.rm = TRUE)
    grid <- list(x = seq(rx[1], rx[2], length.out = 60),
                 y = seq(ry[1], ry[2], length.out = 60))
    tmax <- max(kers$time_ns, na.rm = TRUE)
    for(tn in times_ns[times_ns <= tmax]){
      rc <- reconstruct_bias_from_kers_2d(kers, grid, tn)
      if(!is.null(rc)){
        rc_off <- offset_global_min(rc)
        conv2_plots[[paste0(tn," ns")]] <- plot_fes2d(rc_off, paste0("Convergence 2D @ ", tn, " ns"))
      }
    }
  }
  
  figs_dir <- file.path(out_dir, "figs")
  dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
  bname <- basename(prefix_without_ext)
  ggsave(file.path(figs_dir, paste0(bname, "_cv_timeseries_allCVs.png")), p_ts,   width=8, height=5, dpi=200)
  ggsave(file.path(figs_dir, paste0(bname, "_cv_hist_allCVs.png")),     p_hsts,  width=7, height=5, dpi=200)
  ggsave(file.path(figs_dir, paste0(bname, "_FES2D.png")),              p_fes2,  width=7.2, height=6.2, dpi=200)
  if(length(conv2_plots)){
    for(nm in names(conv2_plots)){
      ggsave(file.path(figs_dir, paste0(bname, "_FES2D_conv_", gsub(" ","",nm), ".png")),
             conv2_plots[[nm]], width=7.2, height=6.2, dpi=200)
    }
  }
  if(nrow(fes2)) readr::write_csv(fes2, file.path(figs_dir, paste0(bname, "_FES2D_offset.csv")))
  
  list(cvseq=cvseq, kers=kers, fes2=fes2)
}

#------------------ Combine replicas (mean ± SEM) -------

combine_fes2d_mean_sem <- function(fes2_list, min_reps = 2){
  fes2_list <- fes2_list[sapply(fes2_list, nrow) > 0]
  if(!length(fes2_list)) return(tibble())
  fes2_list <- lapply(fes2_list, \(d) if("F_off" %in% names(d)) d else d %>% mutate(F_off = F - min(F,na.rm=TRUE)))
  if(!all(sapply(fes2_list, is_regular_grid)))
    stop("At least one replica FES is not on a regular grid; cannot combine safely.")
  grid <- common_grid_2d_intersection(fes2_list)
  mats <- lapply(fes2_list, resample_to_grid, grid = grid)
  df <- reduce(mats, full_join, by=c("cv1","cv2"), suffix = c(".r1",".r2"))
  energy_cols <- setdiff(names(df), c("cv1","cv2"))
  long <- df %>% pivot_longer(all_of(energy_cols), names_to="rep", values_to="val")
  agg <- long %>% group_by(cv1,cv2) %>%
    summarise(n = sum(is.finite(val)),
              F_mean = ifelse(n >= min_reps, mean(val, na.rm=TRUE), NA_real_),
              F_sem  = ifelse(n >= min_reps, sd(val, na.rm=TRUE)/sqrt(n), NA_real_),
              .groups="drop")
  agg %>% filter(is.finite(F_mean))
}

projections_from_mean <- function(fes2_mean){
  list(
    proj_cv1 = project_1d_from_2d(fes2_mean, "cv1"),
    proj_cv2 = project_1d_from_2d(fes2_mean, "cv2")
  )
}

#------------------ Per-ligand pipeline -----------------

run_ligand <- function(
    ligand_tag,
    prefix_base = NULL,           # if NULL → "WtMetD_100ns_<LIG>_2CV_rep"
    reps = 1:3,
    base_dir = file.path(getwd(), ligand_tag),
    times_ns_fixed = c(25,50,75,100)
){
  if(is.null(prefix_base)) prefix_base <- sprintf("WtMetD_100ns_%s_2CV_rep", ligand_tag)
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  
  results <- vector("list", length(reps)); names(results) <- paste0("rep", reps)
  for(r in reps){
    prefix <- file.path(base_dir, sprintf("%s%d", prefix_base, r))
    outdir <- file.path(base_dir, sprintf("rep%d", r))
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    ktmp <- read_kers_nd(paste0(prefix, ".kerseq"))
    tmax <- if(nrow(ktmp)) max(ktmp$time_ns, na.rm = TRUE) else 0
    tchecks <- times_ns_fixed[times_ns_fixed <= tmax]
    if(!length(tchecks) && is.finite(tmax) && tmax>0) tchecks <- round(seq(tmax/4, tmax, length.out=4),1)
    
    results[[paste0("rep",r)]] <- analyze_replica_2cv(prefix, times_ns = tchecks, out_dir = outdir)
  }
  
  # Combine replicas for ligand
  alldir  <- file.path(base_dir, "All"); dir.create(alldir, showWarnings = FALSE, recursive = TRUE)
  figsAll <- file.path(alldir, "figs");  dir.create(figsAll, showWarnings = FALSE, recursive = TRUE)
  
  fes2_list <- lapply(results, `[[`, "fes2")
  fes2_mean <- combine_fes2d_mean_sem(fes2_list, min_reps = 2)
  p_mean    <- plot_fes2d_mean(fes2_mean, sprintf("FES 2D (mean ± SEM) — %s", ligand_tag))
  ggsave(file.path(figsAll, sprintf("%s_All_FES2D_mean.png", ligand_tag)), p_mean, width=7.5, height=6.5, dpi=200)
  if(nrow(fes2_mean)) readr::write_csv(fes2_mean, file.path(figsAll, sprintf("%s_All_FES2D_mean.csv", ligand_tag)))
  
  projs <- projections_from_mean(fes2_mean)
  if(nrow(projs$proj_cv1)){
    ggsave(file.path(figsAll, sprintf("%s_All_proj_CV1.png", ligand_tag)),
           plot_projection_1d(projs$proj_cv1, sprintf("Projection on CV1 — %s", ligand_tag), "CV1 (Å)"),
           width=7.5, height=4.8, dpi=200)
    readr::write_csv(projs$proj_cv1, file.path(figsAll, sprintf("%s_All_proj_CV1.csv", ligand_tag)))
  }
  if(nrow(projs$proj_cv2)){
    ggsave(file.path(figsAll, sprintf("%s_All_proj_CV2.png", ligand_tag)),
           plot_projection_1d(projs$proj_cv2, sprintf("Projection on CV2 — %s", ligand_tag), "CV2 (Å)"),
           width=7.5, height=4.8, dpi=200)
    readr::write_csv(projs$proj_cv2, file.path(figsAll, sprintf("%s_All_proj_CV2.csv", ligand_tag)))
  }
  
  # ===== Promedio por réplica (CV hist & time series) =====
  cvseq_list <- lapply(results, `[[`, "cvseq")
  
  # Hist mean±SEM (CV1 y CV2)
  h_cv1 <- mean_sem_hist(cvseq_list, cv_col = "cv1", nbins = 60)
  h_cv2 <- mean_sem_hist(cvseq_list, cv_col = "cv2", nbins = 60)
  if(nrow(h_cv1) && nrow(h_cv2)){
    p_hist_avg <- plot_cv_hist_avg_two(h_cv1, h_cv2,
                                       sprintf("CV histograms (mean ± SEM) — %s", ligand_tag))
    ggsave(file.path(figsAll, sprintf("%s_All_CV_hist_avg.png", ligand_tag)),
           p_hist_avg, width=7.5, height=5.2, dpi=200)
    readr::write_csv(h_cv1, file.path(figsAll, sprintf("%s_All_CV1_hist_mean_sem.csv", ligand_tag)))
    readr::write_csv(h_cv2, file.path(figsAll, sprintf("%s_All_CV2_hist_mean_sem.csv", ligand_tag)))
  }
  
  # Time series mean±SEM (CV1 y CV2) – grilla común de tiempo (dt=0.1 ns)
  ts_cv1 <- mean_sem_timeseries(cvseq_list, cv_col = "cv1", dt = 0.1)
  ts_cv2 <- mean_sem_timeseries(cvseq_list, cv_col = "cv2", dt = 0.1)
  if(nrow(ts_cv1) && nrow(ts_cv2)){
    p_ts_avg <- plot_cv_timeseries_avg_two(ts_cv1, ts_cv2,
                                           sprintf("CV time series (mean ± SEM) — %s", ligand_tag))
    ggsave(file.path(figsAll, sprintf("%s_All_CV_timeseries_avg.png", ligand_tag)),
           p_ts_avg, width=8.0, height=5.2, dpi=200)
    readr::write_csv(ts_cv1, file.path(figsAll, sprintf("%s_All_CV1_timeseries_mean_sem.csv", ligand_tag)))
    readr::write_csv(ts_cv2, file.path(figsAll, sprintf("%s_All_CV2_timeseries_mean_sem.csv", ligand_tag)))
  }
  
  invisible(list(results=results, fes2_mean=fes2_mean, projs=projs))
}

#------------------ Cross-ligand comparisons ------------

compare_across_ligands <- function(ligands = c("AHA","DJM","BME","CA1","CA3","CA6"),
                                   base_parent = getwd()){
  outdir <- file.path(base_parent, "AllLigands"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  figs   <- file.path(outdir, "figs");         dir.create(figs,     showWarnings = FALSE, recursive = TRUE)
  
  fes_list <- list()
  for(L in ligands){
    fcsv <- file.path(base_parent, L, "All", "figs", sprintf("%s_All_FES2D_mean.csv", L))
    if(file.exists(fcsv)){
      fes_list[[L]] <- readr::read_csv(fcsv, show_col_types = FALSE)
    }
  }
  if(!length(fes_list)) { warning("No per-ligand mean FES found to compare."); return(invisible()) }
  
  # 1D Boltzmann projections for all ligands
  proj_cv1 <- bind_rows(lapply(names(fes_list), function(L){
    x <- project_1d_from_2d(fes_list[[L]] %>% rename(F_mean = F_mean), which = "cv1")
    mutate(x, ligand = L)
  }))
  proj_cv2 <- bind_rows(lapply(names(fes_list), function(L){
    x <- project_1d_from_2d(fes_list[[L]] %>% rename(F_mean = F_mean), which = "cv2")
    mutate(x, ligand = L)
  }))
  
  if(nrow(proj_cv1)){
    p1 <- ggplot(proj_cv1, aes(cv, F_mean, color=ligand)) + geom_line(linewidth=1) +
      labs(x="CV1 (Å)", y="Free energy (offset to 0)", title="Projection on CV1 — all ligands")
    ggsave(file.path(figs, "AllLigands_proj_CV1.png"), p1, width=8, height=5, dpi=200)
    readr::write_csv(proj_cv1, file.path(figs, "AllLigands_proj_CV1.csv"))
  }
  if(nrow(proj_cv2)){
    p2 <- ggplot(proj_cv2, aes(cv, F_mean, color=ligand)) + geom_line(linewidth=1) +
      labs(x="CV2 (Å)", y="Free energy (offset to 0)", title="Projection on CV2 — all ligands")
    ggsave(file.path(figs, "AllLigands_proj_CV2.png"), p2, width=8, height=5, dpi=200)
    readr::write_csv(proj_cv2, file.path(figs, "AllLigands_proj_CV2.csv"))
  }
  
  # Mean FES panel (facets)
  fes_bind <- bind_rows(lapply(names(fes_list), function(L) mutate(fes_list[[L]], ligand=L)))
  if(nrow(fes_bind)){
    pgrid <- ggplot(fes_bind, aes(cv1, cv2, fill = F_mean)) +
      geom_raster(interpolate = TRUE) +
      stat_contour(aes(z=F_mean), color="white", linewidth=0.25) +
      scale_fill_viridis_c() + facet_wrap(~ligand, ncol = 3, scales = "free") +
      labs(x="CV1 (Å)", y="CV2 (Å)", fill="F_mean (kcal/mol)", title="FES 2D (mean) — all ligands")
    ggsave(file.path(figs, "AllLigands_FES2D_mean_grid.png"), pgrid, width=12, height=8, dpi=180)
  }
  
  # ===== Cross-ligand CV histograms (replicate-averaged) =====
  hist_cv1_all <- bind_rows(lapply(ligands, function(L){
    f <- file.path(base_parent, L, "All", "figs", sprintf("%s_All_CV1_hist_mean_sem.csv", L))
    if(file.exists(f)) readr::read_csv(f, show_col_types = FALSE) %>% mutate(ligand = L) else NULL
  })) %>% bind_rows()
  
  hist_cv2_all <- bind_rows(lapply(ligands, function(L){
    f <- file.path(base_parent, L, "All", "figs", sprintf("%s_All_CV2_hist_mean_sem.csv", L))
    if(file.exists(f)) readr::read_csv(f, show_col_types = FALSE) %>% mutate(ligand = L) else NULL
  })) %>% bind_rows()
  
  if(nrow(hist_cv1_all)){
    p_h1 <- ggplot(hist_cv1_all, aes(mid, dens_mean, color=ligand)) +
      geom_line(linewidth=1) +
      labs(x="CV1 (Å)", y="Probability density", title="CV1 histogram — all ligands (replicate-averaged)")
    ggsave(file.path(figs, "AllLigands_CV1_hist_avg.png"), p_h1, width=8, height=5, dpi=200)
  }
  if(nrow(hist_cv2_all)){
    p_h2 <- ggplot(hist_cv2_all, aes(mid, dens_mean, color=ligand)) +
      geom_line(linewidth=1) +
      labs(x="CV2 (Å)", y="Probability density", title="CV2 histogram — all ligands (replicate-averaged)")
    ggsave(file.path(figs, "AllLigands_CV2_hist_avg.png"), p_h2, width=8, height=5, dpi=200)
  }
  
  # ===== Cross-ligand CV time series (replicate-averaged) =====
  ts_cv1_all <- bind_rows(lapply(ligands, function(L){
    f <- file.path(base_parent, L, "All", "figs", sprintf("%s_All_CV1_timeseries_mean_sem.csv", L))
    if(file.exists(f)) readr::read_csv(f, show_col_types = FALSE) %>% mutate(ligand = L) else NULL
  })) %>% bind_rows()
  
  ts_cv2_all <- bind_rows(lapply(ligands, function(L){
    f <- file.path(base_parent, L, "All", "figs", sprintf("%s_All_CV2_timeseries_mean_sem.csv", L))
    if(file.exists(f)) readr::read_csv(f, show_col_types = FALSE) %>% mutate(ligand = L) else NULL
  })) %>% bind_rows()
  
  if(nrow(ts_cv1_all)){
    p_t1 <- ggplot(ts_cv1_all, aes(time_ns, mean, color=ligand)) +
      geom_line(linewidth=0.9) +
      labs(x="Time (ns)", y="CV1 (Å)", title="CV1 time series — all ligands (replicate-averaged)")
    ggsave(file.path(figs, "AllLigands_CV1_timeseries_avg.png"), p_t1, width=8.5, height=5, dpi=200)
  }
  if(nrow(ts_cv2_all)){
    p_t2 <- ggplot(ts_cv2_all, aes(time_ns, mean, color=ligand)) +
      geom_line(linewidth=0.9) +
      labs(x="Time (ns)", y="CV2 (Å)", title="CV2 time series — all ligands (replicate-averaged)")
    ggsave(file.path(figs, "AllLigands_CV2_timeseries_avg.png"), p_t2, width=8.5, height=5, dpi=200)
  }
  
  invisible(list(fes=fes_list, proj_cv1=proj_cv1, proj_cv2=proj_cv2))
}

#==================== EXECUTION =========================
# Set your working directory to the folder containing AHA/, DJM/, BME/, CA1/, CA3/, CA6/
# Example:
setwd("~/Library/CloudStorage/Dropbox/Urease_2024/Pharm_ML_VirtualScreening_Triazoles/Urease_Triazol/09_MetaDynamics/Outputs")
ligs <- c("AHA","DJM","BME","CA1","CA3","CA6")
for(L in ligs) run_ligand(L, reps = 1:3, base_dir = file.path(getwd(), L))
compare_across_ligands(ligs, base_parent = getwd())