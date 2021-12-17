#' Chi - Squared test for hydrological frequency analysis
#' @export
#' @param df Data frame with observations and distributions adjusted
#' @param nc Number of classes
#' @param dist_param Number of parameters of the distributions in df
#' @param alpha Significance level
#' @import grDevices ggplot2 reshape2 stats ggpubr

hydroCHISQR <- function(df, nc, dist_param, alpha) {
  # Number of classes for histograms
  # nc = 1 by default
  if (nc == 3) {
    nbreaks <- grDevices::nclass.FD(df$Obs)
  } else if (nc == 2) {
    nbreaks <- grDevices::nclass.scott(df$Obs)
  } else {
    nbreaks <- grDevices::nclass.Sturges(df$Obs)
  }

  minval <- min(df$Obs)
  maxval <- max(df$Obs)
  int <- round((maxval - minval) / nbreaks, 3) # round to 1 digit more than obs or sim
  brks <- c(minval, minval + seq(nbreaks-1) * int, maxval)

  # Intervals
  factor_obs <- cut(df$Obs, breaks=brks, include.lowest=TRUE)
  obs_out <- as.data.frame(table(factor_obs))
  obs_out <- transform(obs_out, cumFreq = cumsum(obs_out$Freq), relative = prop.table(obs_out$Freq))

  sim_factor <- list()
  for(i in 1:(dim(df)[2]-1)) {
    factor_sim <- cut(df[,i+1], breaks=brks, include.lowest=TRUE)
    sim_out <- as.data.frame(table(factor_sim))
    sim_factor[[i]] <- transform(sim_out, cumFreq = cumsum(sim_out$Freq), relative = prop.table(sim_out$Freq))
  }

  Frec_rel_sim <- list()
  for(i in 1:(dim(df)[2]-1)) {
    Frec_rel_sim[[i]] <- sim_factor[[i]][,'Freq']/sum(obs_out$Freq)
  }

  sim_out <- list()
  for(i in 1:(dim(df)[2]-1)) {
    sim_out[[i]] <- (obs_out[,'relative'] - Frec_rel_sim[[i]])^2/Frec_rel_sim[[i]]
    sim_out[[i]][is.nan(sim_out[[i]])] = 0
  }

  chi_calc <- list()
  for(i in 1:(dim(df)[2]-1)) {
    chi_calc[[i]] <- sum(obs_out[,'Freq'])*sum(sim_out[[i]],na.rm = TRUE)
  }

  # Level of significance alpha (tipical 5%).
  # Confidence level is 1 - alpha
  degree_freedom <- nbreaks - dist_param - 1
  theoretical_chi <- stats::qchisq(1 - alpha, degree_freedom)

  test_chi_hidrol <- list()
  for(i in 1:(dim(df)[2]-1)) {
    test_chi_hidrol[[i]] <- chi_calc[[i]] <= theoretical_chi
  }

  Chi_C_Chi_T <- list()
  for(i in 1:(dim(df)[2]-1)) {
    Chi_C_Chi_T[[i]] <- as.numeric(chi_calc[i])/theoretical_chi
  }

  distrib_selec <- min(as.numeric(chi_calc)/theoretical_chi)
  aux1 <- colnames(df)[which(Chi_C_Chi_T[]==distrib_selec)+1]
  print(paste0("Distribution of minimal Chi-Calculated / Chi-Theorical = ",aux1))

  aux2 <- Chi_C_Chi_T[Chi_C_Chi_T[]==distrib_selec]
  print(paste0("Minimal Chi-Calculated / Chi-Theorical = ",aux2))

  distrib_calc_theo <- as.numeric(chi_calc)/theoretical_chi
  distrib_calc <- as.numeric(chi_calc)
  distrib_theo <- rep(theoretical_chi,length(chi_calc))
  distrib_calc_theo <- data.frame(cbind(Distribution = colnames(df)[2:dim(df)[2]],
                                        chi_calc = round(distrib_calc,3),
                                        chi_theo = round(distrib_theo,3),
                                        calc_theo = round(distrib_calc_theo,3)))

  print("Chi-Calculated / Chi-Theorical for each distribution = ")

  print(distrib_calc_theo)

  marca_de_clase <- brks + int/2

  # Create data
  data <- data.frame(
    name = marca_de_clase[1:length(marca_de_clase) - 1] ,
    value = obs_out[,"Freq"]
  )

  df_aux <- data.frame(matrix(unlist(sim_factor), ncol=length(sim_factor), byrow=F))[(nbreaks+1):(2*nbreaks),]
  colnames(df_aux) <- colnames(df)[2:(dim(df)[2])]
  rownames(df_aux) <- seq(1,nbreaks,1)
  df_aux <- data.frame(name = data$name, df_aux)
  df_aux <- reshape2::melt(df_aux,  id.vars = 'name', variable.name = 'Distrib')

  # Barplot
  plot1 <- ggplot2::ggplot(data, aes(x=data$name, y=data$value, fill="gold2")) +
    geom_bar(colour="black", stat = "identity", width=5) +
    xlab("Variable") + ylab("Absolute Freq.") +
    theme_bw() + scale_fill_manual(values="gold2", labels="Real Frec.") +
    geom_line(data = df, aes(colour = df$Distrib), lwd = 1) +
    geom_point(data = df, aes(colour = df$Distrib), shape = 5, size = 1.5)

  table_chi <- t(distrib_calc_theo)
  tbl <- ggtexttable(table_chi)

  ggarrange(plot1, tbl,
            ncol = 1, nrow = 2,
            heights = c(1.8,1))

}
