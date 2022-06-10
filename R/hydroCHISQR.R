#' Chi - Squared test for hydrological frequency analysis
#' @export
#' @param df Data frame with observations and distributions adjusted
#' @param nc Number of classes, default 1 for Sturges, 2 for Scott, 3 for Freedman-Diaconis
#' @param dist_param Number of parameters of the distributions in df
#' @param alpha Significance level, typical 0.05
#' @param c_test Statistical result or plot, default 1 for statistical and for 2 plot
#' @param c_method Chi - Squared test method: 1 (Varas & Bois, 1998) and 2 (Chow, 1949): e(i) = F(i-1) - F(i), to perform Chi(i) = (e(i) - fi)^2/e(i)
#' @import grDevices ggplot2 reshape2 stats ggpubr FAmle

hydroCHISQR <- function(df, nc, dist_param, alpha, c_test, c_method) {
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


  if (c_method == 1) {

  # Chi-Squared method 1
  sim_out <- list()
  for(i in 1:(dim(df)[2]-1)) {
    sim_out[[i]] <- (obs_out[,'relative'] - Frec_rel_sim[[i]])^2/Frec_rel_sim[[i]]
    sim_out[[i]][is.nan(sim_out[[i]])] = 0
  }

  chi_calc <- list()
  for(i in 1:(dim(df)[2]-1)) {
    chi_calc[[i]] <- sum(obs_out[,'Freq'])*sum(sim_out[[i]],na.rm = TRUE)
  }

  } else {

    # Chi-Squared method 2
    # F(i)
    F_i <- list()
    for(i in 1:(dim(df)[2]-1)) F_i[[i]] <- FAmle::distr(x=brks,model = colnames(df)[2:(dim(df)[2])] ,type='p')

    e_i <- as.data.frame(F_i)
    e_i <- e_i[2:dim(e_i)[1],] - e_i[1:dim(e_i)[1]-1,]

    sim_out <- list()
    for(i in 1:(dim(df)[2]-1)) {
      sim_out[[i]] <- (obs_out[,'relative'] - e_i[[i]])^2/e_i[[i]]
      sim_out[[i]][is.nan(sim_out[[i]])] = 0
    }

    chi_calc <- list()
    for(i in 1:(dim(df)[2]-1)) {
      chi_calc[[i]] <- sum(obs_out[,'Freq'])*sum(sim_out[[i]],na.rm = TRUE)
    }

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

  aux2 <- Chi_C_Chi_T[Chi_C_Chi_T[]==distrib_selec]

  distrib_calc_theo <- as.numeric(chi_calc)/theoretical_chi
  distrib_calc <- as.numeric(chi_calc)
  distrib_theo <- rep(theoretical_chi,length(chi_calc))
  distrib_calc_theo <- data.frame(cbind(Distribution = colnames(df)[2:dim(df)[2]],
                                        chi_calc = round(distrib_calc,3),
                                        chi_theo = round(distrib_theo,3),
                                        calc_theo = round(distrib_calc_theo,3)))

  distrib_calc_theo$test <- with(distrib_calc_theo, ifelse(calc_theo < 1, "Accepted", "Rejected"))

  marca_de_clase <- brks + int/2

  if (c_test == 2) {
    # Create data
    data <- data.frame(
      name = marca_de_clase[1:length(marca_de_clase) - 1] ,
      value = obs_out[,"Freq"]
      )

    df_aux <- data.frame(matrix(unlist(sim_factor), ncol=length(sim_factor), byrow=F))[(nbreaks+1):(2*nbreaks),]
    colnames(df_aux) <- colnames(df)[2:(dim(df)[2])]
    rownames(df_aux) <- seq(1,nbreaks,1)
    df_aux <- data.frame(name = data$name, df_aux)
    df_aux <- reshape2::melt(df_aux, id.vars = 'name', variable.name = 'Distrib')

    # Barplot
    plot1 <- ggplot2::ggplot(data, aes(x=name, y=value, fill="gold2")) +
      geom_bar(colour="black", stat = "identity", width=5) +
      xlab("Variable") + ylab("Absolute Freq.") +
      theme_bw() + scale_fill_manual(values="gold2", labels="Real Frec.") +
      geom_line(data = df_aux, aes(colour = Distrib), lwd = 1) +
      geom_point(data = df_aux, aes(colour = Distrib), shape = 5, size = 1.5)

    table_chi <- t(distrib_calc_theo)
    tbl <- ggtexttable(table_chi)

    ggarrange(plot1, tbl,
            ncol = 1, nrow = 2,
            heights = c(1.8,1))
  } else {
    print("Chi-Calculated / Chi-Theorical for each distribution = ")
    print(distrib_calc_theo)
    print(paste0("Distribution of minimal Chi-Calculated / Chi-Theorical = ",aux1))
    print(paste0("Minimal Chi-Calculated / Chi-Theorical = ",aux2))
  }

}
