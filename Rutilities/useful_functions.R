numerize <- function(data, columNames){
  for (col in columNames){
    data[[which(names(data)==col)]] <- as.character(data[[which(names(data)==col)]])
    data[[which(names(data)==col)]] <- as.numeric(data[[which(names(data)==col)]])
  }
  return(data)
}

deFactorize <- function(data, columNames){
  for (col in columNames){
    data[[which(names(data)==col)]] <- as.character(data[[which(names(data)==col)]])
  }
  return(data)
}

compare<-function(data,x,y){
  a = as.character(data[[x]]) 
  b = as.character(data[[y]])
  comparison <- ifelse (is.na(a) & is.na(b), "both-missing", 
                        ifelse(is.na(a), "missing-in-A",
                          ifelse(is.na(b), "missing-in-B",
                          ifelse(a==b, "concordant",
                               ifelse((a=="Heterozygous" & b=="Homozygous") | (a=="Homozygous" & b=="Heterozygous"), "discordant - het/homo",
                                      ifelse((a=="HeterozygousNonRef" & b=="Heterozygous") | (a=="Heterozygous" & b=="HeterozygousNonRef") | (a=="HeterozygousNonRef" & b=="Homozygous") | (a=="Homozygous" & b=="HeterozygousNonRef"), "discordant - nonRef", 
                                             "discordant - other")
                               )))))
  return(comparison)
}


## standard error for bar plots
stde <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))
## NA become zero
na.zero <- function(x) {replace(x, is.na(x), 0)}
## Get Tidy from Giulio
get_tidy <- function(model, term, exponentiate = FALSE) {
  estimate <- unname(coef(model)[term])
  std.error <- unname(sqrt(diag(vcov(model)))[term])
  statistic <- unname(coef(model)[term] / sqrt(diag(vcov(model)))[term])
  p.value <- min(pnorm(statistic), pnorm(statistic, lower.tail = FALSE)) * 2
  conf.low <- unname(estimate - 1.96 * std.error)
  conf.high <- unname(estimate + 1.96 * std.error)
  if (exponentiate) {
    return(list(term = term, estimate = exp(estimate), std.error = std.error, statistic = statistic, p.value = p.value, conf.low = exp(conf.low), conf.high = exp(conf.high)))
  } else {
    return(list(term = term, estimate = estimate, std.error = std.error, statistic = statistic, p.value = p.value, conf.low = conf.low, conf.high = conf.high))
  }
}

writeLines("NB: for the variant recovery
           the first column needs to be the tissue, as the recovery is computed from it")

variantRecovery <- function(data,x,y){
  a = as.character(data[[x]]) 
  b = as.character(data[[y]])
  recovery = ifelse(is.na(a) & is.na(b), "both-missing",
                    ifelse(is.na(a), "missing-in-tissue",
                      ifelse(is.na(b), "not-recovered_(missing-in-plasma)",
                        ifelse(a==b, "recovered_(concordant)",
                          ifelse((a=="Heterozygous" & b=="Homozygous") | (a=="Homozygous" & b=="Heterozygous"), "recovered_(discordant-het/homo)",
                            ifelse((a=="HeterozygousNonRef" & b=="Heterozygous") | (a=="Heterozygous" & b=="HeterozygousNonRef") | (a=="HeterozygousNonRef" & b=="Homozygous") | (a=="Homozygous" & b=="HeterozygousNonRef"), "recovered_(discordant-nonRef)",
                              "recovered_(discordant-other)"
                            )
                          )
                        )
                      )
                    )
                  )
}



compareSoftware<-function(data,x,y){
  a = as.character(data[[x]]) 
  b = as.character(data[[y]])
  comparison <- ifelse(is.na(a) & is.na(b), "both-missing",
                       ifelse(is.na(a), "missing-in-geneGlobe",
                         ifelse(is.na(b), "missing-in-BxWB",
                           ifelse(a==b, "called-in-both","other")
                         )
                       )
                )
  return(comparison)
}
