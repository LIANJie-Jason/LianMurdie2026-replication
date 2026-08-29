# Private data preparation and long-result extractors for C21-C24 fits.

prepare_robustness_data <- function(P) {
  suppressPackageStartupMessages({library(dplyr); library(logistf); library(brglm2)})
  df <- read.csv(require_file(file.path(P$data, "df_final.csv")))
  d21 <- read.csv(require_file(file.path(P$data, "df_navco21_panel.csv")))
  needed13 <- c("CAMPAIGN","LOCATION","success","limited","ongoing","nonviolent_camp","REGCHANGE","v2csreprss",
    "lnlengthofcam","lnnum_image_sum","lnpop","colonized_english","eng_prob_general","v2svdomaut","v2svstterr",
    "wdi_fdiin","fi_ftradeint_pd","aid_crnio","lnaid","lngdp","p_polity2","v2x_civlib","v2x_clpol","success3")
  needed21 <- c("navco21_id","Year","nonviolent","camp_goals","Num_Image","progress","concession","time_in_campaign",
    "eng_prop_year","v2csreprss","lngdp_yr","aid_crnio","v2svdomaut","v2svstterr","wdi_fdiin","fi_ftradeint_pd",
    "lnaid_yr","p_polity2","v2x_civlib","v2x_clpol")
  miss13 <- setdiff(needed13,names(df)); miss21 <- setdiff(needed21,names(d21))
  if(length(miss13)) stop("df_final missing: ",paste(miss13,collapse=", "),call.=FALSE)
  if(length(miss21)) stop("df_navco21_panel missing: ",paste(miss21,collapse=", "),call.=FALSE)
  if(anyDuplicated(df[c("CAMPAIGN","LOCATION")])) stop("Duplicate NAVCO 1.3 campaign-location keys.",call.=FALSE)
  if(anyDuplicated(d21[c("navco21_id","Year")])) stop("Duplicate NAVCO 2.1 id-year keys.",call.=FALSE)
  df$nonviolent_camp <- as.factor(df$nonviolent_camp); df$REGCHANGE <- as.numeric(as.character(df$REGCHANGE))
  df$goals <- as.factor(df$goals); df$colonized_english <- as.factor(df$colonized_english); df$success3 <- ordered(df$success3)
  df$concession_bin <- as.integer(as.numeric(as.character(df$success))==1 | as.numeric(as.character(df$limited))==1)
  df$concession <- df$concession_bin; df$success_bin <- as.integer(as.numeric(as.character(df$success))==1)
  centers13 <- c(eng_c="eng_prob_general",domaut_c="v2svdomaut",stterr_c="v2svstterr",fdiin_c="wdi_fdiin",
    trade_c="fi_ftradeint_pd",lnaid_c="lnaid",iodonor_c="aid_crnio",repress_c="v2csreprss",civlib_c="v2x_civlib",
    clpol_c="v2x_clpol",polity_c="p_polity2")
  for(nm in names(centers13)) df[[nm]] <- df[[centers13[[nm]]]]-mean(df[[centers13[[nm]]]],na.rm=TRUE)
  d21$nonviolent <- as.factor(d21$nonviolent); d21$REGCHANGE <- as.numeric(d21$camp_goals==0)
  d21$concession_bin <- as.numeric(d21$concession); d21$lnnum_image_sum <- log(d21$Num_Image)
  centers21 <- c(eng_c="eng_prop_year",domaut_c="v2svdomaut",stterr_c="v2svstterr",fdiin_c="wdi_fdiin",
    trade_c="fi_ftradeint_pd",lnaid_c="lnaid_yr",iodonor_c="aid_crnio",repress_c="v2csreprss",civlib_c="v2x_civlib",
    clpol_c="v2x_clpol",polity_c="p_polity2")
  for(nm in names(centers21)) d21[[nm]] <- d21[[centers21[[nm]]]]-mean(d21[[centers21[[nm]]]],na.rm=TRUE)
  list(df=df,d21=d21)
}

long_logistf <- function(model,id,moderator,dv) {
  se <- sqrt(diag(vcov(model))); nms <- names(coef(model))
  data.frame(Model=id,Estimator="Firth PML",DV=dv,Moderator=moderator,Variable=nms,Estimate=unname(coef(model)),
    SE=unname(se),p=unname(model$prob),CI_lo=unname(model$ci.lower),CI_hi=unname(model$ci.upper),N=model$n,
    Converged=NA,Fit_Status=NA_character_,stringsAsFactors=FALSE)
}

safe_bracl <- function(formula,data) tryCatch(suppressWarnings(bracl(formula,data=data,type="MPL_Jeffreys",link="logit",parallel=TRUE,maxit=500)),error=function(e) NULL)
long_bracl <- function(model,id,moderator,n,dv) {
  if(is.null(model)) return(data.frame(Model=id,Estimator="bracl (MPL_Jeffreys)",DV=dv,Moderator=moderator,Variable="(FAILED)",Estimate=NA,SE=NA,p=NA,CI_lo=NA,CI_hi=NA,N=n,Converged=FALSE,Fit_Status=NA_character_))
  tab <- summary(model)$coefficients; ci <- confint(model,level=.95,type="Wald"); ok <- isTRUE(model$converged)
  data.frame(Model=id,Estimator="bracl (MPL_Jeffreys)",DV=dv,Moderator=moderator,Variable=rownames(tab),
    Estimate=if(ok) tab[,"Estimate"] else NA,SE=if(ok) tab[,"Std. Error"] else NA,p=if(ok) tab[,"Pr(>|z|)"] else NA,
    CI_lo=if(ok) ci[rownames(tab),1] else NA,CI_hi=if(ok) ci[rownames(tab),2] else NA,N=n,Converged=ok,Fit_Status=NA_character_,stringsAsFactors=FALSE)
}

fit_firth <- function(formula,data) logistf(formula,data=data,na.action=na.omit,
  control=logistf.control(maxit=1000,maxstep=.5),plcontrol=logistpl.control(maxit=1000),pl=TRUE)
