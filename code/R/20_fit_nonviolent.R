# Fit the accepted C22 nonviolent-only battery, including NVH1 exactly once.
.root <- Sys.getenv("REPLICATION_ROOT", unset=""); if(!nzchar(.root)){.a<-grep("^--file=",commandArgs(FALSE),value=TRUE);.root<-file.path(dirname(normalizePath(sub("^--file=","",.a[1]))),"..","..")}
source(file.path(.root,"code","R","_table_helpers.R")); P <- replication_bootstrap()
source(file.path(P$code,"_robustness_helpers.R"))
suppressPackageStartupMessages(library(geepack)); D <- prepare_robustness_data(P); df<-D$df; d21<-D$d21
df_nv <- df[as.numeric(as.character(df$nonviolent_camp))==1,]; d21_nv <- d21[as.numeric(as.character(d21$nonviolent))==1,]
rows <- list(); add <- function(x) rows[[length(rows)+1L]] <<- x

gee_long <- function(model,id,moderator) {
  tab<-summary(model)$coefficients; err<-model$geese$error; ok<-!is.null(err)&&as.integer(err)==0L
  vals<-c(coef(model),tab[,"Std.err"]); status<-if(!ok) "gee_error_not_reported" else if(!all(is.finite(vals))) "nonfinite_gee_not_reported" else if(max(abs(vals))>=1e4) "unstable_gee_not_reported" else "reported"
  report<-identical(status,"reported"); est<-unname(coef(model)); se<-tab[names(coef(model)),"Std.err"]; p<-tab[names(coef(model)),"Pr(>|W|)"]
  if(!report){est[]<-NA;se[]<-NA;p[]<-NA}
  data.frame(Model=id,Estimator="GEE logit",DV="concession (nonviolent, 2.1 clustered GEE companion)",Moderator=moderator,
    Variable=names(coef(model)),Estimate=est,SE=se,p=p,CI_lo=est-1.96*se,CI_hi=est+1.96*se,
    N=if(!is.null(model$model))nrow(model$model)else length(model$y),Converged=ok,Fit_Status=status,stringsAsFactors=FALSE)
}
fit_gee <- function(f){v<-unique(c(all.vars(f),"navco21_id","time_in_campaign"));z<-d21_nv[complete.cases(d21_nv[,v]),v];z<-z[order(z$navco21_id,z$time_in_campaign),];geeglm(f,data=z,id=navco21_id,family=binomial,corstr="exchangeable")}

# The sole nonviolent H1 focal row.
nvh1<-fit_firth(concession_bin~eng_c+I(eng_c^2)+REGCHANGE+v2csreprss+lnlengthofcam+lnnum_image_sum+lnpop+colonized_english,df_nv)
add(long_logistf(nvh1,"NVH1","English share squared","concession (nonviolent)"))

pol <- list(list("NV1","NV3","domaut_c","domaut"),list("NV2","NV4","stterr_c","stterr"))
for(s in pol){f13<-as.formula(paste("concession_bin~eng_c*",s[[3]],"+REGCHANGE+v2csreprss+lnlengthofcam+lnnum_image_sum"));f21<-as.formula(paste("concession_bin~eng_c*",s[[3]],"+REGCHANGE+v2csreprss+time_in_campaign+lnnum_image_sum"));add(long_logistf(fit_firth(f13,df_nv),s[[1]],s[[4]],"concession (nonviolent)"));add(long_logistf(fit_firth(f21,d21_nv),s[[2]],s[[4]],"concession (nonviolent, 2.1 pooled-year Firth diagnostic)"));add(gee_long(fit_gee(f21),paste0(s[[2]],"_GEE"),s[[4]]))}

econ <- list(list("fdiin_c","FDI inflow","full"),list("trade_c","Trade openness","full"),list("lnaid_c","ln(Aid)","reduced"),list("iodonor_c","IO donor count","reduced"))
for(j in seq_along(econ)){s<-econ[[j]];ctrl13<-if(s[[3]]=="full")"goals+v2csreprss+lngdp+lnpop+colonized_english+lnlengthofcam+lnnum_image_sum" else "REGCHANGE+v2csreprss+lngdp+lnlengthofcam+lnnum_image_sum";f13<-as.formula(paste("concession_bin~eng_c*",s[[1]],"+",ctrl13));f21<-as.formula(paste("concession_bin~eng_c*",s[[1]],"+REGCHANGE+v2csreprss+lngdp_yr+time_in_campaign+lnnum_image_sum"));add(long_logistf(fit_firth(f13,df_nv),paste0("NV",j+4),s[[2]],"concession (nonviolent)"));add(long_logistf(fit_firth(f21,d21_nv),paste0("NV",j+8),s[[2]],"concession (nonviolent, 2.1 pooled-year Firth diagnostic)"));add(gee_long(fit_gee(f21),paste0("NV",j+8,"_GEE"),s[[2]]))}

dom<-list(list("repress_c","CSO repression",TRUE),list("civlib_c","Civil liberties",FALSE),list("clpol_c","Political liberty",FALSE),list("polity_c","Polity",FALSE))
for(j in seq_along(dom)){s<-dom[[j]];ctrl13<-if(s[[3]])"goals+p_polity2+lngdp+lnpop+colonized_english+lnlengthofcam+lnnum_image_sum" else "goals+v2csreprss+lngdp+lnpop+colonized_english+lnlengthofcam+lnnum_image_sum";ctrl21<-if(s[[3]])"REGCHANGE+lngdp_yr+time_in_campaign+lnnum_image_sum" else "REGCHANGE+v2csreprss+lngdp_yr+time_in_campaign+lnnum_image_sum";f13<-as.formula(paste("concession_bin~eng_c*",s[[1]],"+",ctrl13));f21<-as.formula(paste("concession_bin~eng_c*",s[[1]],"+",ctrl21));add(long_logistf(fit_firth(f13,df_nv),paste0("NV",j+12),s[[2]],"concession (nonviolent)"));add(long_logistf(fit_firth(f21,d21_nv),paste0("NV",j+16),s[[2]],"concession (nonviolent, 2.1 pooled-year Firth diagnostic)"));add(gee_long(fit_gee(f21),paste0("NV",j+16,"_GEE"),s[[2]]))}

out<-dplyr::bind_rows(rows)
model_order<-c("NVH1","NV1","NV2","NV3","NV3_GEE","NV4","NV4_GEE","NV5","NV9","NV9_GEE","NV6","NV10","NV10_GEE","NV7","NV11","NV11_GEE","NV8","NV12","NV12_GEE","NV13","NV17","NV17_GEE","NV14","NV18","NV18_GEE","NV15","NV19","NV19_GEE","NV16","NV20","NV20_GEE")
out$.order<-match(out$Model,model_order);out<-out[order(out$.order,seq_len(nrow(out))),];out$.order<-NULL
if(nrow(out)!=289L)stop("C22 long results must have 289 rows; found ",nrow(out),call.=FALSE)
if(sum(out$Model=="NVH1" & out$Variable=="I(eng_c^2)",na.rm=TRUE)!=1L)stop("Accepted C22 contract violation: NVH1 focal row must occur exactly once.",call.=FALSE)
write.csv(out,file.path(P$estimates,"appendix_nonviolent_subsample_results.csv"),row.names=FALSE,na="")
saveRDS(list(NVH1=nvh1),file.path(P$cache,"nonviolent_focal_models.rds"))
