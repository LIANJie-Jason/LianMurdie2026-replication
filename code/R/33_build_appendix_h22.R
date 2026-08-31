# Rebuild Appendix C7-C14, including stable-ID C9/C12 quadratic refits.
.root<-Sys.getenv("REPLICATION_ROOT",unset="");if(!nzchar(.root)){.a<-grep("^--file=",commandArgs(FALSE),value=TRUE);.root<-file.path(dirname(normalizePath(sub("^--file=","",.a[1]))),"..","..")};source(file.path(.root,"code","R","_table_helpers.R"));P<-replication_bootstrap();source(file.path(P$code,"_robustness_helpers.R"))
suppressPackageStartupMessages({library(survival);library(logistf)});spec<-function(n)file.path(P$code,"specs","appendix",n);out<-function(n)file.path(P$appendix_tables,n)
m<-read_models(P,"H22");s<-read.csv(file.path(P$estimates,"H22_results.csv"))
emit<-function(name,ids,n,scaffold=name)write_body(hydrate_scaffold(spec(scaffold),m[ids],s),out(name),n)
emit("H22_concession_campaign_table_body.csv",c("M1","M5","M2","M6","M3","M7","M4","M8"),34L)
emit("H22_strict_success_campaign_table_body.csv",c("M1s","M5s","M2s","M6s","M3s","M7s","M4s","M8s"),34L,"H22_concession_campaign_table_body.csv")
emit("H22_concession_campaign_year_table_body.csv",c("M13","M17","M14","M18","M15","M19","M16","M20"),34L)
emit("H22_strict_success_campaign_year_table_body.csv",c("M13s","M17s","M14s","M18s","M15s","M19s","M16s","M20s"),34L,"H22_concession_campaign_year_table_body.csv")
emit("H22_rest_ordered_table_body.csv",c("M9","M21","M10","M22","M11","M23","M12","M24"),34L)
emit("H22_rest_robustness_table_body.csv",c("M13_AG","M17_GEE","M14_AG","M18_GEE","M15_AG","M19_GEE","M16_AG","M20_GEE"),34L,"H22_concession_campaign_year_table_body.csv")

D<-prepare_robustness_data(P);df<-D$df;d21<-D$d21;df$duration<-df$EYEAR-df$BYEAR+1;df$event_concession<-df$concession_bin;d21$tstart<-d21$time_in_campaign-1;d21$tstop<-d21$time_in_campaign;d21$event_concession<-d21$concession_bin
d21_first<-d21[order(d21$navco21_id,d21$time_in_campaign),];d21_first<-d21_first|>dplyr::group_by(navco21_id)|>dplyr::mutate(.cum=cumsum(event_concession),.prior=dplyr::lag(.cum,default=0L))|>dplyr::filter(.prior==0L)|>dplyr::mutate(event_concession=as.integer(event_concession==1&.cum==1L))|>dplyr::ungroup()
mods<-c(FDI="fdiin_c",Trade="trade_c",Aid="lnaid_c",IODonor="iodonor_c");q9<-list();q12<-list()
for(k in seq_along(mods)){v<-unname(mods[k]);rhs<-paste0("eng_c+I(eng_c^2)+",v,"+eng_c:",v,"+I(eng_c^2):",v)
 fw<-as.formula(paste("Surv(duration,event_concession)~",rhs,"+nonviolent_camp+REGCHANGE+v2csreprss+lnnum_image_sum+lngdp"));ff<-as.formula(paste("concession_bin~",rhs,"+nonviolent_camp+REGCHANGE+v2csreprss+lnlengthofcam+lnnum_image_sum+lngdp"));fc<-as.formula(paste("Surv(tstart,tstop,event_concession)~",rhs,"+nonviolent+REGCHANGE+v2csreprss+lnnum_image_sum+lngdp_yr+cluster(navco21_id)"));fy<-as.formula(paste("concession_bin~",rhs,"+nonviolent+REGCHANGE+v2csreprss+time_in_campaign+lnnum_image_sum+lngdp_yr"))
 q9[[sprintf("C9_%02d_W",k)]]<-survreg(fw,data=df,dist="weibull");q9[[sprintf("C9_%02d_F",k)]]<-logistf(ff,data=df,control=logistf.control(maxit=5000,maxstep=.25),plcontrol=logistpl.control(maxit=5000),pl=TRUE);q12[[sprintf("C12_%02d_C",k)]]<-coxph(fc,data=d21_first);q12[[sprintf("C12_%02d_F",k)]]<-logistf(fy,data=d21,control=logistf.control(maxit=5000,maxstep=.25),plcontrol=logistpl.control(maxit=5000),pl=TRUE)}
saveRDS(c(q9,q12),file.path(P$cache,"H22_quadratic_models.rds"),version=3)
write.csv(model_terms_long(q9,"C9"),file.path(P$estimates,"H22_quadratic_campaign_results.csv"),row.names=FALSE,na="")
write.csv(model_terms_long(q12,"C12"),file.path(P$estimates,"H22_quadratic_campaign_year_results.csv"),row.names=FALSE,na="")
write_body(hydrate_scaffold(spec("H22_quadratic_campaign_table_body.csv"),q9),out("H22_quadratic_campaign_table_body.csv"),44L)
write_body(hydrate_scaffold(spec("H22_quadratic_campaign_year_table_body.csv"),q12),out("H22_quadratic_campaign_year_table_body.csv"),44L)
