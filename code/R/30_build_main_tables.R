# Rebuild accepted main-paper Tables 1-4 from freshly fitted model objects.
.root <- Sys.getenv("REPLICATION_ROOT", unset=""); if(!nzchar(.root)){.a<-grep("^--file=",commandArgs(FALSE),value=TRUE);.root<-file.path(dirname(normalizePath(sub("^--file=","",.a[1]))),"..","..")}
source(file.path(.root,"code","R","_table_helpers.R"))
P <- replication_bootstrap()

make_spec <- function(labels, columns) {
  out <- data.frame(label = labels, stringsAsFactors = FALSE)
  for (nm in names(columns)) out[[nm]] <- columns[[nm]]
  out
}

h1 <- read_models(P, "H1")
h1_ids <- c("M1", "M2", "M4", "M5")
h1_labels <- c("Proportion of English-Language Signs Squared", "Proportion of English-Language Signs",
  "Campaign Primarily Nonviolent", "Campaign Purpose Regime Change", "Campaign Location V-Dem CSO Repression",
  "Campaign Length, in Years (ln)", "Time in Campaign", "Number of Images on Campaign (ln)",
  "Number of Images in Campaign-Year (ln)", "Campaign Location Population (ln)", "Campaign-Year Population (ln)",
  "English Colonial Legacy", "Constant")
h1_spec <- make_spec(h1_labels, list(
  M1=c("I(eng_prob_general^2)","eng_prob_general","nonviolent_camp1","REGCHANGE","v2csreprss",NA,NA,"lnnum_image_sum",NA,"lnpop",NA,"colonized_english1","(Intercept)"),
  M2=c("I(eng_prob_general^2)","eng_prob_general","nonviolent_camp1","REGCHANGE","v2csreprss","lnlengthofcam",NA,"lnnum_image_sum",NA,"lnpop",NA,"colonized_english1","(Intercept)"),
  M4=c("I(eng_prop_year^2)","eng_prop_year","nonviolent1","REGCHANGE","v2csreprss",NA,NA,NA,"lnnum_image_sum",NA,"lnpop_yr","colonized_english1",NA),
  M5=c("I(eng_prop_year^2)","eng_prop_year","nonviolent1","REGCHANGE","v2csreprss",NA,"time_in_campaign",NA,"lnnum_image_sum",NA,"lnpop_yr","colonized_english1","(Intercept)")))
write_body(build_coefficient_body(h1_spec, h1[h1_ids]), file.path(P$main_tables, "H1_main_table_body.csv"), 28L)

h21 <- read_models(P, "H21")
h21_ids <- c("M1","M2","M3","M4","M7","M8","M9","M10")
h21_labels <- c("Proportion of English-Language Signs","Domestic Autonomy","Proportion English-Language Signs X Domestic Autonomy",
  "State Authority Over Territory","Proportion English-Language Signs X State Authority Over Territory","Campaign Length, in Years (ln)",
  "Time in Campaign (years)","Campaign Primarily Nonviolent","Campaign Purpose Regime Change","Campaign Location V-Dem CSO Repression",
  "Number of Images on Campaign (ln)","Number of Images in Campaign-Year (ln)","Constant")
base_h21 <- list(
 M1=c("eng_c","domaut_c","eng_c:domaut_c",NA,NA,NA,NA,"nonviolent_camp1","REGCHANGE","v2csreprss","lnnum_image_sum",NA,"(Intercept)"),
 M2=c("eng_c",NA,NA,"stterr_c","eng_c:stterr_c",NA,NA,"nonviolent_camp1","REGCHANGE","v2csreprss","lnnum_image_sum",NA,"(Intercept)"),
 M3=c("eng_c","domaut_c","eng_c:domaut_c",NA,NA,"lnlengthofcam",NA,"nonviolent_camp1","REGCHANGE","v2csreprss","lnnum_image_sum",NA,"(Intercept)"),
 M4=c("eng_c",NA,NA,"stterr_c","eng_c:stterr_c","lnlengthofcam",NA,"nonviolent_camp1","REGCHANGE","v2csreprss","lnnum_image_sum",NA,"(Intercept)"),
 M7=c("eng_c","domaut_c","eng_c:domaut_c",NA,NA,NA,NA,"nonviolent1","REGCHANGE","v2csreprss",NA,"lnnum_image_sum",NA),
 M8=c("eng_c",NA,NA,"stterr_c","eng_c:stterr_c",NA,NA,"nonviolent1","REGCHANGE","v2csreprss",NA,"lnnum_image_sum",NA),
 M9=c("eng_c","domaut_c","eng_c:domaut_c",NA,NA,NA,"time_in_campaign","nonviolent1","REGCHANGE","v2csreprss",NA,"lnnum_image_sum","(Intercept)"),
 M10=c("eng_c",NA,NA,"stterr_c","eng_c:stterr_c",NA,"time_in_campaign","nonviolent1","REGCHANGE","v2csreprss",NA,"lnnum_image_sum","(Intercept)"))
write_body(build_coefficient_body(make_spec(h21_labels, base_h21), h21[h21_ids]), file.path(P$main_tables, "H21_main_table_body.csv"), 28L)

econ_labels <- c("Proportion of English-Language Signs","Campaign Location FDI Inflows (% of GDP)","Proportion English-Language Signs X FDI Inflows (% of GDP)",
 "International Trade Freedom","Proportion English-Language Signs X International Trade Freedom","Foreign Aid (ln)",
 "Proportion English-Language Signs X Foreign Aid (ln)","Number of IO Donors","Proportion English-Language Signs X Number of IO Donors",
 "Campaign Length, in Years (ln)","Campaign Primarily Nonviolent","Campaign Purpose Regime Change","Campaign Location V-Dem CSO Repression",
 "Number of Images on Campaign (ln)","Campaign Location GDP (ln)","Constant")
econ_col <- function(mod, int, binary=FALSE) c("eng_c", if(mod=="fdiin_c") mod else NA, if(mod=="fdiin_c") int else NA,
 if(mod=="trade_c") mod else NA,if(mod=="trade_c") int else NA,if(mod=="lnaid_c") mod else NA,if(mod=="lnaid_c") int else NA,
 if(mod=="iodonor_c") mod else NA,if(mod=="iodonor_c") int else NA,if(binary) "lnlengthofcam" else NA,"nonviolent_camp1","REGCHANGE","v2csreprss","lnnum_image_sum","lngdp","(Intercept)")
h22 <- read_models(P, "H22"); h22_ids <- c("M1","M5","M2","M6","M3","M7","M4","M8")
h22_cols <- list(M1=econ_col("fdiin_c","eng_c:fdiin_c"),M5=econ_col("fdiin_c","eng_c:fdiin_c",TRUE),M2=econ_col("trade_c","eng_c:trade_c"),M6=econ_col("trade_c","eng_c:trade_c",TRUE),M3=econ_col("lnaid_c","eng_c:lnaid_c"),M7=econ_col("lnaid_c","eng_c:lnaid_c",TRUE),M4=econ_col("iodonor_c","eng_c:iodonor_c"),M8=econ_col("iodonor_c","eng_c:iodonor_c",TRUE))
write_body(build_coefficient_body(make_spec(econ_labels,h22_cols),h22[h22_ids]),file.path(P$main_tables,"H22_main_table_body.csv"),34L)

dom_labels <- c("Proportion of English-Language Signs","Campaign Location V-Dem CSO Repression","Proportion English-Language Signs X V-Dem CSO Repression",
 "V-Dem Civil Liberty Index","Proportion English-Language Signs X V-Dem Civil Liberty Index","V-Dem Political Liberty Index",
 "Proportion English-Language Signs X V-Dem Political Liberty","Campaign Location Polity Score","Proportion English-Language Signs X Polity",
 "Campaign Length, in Years (ln)","Campaign Primarily Nonviolent","Campaign Purpose Regime Change","Number of Images on Campaign (ln)","Campaign Location GDP (ln)","Constant")
dom_col <- function(mod,int,binary=FALSE) c("eng_c",if(mod=="repress_c") mod else "v2csreprss",if(mod=="repress_c") int else NA,
 if(mod=="civlib_c") mod else NA,if(mod=="civlib_c") int else NA,if(mod=="clpol_c") mod else NA,if(mod=="clpol_c") int else NA,
 if(mod=="polity_c") mod else if(mod=="repress_c") "p_polity2" else NA,if(mod=="polity_c") int else NA,
 if(binary) "lnlengthofcam" else NA,"nonviolent_camp1","REGCHANGE","lnnum_image_sum","lngdp","(Intercept)")
h3 <- read_models(P,"H3"); h3_ids <- c("M1","M5","M2","M6","M3","M7","M4","M8")
h3_cols <- list(M1=dom_col("repress_c","eng_c:repress_c"),M5=dom_col("repress_c","eng_c:repress_c",TRUE),M2=dom_col("civlib_c","eng_c:civlib_c"),M6=dom_col("civlib_c","eng_c:civlib_c",TRUE),M3=dom_col("clpol_c","eng_c:clpol_c"),M7=dom_col("clpol_c","eng_c:clpol_c",TRUE),M4=dom_col("polity_c","eng_c:polity_c"),M8=dom_col("polity_c","eng_c:polity_c",TRUE))
write_body(build_coefficient_body(make_spec(dom_labels,h3_cols),h3[h3_ids]),file.path(P$main_tables,"H3_main_table_body.csv"),32L)
