SEER Sigmoidoscopy Project
================
<david.hein@utsouthwestern.edu>
2023-02-19

- <a href="#1-load-packages-and-raw-data"
  id="toc-1-load-packages-and-raw-data">1 Load packages and raw data</a>
- <a href="#2-data-cleaning" id="toc-2-data-cleaning">2 Data cleaning</a>
- <a href="#3-setting-reference-groups"
  id="toc-3-setting-reference-groups">3 Setting reference groups</a>
- <a href="#4-logistic-regression" id="toc-4-logistic-regression">4
  Logistic Regression</a>
  - <a href="#41-logistic-results-sex" id="toc-41-logistic-results-sex">4.1
    Logistic results: sex</a>
  - <a href="#42-logistic-results-age" id="toc-42-logistic-results-age">4.2
    Logistic results: age</a>
  - <a href="#43-logistic-results-stage"
    id="toc-43-logistic-results-stage">4.3 Logistic results: stage</a>
  - <a href="#44-logistic-results-year"
    id="toc-44-logistic-results-year">4.4 Logistic results: year</a>
  - <a href="#45-logistic-results-race"
    id="toc-45-logistic-results-race">4.5 Logistic results: race</a>
  - <a href="#46-logistic-results-multivariate"
    id="toc-46-logistic-results-multivariate">4.6 Logistic results:
    multivariate</a>
- <a href="#5-survival-exploration-wcss"
  id="toc-5-survival-exploration-wcss">5 Survival exploration w/CSS</a>
  - <a href="#51-km-plots-can-see-on-sigmoidoscopy"
    id="toc-51-km-plots-can-see-on-sigmoidoscopy">5.1 KM plots: Can see on
    sigmoidoscopy</a>
  - <a href="#52-km-plots-stage" id="toc-52-km-plots-stage">5.2 KM plots:
    Stage</a>
  - <a href="#53-km-plots-year-of-diag"
    id="toc-53-km-plots-year-of-diag">5.3 KM plots: year of diag</a>
  - <a href="#54-km-plots-sex" id="toc-54-km-plots-sex">5.4 KM plots:
    sex</a>
  - <a href="#55-km-plots-age-group" id="toc-55-km-plots-age-group">5.5 KM
    plots: age group</a>
  - <a href="#56-km-plots-raceeth" id="toc-56-km-plots-raceeth">5.6 KM
    plots: race/eth</a>
  - <a href="#57-km-plots-site" id="toc-57-km-plots-site">5.7 KM plots:
    Site</a>
- <a href="#6-survival-analysis" id="toc-6-survival-analysis">6 Survival
  Analysis</a>
- <a href="#7-km-plots-with-better-95-cis"
  id="toc-7-km-plots-with-better-95-cis">7 KM Plots with better 95%
  CIs</a>
- <a href="#8-css-plots" id="toc-8-css-plots">8 CSS Plots</a>
- <a href="#9-os-plots" id="toc-9-os-plots">9 OS plots</a>
- <a href="#10-session-info" id="toc-10-session-info">10 Session Info</a>

<br>

# 1 Load packages and raw data

``` r
library(broom)
library(MASS)
library(survival)
library(survminer)
library(tidyverse)

library(survRM2)
library(km.ci)
set.seed(2023)

seer_data <- read_delim("seer_data.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
```

<br>

# 2 Data cleaning

``` r
# Get only variables we want
new_data<-seer_data%>%dplyr::select(Sex, `Year of diagnosis`,
                                `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`,
                                `Primary Site`,`Grade (thru 2017)`,
                                `Grade Pathological (2018+)`,
                                `Combined Summary Stage (2004+)`,
                                `Summary stage 2000 (1998-2017)`,
                                `SEER cause-specific death classification`,
                                `Survival months`,
                                `Age recode with single ages and 100+`,
                                `Survival months flag`,
                                `Vital status recode (study cutoff used)`,
                                `Grade Clinical (2018+)`,
                                `Histologic Type ICD-O-3`)

# need to take out 181 appendix, 260, 188, and 189 becuase location unknown
new_data <-new_data %>% filter( !`Primary Site` %in% c(260,188,189,181))

# only adenocarcinoma
new_data <- new_data %>%filter(`Histologic Type ICD-O-3`==8140)

# make a new variable for site being seen on colonsocopy or sigmoid, should include 186 descending, 187 sigmoid, 199 rectosigmoid junction, 209 rectum NOS
new_data <- new_data %>% mutate(can_see_sigmoid = ifelse( `Primary Site` %in% c(186,187,199,209),1,0))

# make new variable for 3 age groups
new_data <- new_data %>% mutate(age_group = as.numeric(str_sub(`Age recode with single ages and 100+`,1,2)))
new_data <- new_data %>% mutate(age_group_final = ifelse(age_group <45,'under45', 'over50'))%>%
  mutate(age_group_final = ifelse(age_group <50 & age_group >= 45,'45-50',age_group_final))

# Combine the two stages
new_data<-new_data%>%mutate(new_stage = ifelse(`Summary stage 2000 (1998-2017)`=="Blank(s)",
                                               `Combined Summary Stage (2004+)`,`Summary stage 2000 (1998-2017)`))

# Filter out unknown stage an in situ
new_data2 <-new_data %>% filter( !new_stage %in%c("Unknown/unstaged","In situ"))

# Make nice grade variable
new_data2 <- new_data2 %>% mutate(final_grade = `Grade Pathological (2018+)`)
new_data2 <- new_data2 %>% mutate(final_grade = ifelse(`Grade (thru 2017)` %in% "Poorly differentiated; Grade III", "3", final_grade))
new_data2 <- new_data2 %>% mutate(final_grade = ifelse(`Grade (thru 2017)` %in% "Undifferentiated; anaplastic; Grade IV", "4", final_grade))
new_data2 <- new_data2 %>% mutate(final_grade = ifelse(`Grade (thru 2017)` %in% "Moderately differentiated; Grade II", "2", final_grade))
new_data2 <- new_data2 %>% mutate(final_grade = ifelse(`Grade (thru 2017)` %in% "Well differentiated; Grade I", "1", final_grade))

# toss em in an unknown/other
new_data2 <- new_data2 %>% mutate(final_grade = ifelse(final_grade%in%c("A", "B", "Blank(s)", "L","C", "H", "9","D"),"unk_other",final_grade))

# make smaller year of diagnosis groups
new_data2 <- new_data2 %>% mutate(year_of_diag = case_when(`Year of diagnosis`%in%c('2000', '2001', '2002' ,'2003' ,'2004') ~ "yod1",
                                                            `Year of diagnosis`%in%c('2005','2006', '2007', '2008' ,'2009') ~ "yod2",
                                                            `Year of diagnosis`%in%c('2010', '2011', '2012' ,'2013' ,'2014') ~ "yod3",
                                                            `Year of diagnosis`%in%c('2015', '2016', '2017' ,'2018' ,'2019') ~ "yod4"))

# make cleaner race/eth column name
new_data2<- new_data2%>%rename(race_eth=`Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`)

# Survival
new_data3<-new_data2%>% mutate(cancer_specific_status = ifelse(`SEER cause-specific death classification`=="Dead (attributable to this cancer dx)",1,0 ))
new_data3<-new_data3%>% mutate(overall_status = ifelse(`Vital status recode (study cutoff used)`=="Alive",0,1 ))
new_data4<-new_data3%>%filter(`Survival months flag`=='Complete dates are available and there are more than 0 days of survival')

# Remove old data sets to save space
rm(seer_data)
rm(new_data)
rm(new_data2)
rm(new_data3)
```

<br>

# 3 Setting reference groups

Makes this easier for filling out tables and for regressions

``` r
# first select only variables we know we will be using to save space
new_data5<-new_data4 %>% dplyr::select(Sex,
                                       race_eth,
                                       can_see_sigmoid,
                                       age_group_final,
                                       new_stage, 
                                       final_grade,
                                       year_of_diag,
                                       cancer_specific_status, 
                                       overall_status,
                                       `Survival months`,
                                       `Primary Site`)%>%
                          rename(survival_months=`Survival months`, primary_site = `Primary Site`)

# set new factor levels, first level is the reference group
new_data5$new_stage <- factor(new_data5$new_stage, levels = c('Localized','Regional','Distant'))
new_data5$age_group_final <- factor(new_data5$age_group_final, levels = c('over50','under45','45-50'))
new_data5$year_of_diag <- factor(new_data5$year_of_diag, levels = c('yod1','yod2','yod3','yod4'))

# we dont need to specify all levels for race that would take too much typing, we just want to set the largest group as the reference
new_data5 <- within(new_data5, race_eth <- relevel(factor(race_eth), ref = 'Non-Hispanic White'))

# The grade variable appears to be unreliable
table(new_data5$final_grade,new_data5$year_of_diag)
```

    ##            
    ##              yod1  yod2  yod3  yod4
    ##   1          4588  4453  4322  6077
    ##   2         45855 47672 48094 50025
    ##   3         12570 12493 10672  9581
    ##   4           426   926  1746  1152
    ##   unk_other  3544  4425  5713 18235

``` r
# now using new_data5
rm(new_data4)

# Print out table 1
table(new_data5$primary_site,new_data5$age_group_final,new_data5$new_stage)
```

    ## , ,  = Localized
    ## 
    ##      
    ##       over50 under45 45-50
    ##   180  13198     389   372
    ##   182  13650     378   367
    ##   183   3115      99    82
    ##   184   5866     256   195
    ##   185   1812     113    81
    ##   186   3521     231   202
    ##   187  17389     987  1048
    ##   199   6519     414   416
    ##   209  18180    1305  1441
    ## 
    ## , ,  = Regional
    ## 
    ##      
    ##       over50 under45 45-50
    ##   180  18453     769   714
    ##   182  15959     778   610
    ##   183   4168     271   166
    ##   184   7985     510   397
    ##   185   3233     264   202
    ##   186   5212     531   405
    ##   187  23032    2213  2041
    ##   199   9707    1012   956
    ##   209  22956    2819  2539
    ## 
    ## , ,  = Distant
    ## 
    ##      
    ##       over50 under45 45-50
    ##   180  10628     507   579
    ##   182   7312     405   421
    ##   183   2149     173   138
    ##   184   3715     321   238
    ##   185   1579     152   138
    ##   186   2641     343   286
    ##   187  14875    1799  1551
    ##   199   6377     736   666
    ##   209  12285    1575  1452

<br>

# 4 Logistic Regression

Calculates odds of tumor in location seen on sigmoidoscopy, first
univariate then multivariate

``` r
# Code below calculates, exps, and does CI, then only saves the output to save space

sex_log<-tidy(glm(can_see_sigmoid~Sex,data=new_data5,family="binomial"),conf.int=TRUE)%>%
          mutate(estimate=exp(estimate), conf.low=exp(conf.low), conf.high=exp(conf.high))

age_log<-tidy(glm(can_see_sigmoid~age_group_final,data=new_data5,family="binomial"),conf.int=TRUE)%>%
          mutate(estimate=exp(estimate), conf.low=exp(conf.low), conf.high=exp(conf.high))

stage_log<-tidy(glm(can_see_sigmoid~new_stage,data=new_data5,family="binomial"),conf.int=TRUE)%>%
            mutate(estimate=exp(estimate), conf.low=exp(conf.low), conf.high=exp(conf.high))

year_log<-tidy(glm(can_see_sigmoid~year_of_diag,data=new_data5,family="binomial"),conf.int=TRUE)%>%
            mutate(estimate=exp(estimate), conf.low=exp(conf.low), conf.high=exp(conf.high))

race_log<-tidy(glm(can_see_sigmoid~race_eth,data=new_data5,family="binomial"),conf.int=TRUE)%>%
            mutate(estimate=exp(estimate), conf.low=exp(conf.low), conf.high=exp(conf.high))

# Multivariate, Check out interactions
multi_log<-tidy(glm(can_see_sigmoid~race_eth+Sex+age_group_final+new_stage+year_of_diag,data=new_data5,family="binomial"),conf.int=TRUE)%>%
  mutate(estimate=exp(estimate),conf.low=exp(conf.low),conf.high=exp(conf.high))
```

<br>

## 4.1 Logistic results: sex

``` r
knitr::kable(sex_log,digits=3)
```

| term        | estimate | std.error | statistic | p.value | conf.low | conf.high |
|:------------|---------:|----------:|----------:|--------:|---------:|----------:|
| (Intercept) |    1.108 |     0.005 |    19.114 |       0 |    1.096 |     1.119 |
| SexMale     |    1.537 |     0.008 |    57.053 |       0 |    1.514 |     1.560 |

## 4.2 Logistic results: age

``` r
knitr::kable(age_log,digits=3)
```

| term                   | estimate | std.error | statistic | p.value | conf.low | conf.high |
|:-----------------------|---------:|----------:|----------:|--------:|---------:|----------:|
| (Intercept)            |    1.265 |     0.004 |    58.960 |       0 |    1.255 |     1.275 |
| age_group_finalunder45 |    2.050 |     0.017 |    43.444 |       0 |    1.985 |     2.118 |
| age_group_final45-50   |    2.187 |     0.017 |    44.779 |       0 |    2.114 |     2.264 |

## 4.3 Logistic results: stage

``` r
knitr::kable(stage_log,digits=3)
```

| term              | estimate | std.error | statistic | p.value | conf.low | conf.high |
|:------------------|---------:|----------:|----------:|--------:|---------:|----------:|
| (Intercept)       |    1.292 |     0.007 |    38.481 |       0 |    1.275 |     1.309 |
| new_stageRegional |    1.043 |     0.009 |     4.816 |       0 |    1.025 |     1.061 |
| new_stageDistant  |    1.213 |     0.010 |    19.090 |       0 |    1.189 |     1.237 |

## 4.4 Logistic results: year

``` r
knitr::kable(year_log,digits=3)
```

| term             | estimate | std.error | statistic | p.value | conf.low | conf.high |
|:-----------------|---------:|----------:|----------:|--------:|---------:|----------:|
| (Intercept)      |    1.339 |     0.008 |    37.335 |   0.000 |    1.318 |     1.359 |
| year_of_diagyod2 |    1.000 |     0.011 |     0.026 |   0.979 |    0.979 |     1.022 |
| year_of_diagyod3 |    1.030 |     0.011 |     2.717 |   0.007 |    1.008 |     1.052 |
| year_of_diagyod4 |    1.085 |     0.010 |     7.774 |   0.000 |    1.063 |     1.107 |

## 4.5 Logistic results: race

``` r
knitr::kable(race_log,digits=3)
```

| term                                               | estimate | std.error | statistic | p.value | conf.low | conf.high |
|:---------------------------------------------------|---------:|----------:|----------:|--------:|---------:|----------:|
| (Intercept)                                        |    1.320 |     0.005 |    60.889 |       0 |    1.308 |     1.332 |
| race_ethHispanic (All Races)                       |    1.269 |     0.012 |    19.679 |       0 |    1.239 |     1.299 |
| race_ethNon-Hispanic American Indian/Alaska Native |    1.228 |     0.044 |     4.642 |       0 |    1.126 |     1.339 |
| race_ethNon-Hispanic Asian or Pacific Islander     |    1.596 |     0.014 |    33.455 |       0 |    1.553 |     1.640 |
| race_ethNon-Hispanic Black                         |    0.803 |     0.012 |   -18.407 |       0 |    0.785 |     0.822 |
| race_ethNon-Hispanic Unknown Race                  |    1.388 |     0.068 |     4.820 |       0 |    1.216 |     1.588 |

## 4.6 Logistic results: multivariate

``` r
knitr::kable(multi_log,digits=3)
```

| term                                               | estimate | std.error | statistic | p.value | conf.low | conf.high |
|:---------------------------------------------------|---------:|----------:|----------:|--------:|---------:|----------:|
| (Intercept)                                        |    0.956 |     0.011 |    -4.279 |   0.000 |    0.936 |     0.976 |
| race_ethHispanic (All Races)                       |    1.182 |     0.012 |    13.566 |   0.000 |    1.154 |     1.211 |
| race_ethNon-Hispanic American Indian/Alaska Native |    1.174 |     0.045 |     3.579 |   0.000 |    1.076 |     1.282 |
| race_ethNon-Hispanic Asian or Pacific Islander     |    1.561 |     0.014 |    31.496 |   0.000 |    1.519 |     1.605 |
| race_ethNon-Hispanic Black                         |    0.773 |     0.012 |   -21.269 |   0.000 |    0.755 |     0.792 |
| race_ethNon-Hispanic Unknown Race                  |    1.345 |     0.069 |     4.295 |   0.000 |    1.176 |     1.541 |
| SexMale                                            |    1.523 |     0.008 |    55.264 |   0.000 |    1.501 |     1.546 |
| age_group_finalunder45                             |    2.013 |     0.017 |    41.810 |   0.000 |    1.948 |     2.080 |
| age_group_final45-50                               |    2.146 |     0.018 |    43.247 |   0.000 |    2.073 |     2.222 |
| new_stageRegional                                  |    1.007 |     0.009 |     0.753 |   0.452 |    0.989 |     1.024 |
| new_stageDistant                                   |    1.161 |     0.010 |    14.497 |   0.000 |    1.138 |     1.185 |
| year_of_diagyod2                                   |    0.977 |     0.011 |    -2.057 |   0.040 |    0.956 |     0.999 |
| year_of_diagyod3                                   |    0.986 |     0.011 |    -1.309 |   0.190 |    0.964 |     1.007 |
| year_of_diagyod4                                   |    1.022 |     0.011 |     2.015 |   0.044 |    1.001 |     1.043 |

<br>

# 5 Survival exploration w/CSS

Makes KM curves using CSS, checks proportional hazards assumptions

## 5.1 KM plots: Can see on sigmoidoscopy

``` r
# CAN SEE
km_cansee<-survfit(Surv(survival_months, cancer_specific_status) ~ can_see_sigmoid, data = new_data5)
plot(km_cansee,fun="cloglog")
ggsurvplot(km_cansee,censor=FALSE) # BIGGEST TIME DEPENDANT!!!
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

## 5.2 KM plots: Stage

``` r
# STAGE
km_new_stage<-survfit(Surv(survival_months, cancer_specific_status) ~new_stage, data = new_data5)
plot(km_new_stage,fun="cloglog")
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## 5.3 KM plots: year of diag

``` r
# YEAR OF DIAG
km_year_of_diag<-survfit(Surv(survival_months, cancer_specific_status) ~ year_of_diag, data = new_data5)
ggsurvplot(km_year_of_diag,censor=FALSE)
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
plot(km_year_of_diag,fun="cloglog")
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

## 5.4 KM plots: sex

``` r
# SEX
km_Sex<-survfit(Surv(survival_months, cancer_specific_status) ~ Sex, data = new_data5)
ggsurvplot(km_Sex,censor=FALSE)
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
plot(km_Sex,fun="cloglog") # lines cross but tiny
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

## 5.5 KM plots: age group

``` r
# AGE GROUP
km_age_group_final<-survfit(Surv(survival_months, cancer_specific_status) ~ age_group_final, data = new_data5)
ggsurvplot(km_age_group_final,censor=FALSE) # lines cross but tiny
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
plot(km_age_group_final,fun="cloglog") 
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

## 5.6 KM plots: race/eth

``` r
# RACE ETH
km_race_eth<-survfit(Surv(survival_months, cancer_specific_status) ~ race_eth, data = new_data5)
ggsurvplot(km_race_eth,censor=FALSE) # lines cross but insignificant 
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
plot(km_race_eth,fun="cloglog") 
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

## 5.7 KM plots: Site

``` r
# SITE (this is just can see sigmoid with more definition)
km_site<-survfit(Surv(survival_months, cancer_specific_status) ~ primary_site, data = new_data5)
ggsurvplot(km_site,censor=FALSE,ylim=c(0.5,1),xlim=c(0,100))
```

    ## Warning: Removed 320 rows containing missing values (`geom_step()`).
    ## Removed 320 rows containing missing values (`geom_step()`).

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
#table(new_data5$primary_site)
```

<br>

# 6 Survival Analysis

``` r
# take out patients with 0 survival months
new_data6<-new_data5%>%filter(survival_months>0)
new_data7<- new_data6%>%select(survival_months,cancer_specific_status,can_see_sigmoid,age_group_final,new_stage,overall_status)

# computes rmst with age group as a covariate, takes dataframe, vector of time points, and the stage as inputs, returns a dataframe with results across time points in columns, rows are the stats we are looking at (add option to switch to OS)
compute_rmst_with_age <- function(data,time_points,stage,css=TRUE){
  
  data<-data %>% filter(new_stage==stage)
  covs<-model.matrix(~age_group_final,data=data)[,2:3]
  
  rms_combined <- data.frame(col1="rmst")
  
  for (t in time_points){
    # calcs rmst with alpha corrected for 3 stages and number of time points 
    if (css) {
       rms_temp <- rmst2(time=data$survival_months,status=data$cancer_specific_status,arm=data$can_see_sigmoid,covariates=covs,tau=t,alpha=0.05/(3*length(time_points)) )
    }else{
       rms_temp <- rmst2(time=data$survival_months,status=data$cancer_specific_status,arm=data$can_see_sigmoid,covariates=covs,tau=t,alpha=0.05/(3*length(time_points)) )
    }
   
    # round and add results to the data frame, each time point will be added as a block of columns 
    rms_temp_res <- data.frame(rms_temp$RMST.difference.adjusted) %>% select(-se.coef., -z ) %>% mutate_all((function(x) round(x,digits=3))) %>% setNames( paste0('m_',t,"_",names(.)) )
    rms_combined <- cbind(rms_combined,rms_temp_res)
  }
  
  return(rms_combined)
}


local_res_combined_css <- compute_rmst_with_age(new_data7,c(12,24,60,120),"Localized")
regional_res_combined_css <- compute_rmst_with_age(new_data7,c(12,24,60,120),"Regional")
distant_res_combined_css <- compute_rmst_with_age(new_data7,c(12,24,60,120),"Distant")

local_res_combined_os <- compute_rmst_with_age(new_data7,c(12,24,60,120),"Localized",FALSE)
regional_res_combined_os <- compute_rmst_with_age(new_data7,c(12,24,60,120),"Regional",FALSE)
distant_res_combined_os <- compute_rmst_with_age(new_data7,c(12,24,60,120),"Distant",FALSE)

knitr::kable(local_res_combined_css)
```

|                        | col1 | m_12_coef | m_12_p | m_12_lower..100 | m_12_upper..100 | m_24_coef | m_24_p | m_24_lower..100 | m_24_upper..100 | m_60_coef | m_60_p | m_60_lower..100 | m_60_upper..100 | m_120_coef | m_120_p | m_120_lower..100 | m_120_upper..100 |
|------------------------|:-----|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|-----------:|--------:|-----------------:|-----------------:|
| intercept              | rmst |    11.720 |  0.000 |          11.697 |          11.743 |    23.117 |      0 |          23.059 |          23.175 |    55.900 |      0 |          55.699 |          56.100 |    107.335 |       0 |          106.811 |          107.859 |
| arm                    | rmst |    -0.017 |  0.105 |          -0.047 |           0.013 |    -0.101 |      0 |          -0.178 |          -0.025 |    -0.989 |      0 |          -1.261 |          -0.718 |     -4.199 |       0 |           -4.916 |           -3.482 |
| age_group_finalunder45 | rmst |     0.235 |  0.000 |           0.199 |           0.271 |     0.727 |      0 |           0.629 |           0.825 |     2.961 |      0 |           2.528 |           3.395 |      8.941 |       0 |            7.510 |           10.372 |
| age_group_final45-50   | rmst |     0.240 |  0.000 |           0.206 |           0.274 |     0.760 |      0 |           0.666 |           0.855 |     3.191 |      0 |           2.773 |           3.609 |      9.040 |       0 |            7.630 |           10.450 |

``` r
knitr::kable(regional_res_combined_css)
```

|                        | col1 | m_12_coef | m_12_p | m_12_lower..100 | m_12_upper..100 | m_24_coef | m_24_p | m_24_lower..100 | m_24_upper..100 | m_60_coef | m_60_p | m_60_lower..100 | m_60_upper..100 | m_120_coef | m_120_p | m_120_lower..100 | m_120_upper..100 |
|------------------------|:-----|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|-----------:|--------:|-----------------:|-----------------:|
| intercept              | rmst |    11.361 |      0 |          11.333 |          11.389 |    21.627 |      0 |          21.552 |          21.702 |    47.923 |      0 |          47.659 |          48.187 |     85.222 |       0 |           84.570 |           85.874 |
| arm                    | rmst |     0.215 |      0 |           0.181 |           0.248 |     0.751 |      0 |           0.661 |           0.841 |     2.361 |      0 |           2.036 |           2.687 |      2.434 |       0 |            1.608 |            3.259 |
| age_group_finalunder45 | rmst |     0.364 |      0 |           0.330 |           0.397 |     1.133 |      0 |           1.025 |           1.240 |     4.712 |      0 |           4.215 |           5.209 |     13.722 |       0 |           12.144 |           15.300 |
| age_group_final45-50   | rmst |     0.366 |      0 |           0.330 |           0.401 |     1.141 |      0 |           1.029 |           1.252 |     4.460 |      0 |           3.930 |           4.990 |     12.140 |       0 |           10.427 |           13.852 |

``` r
knitr::kable(distant_res_combined_css)
```

|                        | col1 | m_12_coef | m_12_p | m_12_lower..100 | m_12_upper..100 | m_24_coef | m_24_p | m_24_lower..100 | m_24_upper..100 | m_60_coef | m_60_p | m_60_lower..100 | m_60_upper..100 | m_120_coef | m_120_p | m_120_lower..100 | m_120_upper..100 |
|------------------------|:-----|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|-----------:|--------:|-----------------:|-----------------:|
| intercept              | rmst |     8.594 |      0 |           8.516 |           8.672 |    13.314 |      0 |          13.147 |          13.481 |    19.169 |      0 |          18.803 |          19.535 |     23.495 |       0 |           22.837 |           24.153 |
| arm                    | rmst |     0.887 |      0 |           0.793 |           0.980 |     2.405 |      0 |           2.198 |           2.612 |     5.569 |      0 |           5.094 |           6.043 |      7.450 |       0 |            6.597 |            8.304 |
| age_group_finalunder45 | rmst |     1.461 |      0 |           1.339 |           1.584 |     3.298 |      0 |           2.982 |           3.614 |     6.087 |      0 |           5.154 |           7.019 |     10.042 |       0 |            7.620 |           12.464 |
| age_group_final45-50   | rmst |     1.220 |      0 |           1.083 |           1.357 |     2.853 |      0 |           2.507 |           3.198 |     5.696 |      0 |           4.697 |           6.696 |      8.554 |       0 |            6.074 |           11.033 |

``` r
knitr::kable(local_res_combined_os)
```

|                        | col1 | m_12_coef | m_12_p | m_12_lower..100 | m_12_upper..100 | m_24_coef | m_24_p | m_24_lower..100 | m_24_upper..100 | m_60_coef | m_60_p | m_60_lower..100 | m_60_upper..100 | m_120_coef | m_120_p | m_120_lower..100 | m_120_upper..100 |
|------------------------|:-----|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|-----------:|--------:|-----------------:|-----------------:|
| intercept              | rmst |    11.720 |  0.000 |          11.697 |          11.743 |    23.117 |      0 |          23.059 |          23.175 |    55.900 |      0 |          55.699 |          56.100 |    107.335 |       0 |          106.811 |          107.859 |
| arm                    | rmst |    -0.017 |  0.105 |          -0.047 |           0.013 |    -0.101 |      0 |          -0.178 |          -0.025 |    -0.989 |      0 |          -1.261 |          -0.718 |     -4.199 |       0 |           -4.916 |           -3.482 |
| age_group_finalunder45 | rmst |     0.235 |  0.000 |           0.199 |           0.271 |     0.727 |      0 |           0.629 |           0.825 |     2.961 |      0 |           2.528 |           3.395 |      8.941 |       0 |            7.510 |           10.372 |
| age_group_final45-50   | rmst |     0.240 |  0.000 |           0.206 |           0.274 |     0.760 |      0 |           0.666 |           0.855 |     3.191 |      0 |           2.773 |           3.609 |      9.040 |       0 |            7.630 |           10.450 |

``` r
knitr::kable(regional_res_combined_os)
```

|                        | col1 | m_12_coef | m_12_p | m_12_lower..100 | m_12_upper..100 | m_24_coef | m_24_p | m_24_lower..100 | m_24_upper..100 | m_60_coef | m_60_p | m_60_lower..100 | m_60_upper..100 | m_120_coef | m_120_p | m_120_lower..100 | m_120_upper..100 |
|------------------------|:-----|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|-----------:|--------:|-----------------:|-----------------:|
| intercept              | rmst |    11.361 |      0 |          11.333 |          11.389 |    21.627 |      0 |          21.552 |          21.702 |    47.923 |      0 |          47.659 |          48.187 |     85.222 |       0 |           84.570 |           85.874 |
| arm                    | rmst |     0.215 |      0 |           0.181 |           0.248 |     0.751 |      0 |           0.661 |           0.841 |     2.361 |      0 |           2.036 |           2.687 |      2.434 |       0 |            1.608 |            3.259 |
| age_group_finalunder45 | rmst |     0.364 |      0 |           0.330 |           0.397 |     1.133 |      0 |           1.025 |           1.240 |     4.712 |      0 |           4.215 |           5.209 |     13.722 |       0 |           12.144 |           15.300 |
| age_group_final45-50   | rmst |     0.366 |      0 |           0.330 |           0.401 |     1.141 |      0 |           1.029 |           1.252 |     4.460 |      0 |           3.930 |           4.990 |     12.140 |       0 |           10.427 |           13.852 |

``` r
knitr::kable(distant_res_combined_os)
```

|                        | col1 | m_12_coef | m_12_p | m_12_lower..100 | m_12_upper..100 | m_24_coef | m_24_p | m_24_lower..100 | m_24_upper..100 | m_60_coef | m_60_p | m_60_lower..100 | m_60_upper..100 | m_120_coef | m_120_p | m_120_lower..100 | m_120_upper..100 |
|------------------------|:-----|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|----------:|-------:|----------------:|----------------:|-----------:|--------:|-----------------:|-----------------:|
| intercept              | rmst |     8.594 |      0 |           8.516 |           8.672 |    13.314 |      0 |          13.147 |          13.481 |    19.169 |      0 |          18.803 |          19.535 |     23.495 |       0 |           22.837 |           24.153 |
| arm                    | rmst |     0.887 |      0 |           0.793 |           0.980 |     2.405 |      0 |           2.198 |           2.612 |     5.569 |      0 |           5.094 |           6.043 |      7.450 |       0 |            6.597 |            8.304 |
| age_group_finalunder45 | rmst |     1.461 |      0 |           1.339 |           1.584 |     3.298 |      0 |           2.982 |           3.614 |     6.087 |      0 |           5.154 |           7.019 |     10.042 |       0 |            7.620 |           12.464 |
| age_group_final45-50   | rmst |     1.220 |      0 |           1.083 |           1.357 |     2.853 |      0 |           2.507 |           3.198 |     5.696 |      0 |           4.697 |           6.696 |      8.554 |       0 |            6.074 |           11.033 |

<br><br>

# 7 KM Plots with better 95% CIs

``` r
# makes a KM fit and calcs logep 95% CIs, returns the data for use with ggplot
make_plot_dat <- function(data_set, age_grp, stage, css = TRUE){
  
  can <- data_set%>%filter(age_group_final==age_grp,new_stage==stage,can_see_sigmoid==1)
  nocan <- data_set%>%filter(age_group_final==age_grp,new_stage==stage,can_see_sigmoid==0)
  
  if (css){
    fit1 <- surv_fit(Surv(survival_months,cancer_specific_status)~1,data=can)
    fit2 <-  surv_fit(Surv(survival_months,cancer_specific_status)~1,data=nocan)
  }else{
    fit1 <- surv_fit(Surv(survival_months,overall_status)~1,data=can)
    fit2 <- surv_fit(Surv(survival_months,overall_status)~1,data=nocan)
  }
  
  fit1 <- km.ci(fit1,conf.level = 0.95,method = "logep")
  fit2 <- km.ci(fit2,conf.level = 0.95,method = "logep")
  
  pltdat1 <- as.data.frame(cbind(fit1$surv,fit1$time,fit1$n.risk,fit1$n.event,fit1$lower,fit1$upper))%>%mutate(type="cansee")
  pltdat2 <- as.data.frame(cbind(fit2$surv,fit2$time,fit2$n.risk,fit2$n.event,fit2$lower,fit2$upper))%>%mutate(type="cant see")
  
  colnames(pltdat1)<-c("surv","time","nrisk","nevent","lower","upper","type")
  colnames(pltdat2)<-c("surv","time","nrisk","nevent","lower","upper","type")

  plt3 <- rbind(pltdat1,pltdat2)
  return(plt3)
}

make_plot <- function(plot_dat,fig_title){# need to add in adjustment based on stage
  
  p <- ggplot(plot_dat ,aes(x=time,y=surv,color=type))+
    geom_step()+
    ylim(if(str_detect(fig_title,"Distant")) c(0,1) else c(0,1))+
    geom_ribbon(aes(ymin=lower,ymax=upper,x=time,fill=type),alpha=.1,size=0.2)+
    ggtitle(fig_title)+theme_test()+
    theme(
      plot.title = element_text(size=8),
      axis.title = element_blank()
      )+
    scale_x_continuous(limits=c(0,240),breaks=c(0,24,60,120,180),labels=c("0","2","5","10","15"))
  return(p)
}
```

<br><br>

# 8 CSS Plots

``` r
# Local
u45_local_plot <- make_plot_dat(new_data7, "under45", "Localized")
a<-make_plot(u45_local_plot,"A: <45 Local")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

``` r
midage_local_plot <- make_plot_dat(new_data7, "45-50", "Localized")
b<-make_plot(midage_local_plot,"B: 45-50 Local")

o50_local_plot <- make_plot_dat(new_data7, "over50", "Localized")
c<-make_plot(o50_local_plot ,"C: >50 Local")


# Regional
u45_reg_plot <- make_plot_dat(new_data7, "under45", "Regional")
d<-make_plot(u45_reg_plot,"D: <45 Regional")

midage_reg_plot <- make_plot_dat(new_data7, "45-50", "Regional")
e<-make_plot(midage_reg_plot,"E: 45-50 Regional")

o50_reg_plot <- make_plot_dat(new_data7, "over50", "Regional")
f<-make_plot(o50_reg_plot,"F: >50 Regional")


# Distant
u45_dist_plot <- make_plot_dat(new_data7, "under45", "Distant")
g<-make_plot(u45_dist_plot,"G:<45 Distant")

midage_dist_plot <- make_plot_dat(new_data7, "45-50", "Distant")
h<-make_plot(midage_dist_plot ,"H: 45-50 Distant")

o50_dist_plot <- make_plot_dat(new_data7, "over50", "Distant")
i<-make_plot(o50_dist_plot ,"I: >50 Distant")

ggarrange(a,b,c,d,e,f,g,h,i,common.legend = TRUE)
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
ggsave(plot=last_plot(),file="css_km.png",height=8, width=11.5)
```

<br><br>

# 9 OS plots

``` r
# Local
u45_local_plot <- make_plot_dat(new_data7, "under45", "Localized",FALSE)
a<-make_plot(u45_local_plot,"A: <45 Local")

midage_local_plot <- make_plot_dat(new_data7, "45-50", "Localized",FALSE)
b<-make_plot(midage_local_plot,"B: 45-50 Local")

o50_local_plot <- make_plot_dat(new_data7, "over50", "Localized",FALSE)
c<-make_plot(o50_local_plot ,"C: >50 Local")


# Regional
u45_reg_plot <- make_plot_dat(new_data7, "under45", "Regional",FALSE)
d<-make_plot(u45_reg_plot,"D: <45 Regional")

midage_reg_plot <- make_plot_dat(new_data7, "45-50", "Regional",FALSE)
e<-make_plot(midage_reg_plot,"E: 45-50 Regional")

o50_reg_plot <- make_plot_dat(new_data7, "over50", "Regional",FALSE)
f<-make_plot(o50_reg_plot,"F: >50 Regional")


# Distant
u45_dist_plot <- make_plot_dat(new_data7, "under45", "Distant",FALSE)
g<-make_plot(u45_dist_plot,"G:<45 Distant")

midage_dist_plot <- make_plot_dat(new_data7, "45-50", "Distant",FALSE)
h<-make_plot(midage_dist_plot ,"H: 45-50 Distant")

o50_dist_plot <- make_plot_dat(new_data7, "over50", "Distant",FALSE)
i<-make_plot(o50_dist_plot ,"I: >50 Distant")

ggarrange(a,b,c,d,e,f,g,h,i,common.legend = TRUE)+ggtitle("OS")
```

![](SEER_CRC_Sigmoidoscopy_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
ggsave(plot=last_plot(),file="os_km.png",height=8, width=11.5)
```

<br><br>

# 10 Session Info

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] km.ci_0.5-6     survRM2_1.0-4   forcats_1.0.0   stringr_1.5.0  
    ##  [5] dplyr_1.1.0     purrr_1.0.1     readr_2.1.3     tidyr_1.3.0    
    ##  [9] tibble_3.1.8    tidyverse_1.3.2 survminer_0.4.9 ggpubr_0.5.0   
    ## [13] ggplot2_3.4.0   survival_3.4-0  MASS_7.3-58.1   broom_1.0.3    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fs_1.6.0            lubridate_1.9.1     bit64_4.0.5        
    ##  [4] httr_1.4.4          tools_4.2.2         backports_1.4.1    
    ##  [7] utf8_1.2.2          R6_2.5.1            DBI_1.1.3          
    ## [10] colorspace_2.1-0    withr_2.5.0         tidyselect_1.2.0   
    ## [13] gridExtra_2.3       bit_4.0.5           compiler_4.2.2     
    ## [16] textshaping_0.3.6   cli_3.4.1           rvest_1.0.3        
    ## [19] xml2_1.3.3          labeling_0.4.2      scales_1.2.1       
    ## [22] survMisc_0.5.6      systemfonts_1.0.4   digest_0.6.31      
    ## [25] rmarkdown_2.20      pkgconfig_2.0.3     htmltools_0.5.4    
    ## [28] dbplyr_2.3.0        fastmap_1.1.0       highr_0.10         
    ## [31] rlang_1.0.6         readxl_1.4.1        rstudioapi_0.14    
    ## [34] generics_0.1.3      farver_2.1.1        zoo_1.8-11         
    ## [37] jsonlite_1.8.4      vroom_1.6.1         car_3.1-1          
    ## [40] googlesheets4_1.0.1 magrittr_2.0.3      Matrix_1.5-1       
    ## [43] munsell_0.5.0       fansi_1.0.3         abind_1.4-5        
    ## [46] lifecycle_1.0.3     stringi_1.7.12      yaml_2.3.7         
    ## [49] carData_3.0-5       grid_4.2.2          parallel_4.2.2     
    ## [52] crayon_1.5.2        lattice_0.20-45     haven_2.5.1        
    ## [55] cowplot_1.1.1       splines_4.2.2       hms_1.1.2          
    ## [58] knitr_1.42          pillar_1.8.1        ggsignif_0.6.4     
    ## [61] reprex_2.0.2        glue_1.6.2          evaluate_0.20      
    ## [64] data.table_1.14.6   modelr_0.1.10       vctrs_0.5.2        
    ## [67] tzdb_0.3.0          cellranger_1.1.0    gtable_0.3.1       
    ## [70] assertthat_0.2.1    xfun_0.37           xtable_1.8-4       
    ## [73] rstatix_0.7.2       ragg_1.2.5          googledrive_2.0.0  
    ## [76] gargle_1.3.0        KMsurv_0.1-5        timechange_0.2.0   
    ## [79] ellipsis_0.3.2
