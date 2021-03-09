****************************************************************************************
************* Step 0: Data Preproceesing and Merge with Census Data ********************
****************************************************************************************
cd "/Volumes/Disk E/Data Plus/data/raw/From REA/Append Data"

use "REA_ASR_3846_final_new", clear
sort state sys_no year
* We propose two ways to measure kwh consumed per residential customer
* 1) For 1938-1940, use (original indicator in report * 12); For 1941-1946, use kwh_billed/consumer connected.
* 2) original indicator in report * 12 (op_kwh_pr*12)
* Note that in the report, kwh_billed is a "flow" variable, while consumer connected is a commulative variable
gen kwh_per_cc1= op_kwh_b/op_ccn
replace kwh_per_cc1= op_kwh_pr*12 if year<1941
replace kwh_per_cc1= op_kwh_pr*12 if year>=1941 & kwh_per_cc1==.
gen kwh_per_cc2= op_kwh_pr*12
label variable kwh_per_cc1 "kwh_bill/number of customers (41-46) & kwh per customer*12 (38-40) "
label variable kwh_per_cc2 "kwh per customer*12 (38-43)"
* since there are some outliers, trim upper way a little for the first indicator at the 1% threshold
winsor2 kwh_per_cc1, replace cut(1 99) trim
winsor2 kwh_per_cc2, replace cut(1 99) trim
* generate log-scale indicators
gen lnstop_dist= ln(min_stop_dist)
gen lnefes_dist= ln(min_efes_dist)
gen ln_kwh1= ln(kwh_per_cc1)
gen ln_kwh2= ln(kwh_per_cc2)
gen ln_mile= ln(op_me)
* generate code indicating state, county or city
encode county, gen(county_code)
encode city, gen(city_code)
encode state, gen(state_code)
encode system, gen(system_code)

/* merge with agriculture, demographic and political dataset to obtain macro-level control variables
gen id_master= _n
reclink state county using "/Volumes/Disk E/Data Plus/data/raw/ICPRS/demo_agri_3045", gen(simiscore) idm(id_master) idu(id_using)
drop simiscore id_using Ustate Ucounty _m
merge m:1 state using "/Volumes/Disk E/Data Plus/data/raw/ICPRS/Political Data/Strength Data_wideformat"
drop if _m==2
drop id_master id_using _m
*/

** interpolate the census data, regarded as constant in 5-year period
* available in 35, 40, 45 data
foreach var in farmsize farm_pop{
gen `var'= .
replace `var'= `var'_35 if year < 1940
replace `var'= `var'_40 if year >= 1940 & year <1945
replace `var'= `var'_45 if year >= 1945
}
* available in 35, 40 data
gen percent_white= .
replace percent_white= per_white_35 if year < 1940
replace percent_white= per_white_40 if year >= 1940
* available in 40, 45 data
foreach i in 39 44{
gen ratio_croptolivestock_`i'= val_crop_`i'/val_livestock_`i'
gen ratio_croptotal_`i'= val_crop_`i'/(val_crop_`i' + val_livestock_`i')
} 
foreach var in ratio_croptolivestock ratio_croptotal{
gen `var'= .
replace `var'= `var'_39 if year < 1944
replace `var'= `var'_44 if year >=1944
}
foreach var in electri_distri_line{
gen `var'= .
replace `var'= `var'_40 if year < 1945
replace `var'= `var'_45 if year >=1945
}
* political data (available in even year and impute in odd year)
gen RCOMPB= .
gen DCOMPB= .
foreach var in RCOMPB DCOMPB{
replace `var'= `var'1938 if year>=1938 & year<=1939
replace `var'= `var'1940 if year>=1940 & year<=1941
replace `var'= `var'1942 if year>=1942 & year<=1943
replace `var'= `var'1944 if year>=1944 & year<=1945
replace `var'= `var'1946 if year==1946
}
* take log transformation of some covariates
foreach var in farmsize_35 farmsize_40 farmsize_45 farm_pop_35 farm_pop_40 farm_pop_45 electri_distri_line_40 electri_distri_line_45 hard_road_40{
gen ln_`var'= ln(`var')
}
** generate the lag term of covariates
* available in 35, 40, 45 data
foreach var in farmsize farm_pop{
gen ln_`var'_lag= .
replace ln_`var'_lag= ln_`var'_35 if year < 1941
replace ln_`var'_lag= ln_`var'_40 if year >= 1941 & year <1946
replace ln_`var'_lag= ln_`var'_45 if year >= 1946
}
* available in 35, 40 data
gen percent_white_lag= .
replace percent_white_lag= per_white_35 if year < 1941
replace percent_white_lag= per_white_40 if year >= 1941
* available in 40, 45 data
foreach var in ratio_croptolivestock ratio_croptotal{
gen `var'_lag= .
replace `var'_lag= `var'_39 if year < 1945
replace `var'_lag= `var'_44 if year >=1945
}
foreach var in ln_electri_distri_line{
gen `var'_lag= .
replace `var'_lag= `var'_40 if year < 1946
replace `var'_lag= `var'_45 if year >= 1946
}
* political data
gen RCOMPB_lag= .
gen DCOMPB_lag= .
foreach var in RCOMPB DCOMPB{
replace `var'_lag= `var'1936 if year==1938
replace `var'_lag= `var'1938 if year>=1939 & year<=1940
replace `var'_lag= `var'1940 if year>=1941 & year<=1942
replace `var'_lag= `var'1942 if year>=1943 & year<=1944
replace `var'_lag= `var'1944 if year>=1945 & year<=1946
}
* compute higher order terms
foreach v2 of varlist dist_tolargecity ln_farmsize_lag ln_farm_pop_lag percent_white_lag ratio_croptolivestock_lag ratio_croptotal_lag ln_electri_distri_line_lag ln_hard_road_40{
	egen tmp=mean(`v2')
	g `v2'_c=`v2'-tmp
	drop tmp
	g `v2'_sq=`v2'_c^2
	drop `v2'_c
	}
*
drop *_29 *_30 *_44 *_45 *46 *48 *36 *38 *40 *42
duplicates drop system_code year, force
sort state sys_no year
compress
save "/Volumes/Disk E/Data Plus/data/PSM/Data/pre_data", replace


******* Note: the following code only includes the distance to Farm Equipment Tour stops

****************************************************************************
************* Step 1: Generate Treat and START Variable ********************
****************************************************************************
cd "/Volumes/Disk E/Data Plus/data/PSM"

// Threshold of treatment group (<=30 miles)
use "Data/pre_data", clear
* 
xtset system_code year, yearly
gen after= 0
replace after= 1 if min_efes_dist <= 30
label variable after "=1, if site i is located within 30 miles of roadshow stop in year t"
* some treatment sites are located within 30 miles of roadshow stop from the start of sample periods (i.e. year 1938), we generate a label for these sites
bys system: egen min_after= min(after)
gen sys_always_treat= 1 if min_after== 1
drop min_after
label variable sys_always_treat "=1, if treatment sites are located within 30 miles of roadshow stop from the start of sample periods"
* drop the sites that only have 1-year observation (since these sites are "young", the data is not realible)
bys system: gen obs_nosite= _N
drop if obs_nosite==1
drop obs_nosite
* drop sites that have gaps on their data (year variable is not consecutive)
xtset system_code year
bysort system_code (year) : drop if _N < (year[_N] - year[1] + 1) 
* drop sites that don't have outcome variable (ln_kwh1) during the whole sample periods
bys system_code: egen maxkwh= max(ln_kwh1)
drop if  maxkwh==.
drop maxkwh
* generate leading variable (outcome variable: p=0, 1, 2, 3, 4, 5)
xtset system_code year
foreach var of varlist kwh_per_cc1 ln_kwh1{
gen l`var'1= l1.`var'
gen f`var'1= f1.`var'
gen f`var'2= f2.`var'
gen f`var'3= f3.`var'
gen f`var'4= f4.`var'
gen f`var'5= f5.`var'
}
* generate the lag of energized miles
bys system: gen lag_mile= ln_mile[_n-1]
* generate a variable which equals to 0 if the site is never located within 30 miles radius of a roadshow stop
gen T= 0
bys system: egen max_after= max(after)
replace T= 1 if max_after==1
drop max_after
label variable T "=0, site is never located within 30 miles radius of a roadshow stop during sample periods"
* generate start variable (START=1 only if site i starts being treated in that year)
xtset system_code year, yearly
sort system_code year
gen after_lag= l.after
gen after_diff= after-after_lag
gen START= 0
replace START= 1 if after_diff== 1 
* !! note: some sites are treated since the start year of sample, whether we should drop these sites?
replace START= 1 if after_lag==. & after==1
drop after_lag after_diff
// theoretically, we should drop the sites that are treated since entering the data
drop if sys_always_treat==1
* check if each site has only one start year
bys system_code: egen temp= total(START)
tab temp
* for treatment site, only keep the records before (including) they start being treated
gen temp2= year if START==1
bys system_code: egen temp3= min(temp2)
replace START= . if year>temp3
drop if START==.
drop temp*
label variable START "=1, site i starts being treated in that year"
* since stop information is not updated since 1941, we drop the records after this year.
drop if year>=1942
save "Data/PSP_Sample_BeforeMatched_t0", replace

*****************************************************************************************************
************* Step 2: Use Logit Model to Estimate Propensity Score for each site ********************
*****************************************************************************************************

* this is the control variables used in propensity score matching and later regression analysis

global syvar "ln_mile"
global svar "dist_tolargecity dist_tolargecity_sq"
global cvar "ratio_croptotal_lag ln_electri_distri_line_lag"
global cyvar "ln_farmsize_lag  percent_white_lag"
global pvar "DCOMPB_lag"

set more off
set matsize 11000

* use logit model to estimate propensity score
use "Data/PSP_Sample_BeforeMatched_t0", clear
logit START $syvar $svar $cvar $cyvar $pvar i.year
est store logit_model
esttab logit_model using "Logit Estimation.rtf" ,replace se pr2 nogap ///
   keep($syvar $svar $cvar $cyvar $pvar) ///
   scalars(N) sfmt(%9.0g %9.3g %9.0g) star(* 0.1 ** 0.05 *** 0.01) b(%6.4f) 
predict phat if e(sample), pr 
drop if phat==.

*****************************************************************************************************
************* Step 3: Match using kernel method year by year ****************************************
*****************************************************************************************************

* rescale estimated propensity score to ensure same-"start"-year control unit is matched with treated unit
egen idyr= group(year)
gen double phat_exact=phat+(idyr-1)*10 /* (idyr-1) part can reflect the difference of year */

* matching
psmatch2 START, kernel pscore(phat_exact) com outcome(ln_kwh1)
drop if _support==0 /* drop treated unit out of support */
drop if _w==. /* drop unmatched sample */
g esample = e(sample) /* e(sample) equals to 1 if the observation is in the estimation sample and 0 otherwise, thus we can know whether observation is in the matching sample */
save "Data/PSP_Sample_Matched_t0", replace

************************************************************************************************
*************************Step 4: Conduct Balance Property Checking******************************
************************************************************************************************

* store mean bias and chi_square joint significance test
putexcel set "pstest_indicator", sh("PSIND_t0") replace
pstest $syvar $svar $cvar $cyvar $pvar
putexcel A1 = "Test Indicator"
putexcel A2 = "Median Bias"
putexcel A3 = "Prob(chi2)"
putexcel B1 = "Value"
putexcel B2 = `r(meanbias)'
putexcel B3 = `r(chiprob)'
putexcel C1 = "Period"
putexcel C2 = "s=0"
putexcel C3 = "s=0"
putexcel close
* store t-value test result (check if there is no systermatic difference of covariates after matching)
putexcel set "pstest_indicator", sh("pst_t0_balance") modify
putexcel A1 = "Variable"
putexcel B1 = "Difference"
putexcel C1 = "T-value"
putexcel D1 = "p-value"
local j=2
/* there are total 8 control variables */
while `j'<=9{
foreach var in $syvar $svar $cvar $cyvar $pvar{
svyset _n [pw=_w]
putexcel A`j' = "`var'"
svy: mean `var', over(START) coeflegend
lincom _b[`var':1]-_b[`var':0] 
putexcel B`j' = `r(estimate)'
putexcel C`j' = `r(t)'
putexcel D`j' = `r(p)'
local j=`j'+1
}
}
putexcel close
* store number of treated and caliper
g treated=1 if START==1
g N=_N
tabstat treated N, s(N mean) c(s) save
tabstatmat all
matrix Stats=all'
xml_tab Stats , save("Baseline_Treated.xml") replace sheet(obsT_t0) 
* store mean of treated and control group
svyset _n [pw=_w]
svy: mean $syvar $svar $cvar $cyvar $pvar if  START==1
est store treated
svy: mean $syvar $svar $cvar $cyvar $pvar if  START==0
est store control
xml_tab treated control , save("Baseline_Treated.xml") replace sheet(PST_t0) 
est clear

************************************************************************************************
****************************** Step 5: Loop over p=-1, 1-5 *************************************
************************************************************************************************

** p=-1
use "Data/PSP_Sample_BeforeMatched_t0", clear

* sample selection
drop if lln_kwh11==. /* drop the sample if the sample loses outcome variable information in leading periods s=i */
xtset system_code year
bysort system_code (year) : drop if _N < (year[_N] - year[1] + 1) 

* use logit model to estimate propensity score
qui: logit START $syvar $svar $cvar $cyvar $pvar i.year
predict phat if e(sample), pr 
drop if phat==.

* rescale estimated propensity score
egen idyr= group(year)
gen double phat_exact=phat+(idyr-1)*10 /* (idyr-1) part can reflect the difference of year */
su phat_exact phat

* kernel matching
psmatch2 START, kernel pscore(phat_exact) com outcome(lln_kwh11)
drop if _support==0
g esample=e(sample) /* e(sample) equals to 1 if the observation is in the estimation sample and 0 otherwise, thus we can know whether observation is in the matching sample */
drop if _w==. /* drop unmatched sample */
save "Data/PSP_Sample_Matched_t-1", replace

* balance property checking
* store mean bias and chi_square joint significance test
putexcel set "pstest_indicator", sh("PSIND_t-1") modify
pstest $syvar $svar $cvar $cyvar $pvar
putexcel A1 = "Test Indicator"
putexcel A2 = "Median Bias"
putexcel A3 = "Prob(chi2)"
putexcel B1 = "Value"
putexcel B2 = `r(meanbias)', nformat(number_d3)
putexcel B3 = `r(chiprob)', nformat(number_d3)
putexcel C1 = "Period"
putexcel C2 = "s=-1"
putexcel C3 = "s=-1"
putexcel close
* store t-value test result
putexcel set "pstest_indicator", sh("pst_t-1_balance") modify
putexcel A1 = "Variable"
putexcel B1 = "Difference"
putexcel C1 = "T-value"
putexcel D1 = "p-value"
local j=2
while `j'<=9{
foreach var in $syvar $svar $cvar $cyvar $pvar{
svyset _n [pw=_w]
putexcel A`j' = "`var'"
svy: mean `var', over(START) coeflegend
lincom _b[`var':1]-_b[`var':0] 
putexcel B`j' = `r(estimate)'
putexcel C`j' = `r(t)'
putexcel D`j' = `r(p)'
local j=`j'+1
}
}
putexcel close

* store number of treated and caliper
g treated=1 if START==1
g N=_N
tabstat treated N, s(N mean) c(s) save
tabstatmat all
matrix Stats=all'
xml_tab Stats , save("Baseline_Treated.xml") append sheet(obsT_t-1) 
	
* store mean of treated and control group
svyset _n [pw=_w]
svy: mean $syvar $svar $cvar $cyvar $pvar if  START==1
est store treated
svy: mean $syvar $svar $cvar $cyvar $pvar if  START==0
est store control
	
xml_tab treated control, save("Baseline_Treated.xml") append sheet(PST_t-1)
est clear


** p=1-5
forvalues i=1/5{

use "Data/PSP_Sample_BeforeMatched_t0", clear

* sample selection
drop if fln_kwh1`i'==. /* drop the sample if the sample loses outcome variable information in leading periods s=i */
xtset system_code year
bysort system_code (year) : drop if _N < (year[_N] - year[1] + 1) 

* use logit model to estimate propensity score
qui: logit START $syvar $svar $cvar $cyvar $pvar i.year
predict phat if e(sample), pr 
drop if phat==.

* rescale estimated propensity score
egen idyr= group(year)
gen double phat_exact=phat+(idyr-1)*10 /* (idyr-1) part can reflect the difference of year */
su phat_exact phat

* kernel matching
psmatch2 START, kernel pscore(phat_exact) com outcome(fln_kwh1`i')
drop if _support==0
g esample=e(sample) /* e(sample) equals to 1 if the observation is in the estimation sample and 0 otherwise, thus we can know whether observation is in the matching sample */
drop if _w==. /* drop unmatched sample */
save "Data/PSP_Sample_Matched_t`i'", replace

* balance property checking
* store mean bias and chi_square joint significance test
putexcel set "pstest_indicator", sh("PSIND_t`i'") modify
pstest $syvar $svar $cvar $cyvar $pvar
putexcel A1 = "Test Indicator"
putexcel A2 = "Median Bias"
putexcel A3 = "Prob(chi2)"
putexcel B1 = "Value"
putexcel B2 = `r(meanbias)', nformat(number_d3)
putexcel B3 = `r(chiprob)', nformat(number_d3)
putexcel C1 = "Period"
putexcel C2 = "s=`i'"
putexcel C3 = "s=`i'"
putexcel close
* store t-value test result
putexcel set "pstest_indicator", sh("pst_t`i'_balance") modify
putexcel A1 = "Variable"
putexcel B1 = "Difference"
putexcel C1 = "T-value"
putexcel D1 = "p-value"
local j=2
while `j'<=9{
foreach var in $syvar $svar $cvar $cyvar $pvar{
svyset _n [pw=_w]
putexcel A`j' = "`var'"
svy: mean `var', over(START) coeflegend
lincom _b[`var':1]-_b[`var':0] 
putexcel B`j' = `r(estimate)'
putexcel C`j' = `r(t)'
putexcel D`j' = `r(p)'
local j=`j'+1
}
}
putexcel close

* store number of treated and caliper
g treated=1 if START==1
g N=_N
tabstat treated N, s(N mean) c(s) save
tabstatmat all
matrix Stats=all'
xml_tab Stats , save("Baseline_Treated.xml") append sheet(obsT_t`i') 
	
* store mean of treated and control group
svyset _n [pw=_w]
svy: mean $syvar $svar $cvar $cyvar $pvar if  START==1
est store treated
svy: mean $syvar $svar $cvar $cyvar $pvar if  START==0
est store control
	
xml_tab treated control, save("Baseline_Treated.xml") append sheet(PST_t`i')
est clear

}

* Compile the balancing property checking results from period=-1 to 5
forvalues i=-1/5{
import excel using "pstest_indicator", sheet("PSIND_t`i'") clear firstrow
save "pst_`i'.dta", replace
}

use "pst_-1.dta", clear
forvalues i=0/5{
append using "pst_`i'.dta"
}
order Period TestIndicator Value 
export excel using "pstest_baseline.xlsx", sheet("PSIND") replace firstrow(variables)

forvalues i=-1/5{
import excel using "pstest_indicator", sheet("pst_t`i'_balance") clear firstrow
gen Period=`i'
save "pst_`i'.dta", replace
}

use "pst_-1.dta", clear
forvalues i=0/5{
append using "pst_`i'.dta"
}
order Period
export excel using "pstest_baseline.xlsx", sheet("T-test") sheetmodify firstrow(variables)

erase "pstest_indicator.xlsx"
forvalues i=-1/5{
erase "pst_`i'.dta"
}

*****************************************************************************************************
************************************ Step 6: Calculate ATT ******************************************
*****************************************************************************************************

// Immediate Effect

* ATT - p=-1
putexcel set "ATT_period", sh("ATT_t-1") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"

putexcel A2 = -1

use "Data/PSP_Sample_Matched_t-1", clear
qui: reghdfe lln_kwh11 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_t_1
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
local numb= e(N)
putexcel G2 = `numb'

putexcel close

* ATT - p=0
putexcel set "ATT_period", sh("ATT_t0") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"

putexcel A2 = 0

use "Data/PSP_Sample_Matched_t0", clear
qui: reghdfe ln_kwh1 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_t0
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
local numb= e(N)
putexcel G2 = `numb'

putexcel close

* ATT - p=1 to p=5
forvalues i=1/5{
putexcel set "ATT_period", sh("ATT_t`i'") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"

putexcel A2 = `i'

use "Data/PSP_Sample_Matched_t`i'", clear
qui: reghdfe fln_kwh1`i' START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_t`i'
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
local numb= e(N)
putexcel G2 = `numb'

putexcel close
}

* compile ATT results from p=-1 to 5
forvalues i=-1/5{
import excel using "ATT_period", sh("ATT_t`i'") firstrow clear
save temp_`i', replace
}
use temp_-1, clear
forvalues i=0/5{
append using temp_`i'
}
export excel using "ATT_period.xlsx", firstrow(variables) replace

* draw ATT dynamic graph (current effect)
import excel using "ATT_period.xlsx", firstrow clear
set scheme s1mono
format beta %9.3f
twoway (connected beta Period, lcolor(black) lpattern(solid) lw(thick) mlcolor(black) mfcolor(black) msize(*1.2) mlabel(beta) mlabp(6) mlabs(small)) ///
	(line ci_l Period, lcolor(black) lpattern(dash) lw(medium)) ///
	(line ci_r Period, lcolor(black) lpattern(dash) lw(medium)) ///
	, graphregion(fcolor(white)) legend(off) title("Roadshow's Dynamic Effect") /// 
	ytitle("Roadshow Effect (%)", size(medium)) yline(0, lp(solid) lw(medthin) lc(red)) ylab(-0.1(0.05)0.2,nogrid) ///
	xtitle("Period", size(medium)) xline(0, lp(dash) lw(thin) lc(red)) xlab(-1(1)5,labsize(medium)) ///
	note(" " " Notes: Dash Lines Represent 95% Confidence Interval of Each Point Estimate.", size(small)) 
graph export "ATT_period.png", as(png) replace


// Cumulative Effect
* p=1
use "Data/PSP_Sample_Matched_t1", clear
putexcel set "ATT_cum", sh("ATT_t1") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"
putexcel A2 = 1
g cum1 = ln_kwh1 + fln_kwh11 
qui: reghdfe cum1 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_cum_t1
local numb= e(N)
putexcel G2 = `numb'
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
putexcel close
* p=2
use "Data/PSP_Sample_Matched_t2", clear
putexcel set "ATT_cum", sh("ATT_t2") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"
putexcel A2 = 2
g cum2 = ln_kwh1 + fln_kwh11 + fln_kwh12
qui: reghdfe cum2 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_cum_t2
local numb= e(N)
putexcel G2 = `numb'
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
putexcel close
* p=3
use "Data/PSP_Sample_Matched_t3", clear
putexcel set "ATT_cum", sh("ATT_t3") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"
putexcel A2 = 3
g cum3 = ln_kwh1 + fln_kwh11 + fln_kwh12 + fln_kwh13
qui: reghdfe cum3 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_cum_t3
local numb= e(N)
putexcel G2 = `numb'
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
putexcel close
* p=4
use "Data/PSP_Sample_Matched_t4", clear
putexcel set "ATT_cum", sh("ATT_t4") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"
putexcel A2 = 4
g cum4 = ln_kwh1 + fln_kwh11 + fln_kwh12 + fln_kwh13 + fln_kwh14
qui: reghdfe cum4 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_cum_t4
local numb= e(N)
putexcel G2 = `numb'
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
putexcel close
* p=5
use "Data/PSP_Sample_Matched_t5", clear
putexcel set "ATT_cum", sh("ATT_t5") modify
putexcel A1 = "Period"
putexcel B1 = "beta"
putexcel C1 = "se"
putexcel D1 = "pvalue"
putexcel E1 = "ci_left"
putexcel F1 = "ci_right"
putexcel G1 = "obs"
putexcel A2 = 5
g cum5 = ln_kwh1 + fln_kwh11 + fln_kwh12 + fln_kwh13 + fln_kwh14 + fln_kwh15
qui: reghdfe cum5 START $syvar $svar $cvar $cyvar $pvar i.state_code#c.year [pw=_weight], absorb(year county) cluster(system_code)
est store psps_cum_t5
local numb= e(N)
putexcel G2 = `numb'
local beta= _b[START]
putexcel B2 = `beta'
local se= _se[START]
putexcel C2 = `se'
local df=  e(df_r)
local p= ttail(`df', abs(`beta'/`se'))*2
putexcel D2 = `p'
local crt= invttail(`df', 0.05)
local ci_l= `beta' - `se'*`crt'
putexcel E2 = `ci_l'
local ci_r= `beta' + `se'*`crt'
putexcel F2 = `ci_r'
putexcel close
	 
* compile results from p=1 to 5
forvalues i=1/5{
import excel using "ATT_cum", sh("ATT_t`i'") firstrow clear
save temp_`i', replace
}
use temp_1, clear
forvalues i=2/5{
append using temp_`i'
}
export excel using "ATT_cum.xlsx", firstrow(variables) replace

* draw ATT dynamic graph (cumulative effect)
import excel using "ATT_cum.xlsx", firstrow clear
set scheme s1mono
format beta %9.3f
twoway (connected beta Period, lcolor(black) lpattern(solid) lw(thick) mlcolor(black) mfcolor(black) msize(*1.2) mlabel(beta) mlabp(12) mlabs(small)) ///
	(line ci_l Period, lcolor(black) lpattern(dash) lw(medium)) ///
	(line ci_r Period, lcolor(black) lpattern(dash) lw(medium)) ///
	, graphregion(fcolor(white)) legend(off) title("Roadshow's Cummulative Effect") /// 
	ytitle("Roadshow Cummulative Effect (%)", size(medium)) yline(0, lp(solid) lw(medthin) lc(red)) ylab(,nogrid format(%9.1f)) ///
	xtitle("Period", size(medium)) xlab(1(1)5,labsize(medium)) ///
	note(" " " Notes: Dash Lines Represent 95% Confidence Interval of Each Point Estimate.", size(small)) 
graph export "ATT_cum.png", as(png) replace

forvalues i=0/5{
erase "temp_`i'.dta"
}

** Our result is robust to cluster on borrower-level and use White robust standard error

** Q: ATT robust to distance thereshold? 
** Answer: It seems that taking 30, 35, 40, 45 miles as threshold, the graph is pretty similar (significant effect in about three periods after the tour);
**         However, when 20 & 50 miles is used as threshold, immediate effect in current periods (p=0) is pretty large and the effect diminishes in one period.

** Q: ATT effect difference between FET stop and all demonstration stops
** Answer: under all distance thresholds, FET stop has larger effect! (e.g. under 30 miles, in p=0, FET 0.114 v.s. all 0.097
