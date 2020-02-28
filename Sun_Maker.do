******************************************************
** Sun et al. (2018) IW estimation implementation   **
** - This is an automatization of the process 		**
**	taking as guidance the exercise carried by 		**
**	Sun et al. for the hospitals example.			**
**													**
** By: David Vargas | davidl.vargas@urosario.edu.co	** 
**								   					**
** This code performs the dynamic FE (as Sun does)	**
** IW estimate (as Sun does)	 					**
** Plot of estimates (as Sun does) 					**
******************************************************

clear all

*** path to files

** data 
global data "C:\Users\David\Dropbox\Mounu\Even Study\"

** output (to save results)
global output "C:\Users\David\Dropbox\Mounu\Even Study\Output"

*** Needed definitions

/** Replicate Sun? If so change global to YES, and be sure to have HRS_long.dta 
on your data folder **/
global replicate_sun "NO"

/** Are you using your data? Fill the following globals **/

** load your data 
use "${data}/ADPREG_TEMP.dta", clear
/* Data preparation */ 

* Centered time var
by rbd: egen min_adp = min(year) if adp == 1
	spread min_adp, by(rbd)
g t = year - min_adp 

* sample selection (since 2012 are some treated)
drop if year < 2011

* ever treated
gen ever_treat = (min_adp != .)

/* End data preparation */

** Dependant variable
global dep "sdWfe_prin"

** Centered time definition (ex. t = -3, -2, -1, 0, 1, 2, 3)
** In the first all should be controls and in the last all treatments
global evt_time "t"
	
** Time var (year, survey wave, etc..) [This are considered consecutive (if not gen a time period)]
global time_var "year"

** Dummy for ever treated
global ever_treat "ever_treat" 

** ID for unit of analysis (for cluster and FE)
global unit_ID "rbd"


********************************************************************************
********************************************************************************
********************************************************************************

/*** In case replicate_sun == YES this will run to load the data and stablish 
the correct globals */

if "${replicate_sun}" == "YES" {

	use "${data}/HRS_long.dta", clear
	
	*Preliminaries
	drop if wave < 7 // keep a balanced sample for wave 7-11
	bys hhidpn: gen N = _N
	keep if N == 5
	
	** Dependant variable
	global dep "oop_spend"

	** Centered time definition (ex. t = -3, -2, -1, 0, 1, 2, 3)
	** In the first all should be controls and in the last all treatments
	global evt_time "evt_time"
		
	** Time var (year, survey wave, etc..) [This are considered consecutive]
	global time_var "wave"

	** Dummy for ever treated
	global ever_treat "ever_hospitalized" 

	** ID for unit of analysis (for cluster and FE)
	global unit_ID "hhidpn"
}

********************
* Some preparation *
********************

set matsize 800

** Keep those that were treated
keep if ${ever_treat} == 1


** Drop those foe which we have no controls
bys ${unit_ID}: egen flag = min(${evt_time})
drop if flag >= 0 & flag != .
drop if flag == . 
drop flag


** Define time of treatment
gen time_t = $time_var if $evt_time == 0
bys ${unit_ID}: egen time_treat = min(time_t)
drop time_t

** Cohort dummies
tab time_treat, gen(time_treat_)

** Calendar dummies
xi i.${time_var}
sum ${time_var}
local max_time_var = r(max)
global max_time_var `max_time_var'
global min_time_var `r(min)'

** save the number of the last wave
sum time_treat
local last_wave_treat = r(max)
global last_wave_treat `last_wave_treat'


** l wave relative to treatment [event time](dummies)
tab $evt_time, gen(evt_time_)
local N_waves = r(r)
global N_waves `N_waves'

** The omited wave dummy is the one at t= -1 (As Sun does)
** define wave to exclude
local counter = 0
foreach v of varlist evt_time_*{
	local counter = `counter' + 1
	qui tab ${evt_time} if `v' == 1, matrow(val)
	if val[1,1] == -1 {
		global time_excluded "`v'"
	}
	if val[1,1] == 0 {
		global time_zero "`v'"
		global evt_zero "`counter'"
	}
}

** global of waves

** all waves 
global waves
foreach v of varlist evt_time_*{
	global waves $waves `v'
}

** w/out excluded waves
* excluding t = -1
local w_excl = subinstr("$waves", "$time_excluded ", "",.) 
* The first wave not controls --> is exclueded
local w_excl = subinstr("`w_excl'", "evt_time_1 ", "",.) 
* global definition
global waves2 "`w_excl'"
local N_waves_e = `N_waves' - 2
global N_waves_e `N_waves_e'


*** Generate counts for each cohort

qui levelsof time_treat
global time_t_levels `r(levels)'
foreach i of global time_t_levels {
	di `i'
	count if time_treat == `i' & ${time_var} == $min_time_var
	global c_`i' = r(N)
}

if "${replicate_sun}" == "YES" {
	/* IMPORTANT NOTICE: For an unclear reason Sun capped the studied 
	sample after calculating the cohort counts for the weights. 
	If you're replicating Sun and want to replicate
	as she does uncomment the next line */
	keep if age_hosp <= 59 //[SUN]
}

*****************
*** Estimation **
*****************

*************** Regular linear two-way dynamic FE **********************
reghdfe ${dep}  ${waves2} _I${time_var}_*  if ${ever_treat} , absorb(${unit_ID}) cluster(${unit_ID})

// store results
mat b = e(b)
mat b = b[.,1..$N_waves_e]
mat V = vecdiag(e(V))
mat V = V[.,1..$N_waves_e]	
mat b2way_hrs = b', V'

*************** Interaction weighted (IW) estimator *******************

** get variance from estimating the share of cohorts

// exclude the last cohort because it will be used as the control units
foreach vv of varlist evt_time_* {
	replace `vv' = 0 if time_treat == $last_wave_treat
}

preserve
foreach vv of varlist evt_time_* {
	replace `vv' = . if time_treat == $last_wave_treat
}  
estimates clear

***** Estimation *****

*** Preparation: 
// global of estimates to fill
global estimates1

* exclude waves at the "corner" waves
// exclude the last wave 
foreach v of varlist evt_time_*{
	local last = "`v'"
}
di "`last'"
local w_excl = subinstr("$waves", "`last'", "",.) 
// exclude first and second wave
local w_excl = subinstr("`w_excl'", "evt_time_1 ", "",.) 
local w_excl = subinstr("`w_excl'", "evt_time_2 ", "",.) 
// global definition
global waves3 "`w_excl'"
local N_waves_e2 = $N_waves_e - 1
global N_waves_e2 `N_waves_e2'
// The excluded time_treat_? is the last:
global times_t 
foreach v of varlist time_treat_* {
	local last = "`v'"
	global times_t "${times_t} `v'"
}
local tt_excl = subinstr("$times_t", "`last'", "",.) 
global times2 "`tt_excl'"


*** Estimation 
local counter = 0
foreach time_v of global times2{
	local counter = `counter' + 1
	reg `time_v' ${waves3} if ${time_var} < ${max_time_var}, nocons
	estimates store est_`counter'
	global estimates1 "${estimates1} est_`counter'"
}

global n_est `counter'

suest ${estimates1}, vce(robust)

// to store variance covariance 
matrix V = e(V)
//Extract matrix, install package "matselrc"
	// Not include constant 
local low = 1
local limit = $N_waves_e2
global mat_select 
local counter = 0
foreach time_v of global times2{
	local counter = `counter' + 1
	if `counter' == 1{ 
	global mat_select "${mat_select} `low'/`limit'" 
	}
	else{
	global mat_select "${mat_select},`low'/`limit'" 
	}
	local low = `limit' + 2
	local limit = `limit' + $N_waves_e2 + 1
}
di "$mat_select"
	// Extract matrix 
matselrc V Vcov, r("$mat_select") c("$mat_select")

restore

** Saturated specification
** estimate a linear two-way FE specification that saturates in relative time indicators with cohort indicators

// determine the posible treatment period range
sum time_treat
local last_treat = r(max)
local first_treat = r(min)
global last_treat `last_treat'
global first_treat `first_treat'

local counter = 0
forval i = `first_treat'/`last_treat' {
local counter = `counter' + 1
}

local re_number_groups = `counter' - 1
global re_number_groups `re_number_groups'

* Gen the interactions and determine subset of factible interactions
// Factible interactions have no zero mean in the times before the last period

// clean macros to avoid errors
global factible_inter 
forval i = `first_treat'/`last_treat' {
	global factible_inter_`i'
}
// Calculate interaction and determine feaseble interactions by treatment time
foreach vv of varlist evt_time_* {
	forval i = `first_treat'/`last_treat' {
		gen `vv'_`i' = `vv' * (time_treat == `i')
		sum `vv'_`i' if ${time_var} < ${max_time_var}
		if "`vv'" != "$time_excluded" {
			local mean = r(mean)
			if `mean' > 0{
				global factible_inter_`i' "${factible_inter_`i'} `vv'_`i'"
			}
		}
	}
}
// join feaseable in one macro
forval i = `first_treat'/`last_treat' {
	global factible_inter "${factible_inter} ${factible_inter_`i'}"
}

*** CATT estimation 
reghdfe ${dep} ${factible_inter} _I${time_var}_*  if ${ever_treat} & ${time_var} < ${max_time_var} , absorb(${unit_ID}) cluster(${unit_ID})

// store results
mat b = e(b)
mat V = vecdiag(e(V))
mat Bb = b
mat Vv = V

// store results a matrix of estimates (variance) per cohort
local jlim = $N_waves - 1
forval i = $first_treat/$last_treat {
global callb
local counter = 0
	forvalues j = 2/`jlim' {
	local vv = "evt_time_`j'"
		local counter = `counter' + 1
		local test = "`vv'_`i'"
		cap di _b["`test'"]
		if _rc == 0 {
			di "`test'"
			local coln = colnumb(Bb, "`test'")
			if "`counter'" == "1"{
				global callb "b[.,`coln'..`coln']"
				global callv "V[.,`coln'..`coln']"
			}
			else {
				global callb "${callb},Bb[.,`coln'..`coln']"
				global callv "${callv},Vv[.,`coln'..`coln']"
			}
		}
		else if "`vv'" == "$time_excluded" {
			if "`counter'" == "1"{
				global callb "0"
				global callv "0"
			}
			else {
				global callb "${callb},0"
				global callv "${callv},0"
			}
		}
		else{
			if "`counter'" == "1"{
				global callb "."
				global callv "."
			}
			else {
				global callb "${callb},."
				global callv "${callv},."
			}
		}
	}
	mat b_`i' = (${callb})
	mat V_`i' = (${callb})
}

**** Weighted averages ****
// take weighted average of estimates from the saturated specification to form IW estimates
// variance estimate for IW estimators follow from Proposition 6 in the appendix to the main text

mat biw_long = b // copy the coefficients from the saturated specification, to add to variance estimates later

* Gen iterables to ease (automate) estimation later on
local jlim = $N_waves - 1
forvalues j = 2/`jlim' {
global coho_`j'
global count_`j'
global group_`j'
global fact_var_`j'
local counter = 0
	local vv = "evt_time_`j'"
	forval i = ${first_treat}/${last_treat}{
	local counter = `counter' + 1
		local test = "`vv'_`i'"
		if strpos("$factible_inter", "`test'"){
			di "`test'"
			global coho_`j' "${coho_`j'}  `test'"
			global group_`j' "${group_`j'} ${c_`i'}"
			if "`i'" != "${last_treat}" {
			global fact_var_`j' "${fact_var_`j'} est_`counter'_mean:`vv'"
			}
			if "${count_`j'}" == ""{
				global count_`j' "${c_`i'}"
			}
			else{
				global count_`j' "${count_`j'} + ${c_`i'}"
			}
		}
	}
	di "${coho_`j'}"
	di "${count_`j'}"
	di "${group_`j'}"
	di "${fact_var_`j'}"
}


/* For each evt_time_? determine the coeffients, 
the variance of each and the weights given the cohorts */
/* IMPORTANT NOTICE: Sun seems to have an error for 
her evt_time_5 IW estimate. Takes the wrong column for the 10th wave.
She takes 8, the correct one is the 9th. */

** select coefficients and variances to withdraw
local jlim = $N_waves - 1
forvalues j = 2/`jlim' {
local counterb = 0
local counter_v = 0
global tempb_call`j'
global tempVcov_call`j'

	foreach coh of global coho_`j' {
		local counterb = `counterb' + 1
		local coln = colnumb(biw_long,"`coh'")
		if "`counterb'" == "1"{
			global tempb_call`j' "`coln'"
		}
		else {
			global tempb_call`j' "${tempb_call`j'},`coln'"
		}
	}
	foreach coh_v of global fact_var_`j' {
		local counter_v = `counter_v' + 1
		local coln_v = colnumb(Vcov, "`coh_v'")
		if "`counter_v'" == "1" {
			global tempVcov_call`j' "`coln_v'"
		}
		else {
			global tempVcov_call`j' "${tempVcov_call`j'}, `coln_v'"
		}
	}
	
di "c(${tempb_call`j'})"
di "${tempVcov_call`j'}"
}

** Define weights for linear combination 
local jlim = $N_waves - 1
forvalues j = 2/`jlim' {
	local i = 0
	foreach size in ${group_`j'} {
		local i = `i' + 1 
		global weights_`j'_`i' "`size'/(${count_`j'})"
		global Ngroups_`j' = `i'
		di "${weights_`j'_`i'}"
	}
	di "${Ngroups_`j'}"
}

** Define linear combination for each time 
local jlim = $N_waves - 1
forvalues j = 2/`jlim' {
	if "${Ngroups_`j'}" == "1" {		
		global l_com_`j' "${coho_`j'}"
	}
	else {
		local i = 0
		foreach coh in ${coho_`j'} {
			local i = `i' + 1
			if "`i'" == "1" {
				global l_com_`j' "${weights_`j'_`i'} * `coh'"
			}
			else {
				global l_com_`j' "${l_com_`j'} + ${weights_`j'_`i'} * `coh'"
			}
		}
	}
	di "${l_com_`j'}"
}

**** Actual IW estimation ****
* Matrix call to save results
global biw_hrs_in
local last_bm = $last_treat - 1
forval i = $first_treat/`last_bm' {
	global biw_hrs_in "${biw_hrs_in},b_`i'',V_`i''"
}

* estimation 
local jlim = $N_waves - 1
forvalues j = 2/`jlim' {
	local vv = "evt_time_`j'"
	di "Estiamte for `vv'" 
	if "`vv'" == "$time_excluded" {
		matrix b = b \ 0
		matrix V = V \ 0
	}
	else{
		if "${Ngroups_`j'}" == "1" {
			lincom "${l_com_`j'}"
			if "`j'" == "2" {
				matrix b = r(estimate)
				matrix V = r(se)^2
			}
			else {
				matrix b = b \ r(estimate)
				matrix V = V \ r(se)^2  
				if "`j'" == "`jlim'"{
					mat biw_hrs = b,V${biw_hrs_in}
				}
			}
		}
		else {
			* extract betas
			matselrc biw_long tempb, c(${tempb_call`j'})
			* extract variance
			matselrc Vcov tempVcov, r(${tempVcov_call`j'}) c(${tempVcov_call`j'})
			* lincom 
			matrix temp = tempb*tempVcov*tempb'
			lincom "${l_com_`j'}"
			matrix b = b \ r(estimate)
			matrix V = V \ (r(se)^2 + temp)
			* save results on a matrix 
			if "`j'" == "`jlim'"{
			di "b,V${biw_hrs_in}"
				mat biw_hrs = b,V${biw_hrs_in}
			}
		}
	}
} 


********* output results *********

local jlim = $N_waves - 1
local onepre = $evt_zero - 1
clear
global n_obs `jlim'
set obs $n_obs
gen t = _n - `onepre'

**** Extract FE results 

* add a zero at the excluded t = -1 in the FE 
local rown = rownumb(b2way_hrs,"$time_zero")
di "`rown'"
local rown_b = `rown' - 1
local rown_a = `rown'

mat b2way_hrs = b2way_hrs[1..`rown_b',.]\(0,0)\b2way_hrs[`rown_a'..${N_waves_e},.]
* call matrix
svmat b2way_hrs
* The first column is the estimate the second the variance
rename b2way_hrs1 FE_b
rename b2way_hrs2 FE_v
gen FE_ci_u = sqrt(FE_v)*1.96+FE_b
gen FE_ci_l = -sqrt(FE_v)*1.96+FE_b
gen FE_sd = sqrt(FE_v)


**** Extract IW 
svmat biw_hrs
rename biw_hrs1 IW_b
rename biw_hrs2 IW_v

gen IW_ci_u = sqrt(IW_v)*1.96+IW_b
gen IW_ci_l = -sqrt(IW_v)*1.96+IW_b
gen IW_sd = sqrt(IW_v)

**** CATT results 
local counter_b = 1
local counter = 0 
local last_bm = $last_treat - 1
forval i = $first_treat/`last_bm' {
	local counter = `counter' + 1
	local counter_b = `counter_b' + 2
	local counter_v = `counter_b' + 1
	rename biw_hrs`counter_b' CATT`counter'_b
	rename biw_hrs`counter_v' CATT`counter'_v
	gen CATT`counter'_sd = sqrt(CATT`counter'_v)
}

* save results to excel
export delimited using "${output}/Sun_Results.csv", replace

*****************************************************************
******************************************************************

**********************
** Estimates Figure **
**********************

* Simple ajustment to make results visible
gen t_FE = t - 0.05
gen t_IW = t + 0.05

* global for plot
program scatterwCI 
	args var pat symb
{
	global `var'_scatterwCI "rcap `var'_ci_u `var'_ci_l t_`var', lw(medthick) lpattern(`pat') || scatter `var'_b t_`var' , c(l) msize(medlarge) m(`symb') "
}
end

* Plot 
scatterwCI FE dash T
scatterwCI IW line O
twoway ${FE_scatterwCI} || ${IW_scatterwCI} ///	
	legend(order(2 "FE" 4 "IW")) ///
	xtitle("Relative Wave") ytitle("Estimate") ///
	graphregion(fcolor(white)) scheme(sj)
	
graph export "${output}/Sun_Graph.pdf", replace
