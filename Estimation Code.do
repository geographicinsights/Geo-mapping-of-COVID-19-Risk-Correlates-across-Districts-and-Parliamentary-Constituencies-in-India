*Stata/MP version 14.1/16.0/16.1
*August 17, 2020
*Geo-mapping of COVID-19 Risk Correlates in India: A Data Byte and Visualization

*NOTE: the execution time of the code is around 1 hour










/* Select project directory (where all files needed should be stored and 
all files created here will be stored) */

global dir ///
"c:\download" // E.G.








/*******************************************************************************
********************************************************************************
********************************************************************************

Downloading the datasets

There should be 2-3 data files in the project directory before running this code

1) IAPR74DT.zip with the household survey (household member recode) from dhsprogram.com
The file is accessible after a simple registration from dhsprogram.com:
https://dhsprogram.com/customcf/legacy/data/download_dataset.cfm?Filename=IAPR74DT.ZIP&Tp=1&Ctry_Code=IA&surv_id=355&dmode=cf

AND (if area.csv is not available)

2a) sdr_subnational_boundaries_[date suffix].zip with the district boundaries for 
DHS India 2015-16 from spatialdata.dhsprogram.com:
https://spatialdata.dhsprogram.com/boundaries/#view=table&countryId=IA

OR (recommended!)

2b) area.csv is available as a supplementary file with this paper, where 
area in square meters for each district has been calculated from file 2 (sdr_subnational_boundaries_[date suffix].zip)
If this file is not in the project directory, a block of code will estimate the 
approximate area from file 2 instead

AND (Optional)

3) "District to PC translator.csv" (which is only used to create estimates for 
Parliamentary Constituencies)



********************************************************************************
********************************************************************************
*******************************************************************************/










qui { // remove this line if you want to see the output

capture noi cd "$dir"
if _rc != 0 {
noi dis "NOTE! Enter a project directory at the top of this do-file"
x
}

set matsize 799 // the maximum number of variables

* These programs are needed
ssc install mat2txt
net install "http://www.stata.com/stb/stb56/dm79"
ssc install shp2dta
net install "http://www.stata-journal.com/software/sj15-2/dm0080.pkg"
ssc install distinct
ssc install reclink


********************************************************************************
********************************************************************************
********************************************************************************


capture confirm file area.csv
if _rc != 0 {

#delimit;	
noi dis "NOTE! If the file area.csv is not present in the project directory ($dir), 
the area of a district will be calculated using this block of code. 
It is, howvever, highly recommended to rather download the area.csv file available 
with this paper, or use a GIS software to create area.csv from 
sdr_subnational_boundaries2.shp, since the fieldarea command in Stata only 
calculates approximate district size. A file titled 'area.csv' containing area 
calculated from sdr_subnational_boundaries2.shp is available as a supplement 
to this paper";
#delimit cr

local file: dir "$dir" files "sdr_subnational_boundaries_*.zip"
copy `file' "boundries.zip", replace
unzipfile boundries.zip, replace
shp2dta using "shps\sdr_subnational_boundaries2.shp", database(dbf)  coordinates(xy) replace
use xy, clear
fieldarea _X _Y, generate(area) id(_ID)
save area.dta , replace
use dbf, clear
merge 1:1 _ID using area.dta, nogen
gen state = OTHREGNA
gen district= REGNAME
}

if _rc == 0  import delimited area , clear
save area.dta, replace

********************************************************************************
/*** Generating sampling weights, household ID, and and indicator for the 
non-pregnant female 15-49 sub-sample */

unzipfile "IAPR74DT.zip" , replace
use IAPR74FL.dta, clear 

gen wt = hv005/1000000
egen hid = group(hhid)
gen biom = hv104 == 2 & inrange(hv105, 15, 49) & ha54 == 0

********************************************************************************
/*** Generating non-clinical risk correlates for the total population */

gen _Nover65 = hv105>65 if hv105 < 98
gen _Nnowash = hv232+hv230b!=2 if hv232+hv230b < .
replace _Nnowash = 1 if hv230a == 2
gen _Ncrowded =  (hv009/hv216>3) | (hv216>3 & hv216 == 0) if  hv009+hv216<.


********************************************************************************
/*** Generating hypertension for non-pregnant females 15-49 years old */

/* systolic blood pressure */
gen sys1 =  shb16s if  inrange(shb16s,60,250) & hv104 == 2
gen sys2 =  shb23s  if inrange(shb23s,60,250) & hv104 == 2
gen sys3 =  shb27s  if inrange(shb27s,60,250) & hv104 == 2
gen msys = (sys1+sys2+sys3)/3
replace msys = (sys1+sys2)/2  if msys ==.
replace msys = (sys1+sys3)/2  if msys ==.
replace msys = (sys2+sys3)/2  if msys ==.
replace msys = sys3 if msys ==.
replace msys = sys2 if msys ==.
replace msys = sys1 if msys ==.

/* diastolic blood pressure */
gen dia1 = shb16d if inrange(shb16d,40,140) & hv104 == 2
gen dia2 = shb23d if inrange(shb23d,40,140) & hv104 == 2
gen dia3 = shb27d if inrange(shb27d,40,140) & hv104 == 2
gen mdia = (dia1+dia2+dia3)/3
replace mdia = (dia1+dia2)/2  if mdia ==.
replace mdia = (dia1+dia3)/2  if mdia ==.
replace mdia = (dia3+dia2)/2  if mdia ==.
replace mdia = dia3 if mdia ==.
replace mdia = dia2 if mdia ==.
replace mdia = dia1 if mdia ==.

/* hypertension from https://dhsprogram.com/pubs/pdf/FR339/FR339.pdf page 4 */
gen _nhyper = msys > 140 | mdia > 90 if msys+mdia<. & ha54 == 0


********************************************************************************
/*** Generating obesity for non-pregnant females 15-49 years old */

gen kg = ha2/10 if inrange(ha11, -600,600)
gen cm = ha3/1000 if inrange(ha5, -600,600)

gen _nsobese = (kg/cm^2) >=30 if kg/cm<. & ha54 == 0


********************************************************************************
/*** Generating diabetes for non-pregnant females 15-49 years old */

gen fasted = shb51 > 8 &  shb52 > 8 if shb51+shb52 <. &  hv104 == 2
gen _ndiab =(shb70>=126 & fasted==1) | (shb70>=200 & fasted==0) if shb70 < 500 ///
& fasted<. &  hv104 == 2 & ha54 == 0


********************************************************************************
********************************************************************************

gen nbiom = _nhyper+_ndiab+_nsobese==.
keep hv024 shdistri hid nbiom  biom wt _Ncrowded _Nover65 _Nnowash _nhyper ///
_ndiab _nsobese hv003
compress
save d.dta , replace


********************************************************************************
/*** Estimating the number of people experiencing each risk factor in each 
district as a proportion of the total NFHS sample. Used for creating count measures
by multiplying each estimate by the total Indian population in mid-year 2015 */

use d.dta , clear
svyset hid [pweight = wt]
foreach var in _Nover65 _Ncrowded _Nnowash {
gen d`var' = shdistri if `var' == 1
replace d`var' = 9999 if `var' == 0
replace d`var' = . if `var' == .
svy: tabulate d`var', se  nolabel
matrix v = e(V)
local rows = rowsof(v)
matrix se = J(`rows', 1,.)
forvalues i = 1/`rows' {
matrix se[`i',1] = sqrt(v[`i',`i']) // this gives the se estimates from the svy: tabulate table
}
matrix z = e(Row)', e(b)', se
mat2txt , matrix(z) saving(dp`var') replace
drop d`var'
}


********************************************************************************
/*** Estimating the number of people in a district with the repective risk
factors as a proportion of the total number of people in a district. Used for 
creating the percentage measures */

use d.dta , clear
foreach var in _Nover65 _Ncrowded _Nnowash _nhyper _ndiab _nsobese {
if "`var'"  == "_nhyper" drop if nbiom == 1
mean `var' [pweight=wt] , cluster(hid) over(shdistri, nolabel) nolegend
matrix e`var' = r(table)'
mat2txt , matrix(e`var') saving(e`var') replace
}


********************************************************************************
/*** Estimating the total number of people in a district as a proportion of the 
total number in the NFHS sample. Used for creating a district population measure 
which was then used to calcualted population density */

use d.dta , clear
keep wt shdistri hid biom
svyset hid [pweight = wt]
svy: tabulate shdistri, se  nolabel
matrix v = e(V)
local rows = rowsof(v)
matrix se = J(`rows', 1,.)
forvalues i = 1/`rows' {
matrix se[`i',1] = sqrt(v[`i',`i']) // this gives the se estimates from the svy: tabulate table
}
matrix z = e(Row)', e(b)', se
clear
svmat2 z, names(shdistri tpd_b tpd_se)
save raw_estimates.dta, replace


********************************************************************************
/*** Combining all estimates into one file */

foreach var in _Nover65 _Ncrowded _Nnowash _nhyper _ndiab _nsobese {
import delimited "e`var'.txt", varnames(1) clear 
local lab =substr("`var'",3,length("`var'"))
for var b se t pvalue ll ul df crit eform: rename X `lab'_e_X
replace v1 = subinstr(v1, "`var'", "",.)
replace v1 = subinstr(v1, ":", "",.)
gen shdistri = ""
egen l = max(length(v1)) 
sum l
forval i = 1/`r(max)' {
replace shdistri=shdistri+substr(v1,`i',1) if real(substr(v1,`i',1))<.
}
destring shdistri, replace
merge 1:1 shdistri using raw_estimates.dta , nogen
save raw_estimates.dta, replace
}

foreach var in _Nover65 _Ncrowded _Nnowash {
import delimited "dp`var'.txt", varnames(1) clear 
rename (c1 y1 v4)(shdistri b se)
local lab =substr("`var'",3,length("`var'"))
for var b se: rename X `lab'_pd_X
merge 1:1 shdistri using raw_estimates.dta , nogen
save raw_estimates.dta, replace
}

********************************************************************************
/*** Summing up sample sizes and missing observations for each measure and for each
district, and adding to the estimates dataset */

use d.dta , clear
bys shdistri: gen  obst =_N
bys shdistri: egen obsf15to49  = total(biom)
egen tag = tag(shdistri hid)  
bys shdistri: egen  obst_hh =sum(tag)
drop tag

foreach var in over65 crowded nowash {
bys shdistri: egen  `var'_m =sum(_N`var'==.)
}
foreach var in hyper diab sobese {
bys shdistri: egen  `var'_m =sum(_n`var'==. & biom == 1)
}

keep shdistri obs* hv024 *_m
duplicates drop shdistri , force
merge 1:1 shdistri using raw_estimates.dta , nogen


********************************************************************************
/*** Adding district size in square meters to the data, to calcualte 
population density */

/* there are some differences in the district and state names between the PR file
and the district boundries file */
decode shdistri, gen(district)
decode hv024 , gen(state)
replace state = trim(subinstr(proper(state), " And " , " and ", .))
replace state = "Himanchal Pradesh"  if state == "Himachal Pradesh" 
replace state = trim(subinstr(state, " and " , " & " ,.)) if strpos(state, "Jammu")
replace district= trim(subinstr(proper(district), " And " , " and ", .))
drop if hv024==.

/* the size of each district is estimated from a shapefile called 
sdr_subnational_boundaries2 downloaded from 
https://spatialdata.dhsprogram.com/boundaries/#view=table&countryId=IA
titled India 2015 DHS */

merge 1:1 state district using area.dta, nogen keepusing(area)





save raw_estimates.dta, replace






********************************************************************************
********************************************************************************
********************************************************************************
/*** Cleaning and adding variable labels to create the 
final estimates file (Supplementary Table S1) */

use raw_estimates.dta, clear
replace area = area/1000^2
drop obs*  *_ul *_ll *_eform *_crit *_t *_pvalue  *_df  v* *_m l
foreach var of varlist *_b *_se {
if strpos("`var'", "_e_")  replace `var' = `var' * 100
}
order *_e_b *_e_se *_pd_b *_pd_se
order district state area crowded* nowash* over65* hyper* diab* sobese*

foreach var of varlist *_pd_b {
gen `var'_2 = `var' * 1310152000
label var `var'_2 "`:variable label `var''"
}

gen tpd_b_2 = tpd_b * 1310152000
gen density = tpd_b_2/area

foreach var of varlist * {
if strpos("`var'" , "density") label var `var' "Density (people/per square km)"
if strpos("`var'" , "over65")  label var `var' "Elderly" 
if strpos("`var'" , "nowash")  label var `var' "No handwashing"
if strpos("`var'" , "hyper")   label var `var' "Hypertension"
if strpos("`var'" , "diab")    label var `var' "Diabetes"
if strpos("`var'" , "sobese")  label var `var' "Obesity"
if strpos("`var'" , "crowded") label var `var' "Crowding"
if strpos("`var'", "_pd_") & strpos("`var'", "_2")==0     label var `var' "Proportion of all '`:variable label `var''' in district"
if strpos("`var'", "_pd_") & strpos("`var'", "_2")!=0     label var `var' "`:variable label `var'' (count)"
if strpos("`var'", "_e_")      label var `var' "`:variable label `var'' (%)"
if strpos("`var'", "_se")      label var `var' "`:variable label `var'' (Standard error)"
}


label var tpd_b_2 "Population in district" 
label var tpd_b "Propotion of sample in district" 
label var tpd_se "Propotion of sample in district (Standard error)" 
label var area "District area in square km"
label var state "State"
label var district "District"
sort state district
save all_estimates_clean.dta, replace 

drop hv024 shdistri
ds, has(type numeric)
foreach var of varlist `r(varlist)' {
tostring `var', replace force
}
set obs `=_N+1'
foreach var of varlist * {
replace `var' = "`:variable label `var''" if _n==_N
}
gen nr = _n
replace nr = -1 if  _N==_n
sort nr
drop nr
export delimited Table_S1_all_district_estimates , replace novarnames 


********************************************************************************
/*** Creating a table showing number of observations and missing observations 
for each variable, total and by district (Supplementary Table S2) */


use raw_estimates.dta, clear
drop area  *_ul *_ll *_eform *_crit *_t *_pvalue *_df hv024 shdistri v* *_b *_se l
order district state obst obst_hh crowded* nowash* over65* obsf* hyper* diab* sobese*

ds, has(type numeric)
sort state district
gen nr = _n
set obs `=_N+1'
replace nr = -1 if nr == .

foreach var of varlist `r(varlist)' {
sum `var' 
replace `var' = `r(sum)' if `var' == .
gen `var's = string(`var', "%9.0fc") 
drop `var'
}

replace district = "Total" if nr == -1
set obs `=_N+1'
replace nr = -2 if nr == .

foreach var of varlist * {
if strpos("`var'" , "obsf")    replace `var'  = "Females 15-49" if nr == -2
if strpos("`var'" , "obsts")   replace  `var' = "Observations" if nr == -2
if strpos("`var'" , "obst_")   replace  `var' = "Households" if nr == -2

if strpos("`var'" , "over65")  replace  `var'= "Elderly"  if nr == -2
if strpos("`var'" , "nowash")  replace  `var' ="No handwashing" if nr == -2
if strpos("`var'" , "hyper")   replace  `var' ="Hypertension" if nr == -2
if strpos("`var'" , "diab")    replace  `var' ="Diabetes" if nr == -2
if strpos("`var'" , "sobese")  replace  `var' ="Obesity" if nr == -2
if strpos("`var'" , "crowded") replace  `var' ="Crowding" if nr == -2



}
replace state ="State" if nr == -2
replace district ="District" if nr == -2
set obs `=_N+1'
replace nr = -3 if nr == .
foreach var of varlist crowded* nowash* over65* hyper* diab* sobese* {
replace `var' = "Missing" if nr == -3
}

sort nr
drop nr 
export delimited Table_S2_observations_and_missingness, replace novarnames 


********************************************************************************
********************************************************************************
********************************************************************************
/*** Making histograms for district estimates */



use all_estimates_clean.dta, clear
set type double
replace density = log10(density)
label var density "log{sub:10} Population density (people/square km)"

local col_density teal
local col_crowded_e_b maroon
local col_nowash_e_b sand
local col_over65_e_b purple
local col_hyper_e_b khaki
local col_diab_e_b sienna
local col_sobese_e_b blue

foreach var in density crowded_e_b nowash_e_b over65_e_b hyper_e_b diab_e_b sobese_e_b {
local lab = subinstr(subinstr(subinstr(subinstr("`: variable label `var''", " (%)", "",.), "log{sub:10}", "log10",.), " (people/square km)" , "",.), " ", "_",.)
sum `var'
local w = (`r(max)'-`r(min)')/round(sqrt(r(N)))+0.00000001 // this will end up being the same as the default stata bin width
local v = string(round(`w', 0.01), "%9.2fc")
local hist
local min = `r(min)'-0.0000000001 // otherwise there is a small error
egen double temp = cut(`var') , at(`min'(`w')101) icodes
local j = 0.1
distinct temp
local d = 1.4/`r(ndistinct)'
levelsof temp, local(levels)
foreach i in  `levels' {
local j = `j' + `d'
local hist  (hist `var' if temp == `i' , fcol(`col_`var''*`j') lcol(`col_`var''*`j') lwidth(vvthin ) freq width(`w') start(`min')) `hist'
}
twoway /* (hist `var', freq fcol(red) lcol(red) lwidth(vvthin )) */   ///
`hist' ///
/* (hist `var',   freq fcol(none) lcol(red) lwidth(vvthin ))  */ ///
 , graphregion(color("white") ///
margin (2 2 2 2)) legend(off)   xlab(#10, labsize(small)) ///
xtitle(, height(5) size(medsmall))  ytitle("Number of districts",  height(5) size(medsmall))  ///
ylab(#10, labsize(small) nogrid angle(0)) b2("Bin width: `v'",  col(gs6) size(2.2) ring(6) pos(6))
graph display, ysize(4) xsize(4)
graph export histogram_`lab'.pdf , replace
drop temp
}



********************************************************************************
********************************************************************************
********************************************************************************
/*** Making PC estimates */

/* A file titled "District to PC translator.csv" should be stored in the project 
folder */

use all_estimates_clean.dta , replace
gen nrm = _n
duplicates drop district state, force
keep district nrm state
gen dist = lower(district)
gen st = lower(state)
replace st =  "andhra pradesh" if st == "telangana"
save my_dist.dta, replace

import delimited "District to PC translator.csv", case(preserve) encoding(UTF-8) clear
gen nrp = _n
duplicates drop GIS_ST_NAME DISTRICT, force
drop if GIS_PC_NAME == " "
keep GIS_ST_NAME nrp DISTRICT
gen dist = lower(DISTRICT)
gen st =   lower(GIS_ST_NAME)
reclink dist st using my_dist.dta , idm(nrp) idu(nrm) gen(score) minscore(0)
sort score
keep state district DISTRICT GIS_ST_NAME
save PCdistkey.dta, replace
import delimited "District to PC translator.csv", case(preserve) encoding(UTF-8) clear
merge m:1 DISTRICT GIS_ST_NAME using PCdistkey.dta, nogen
merge m:m state district using all_estimates_clean.dta
drop if _merge == 1
drop _merge

foreach var in  crowded nowash  over65 {
bys ST_CODE PC_CODE: egen num = sum(`var'_e_b*Tot_pop*Pct_pop)
bys ST_CODE PC_CODE: egen den = sum(POP)
gen _`var'_b = num/den
label var _`var'_b "`:variable label `var'_e_b'"
drop num den
}

duplicates drop ST_CODE PC_CODE, force
keep _* GIS_PC_NAME state
save PC_estimates.dta, replace
sort state GIS_PC_NAME
drop if GIS_PC_NAME == ""

gen nr = _n
set obs `=_N+1'
replace nr = 0 if nr ==.
foreach var of varlist _* {
gen `var's = string(`var')
replace `var's = "`:variable label `var''" if nr == 0
drop `var'
}
replace GIS_PC_NAME = "Parliamentary Constituency" if nr == 0
replace state = "State" if nr == 0
sort nr
drop nr

export delimited Table_S5_PC_estimates , replace novarnames 


}
