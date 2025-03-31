* Stata codes for figures and endogenous kink model publication bias test

* Import data 
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


* Codes for funnel plots and box plots: asymmetric estimates

twoway (scatter inv_se beta, mcolor(dknavy) msymbol(circle_hollow)) if hike == 1 & inter_overnight == 1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% overnight rate hike"') xline(0.600, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.374, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Upward pass-through: overnight interbank rate"') xsize(6.5) ysize(5.0)


twoway (scatter inv_se beta, mcolor(dknavy) msymbol(circle_hollow)) if cut == 1 & inter_overnight == 1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% overnight rate cut"') xline(0.376, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.235, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Downward pass-through: overnight interbank rate"') xsize(6.5) ysize(5.0)



twoway (scatter inv_se beta, mcolor(dknavy) msymbol(circle_hollow)) if hike == 1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% reference rate hike"') xline(0.443, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.348, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Upward pass-through: full sample"') xsize(6.5) ysize(5.0)


twoway (scatter inv_se beta, mcolor(dknavy) msymbol(circle_hollow)) if cut == 1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% reference rate cut"') xline(0.211, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.21, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Downward pass-through: full sample"') xsize(6.5) ysize(5.0)





graph hbox beta, over(Study) xsize(9) ysize(11) scale(0.4) yline(0.435, lpattern(solid) lcolor(red) lwidth(vthin))


graph hbox beta if inter_overnight == 1, over(Study) xsize(9) ysize(11) scale(0.5) yline(0.375, lpattern(solid) lcolor(red) lwidth(vthin))




* Declare panel dataset
xtset study_id



*** Overnight interbank reference rate sample

** OLS model

* Upward pass-through:
regress t inv_se if hike == 1 & inter_overnight == 1, vce(cluster study_id) 

* Downward pass-through:
regress t inv_se if cut == 1 & inter_overnight == 1, vce(cluster study_id) 



** Fixed effects model

* Upward pass-through:
xtreg t inv_se if hike == 1 & inter_overnight == 1, fe vce(cluster study_id)

* Downward pass-through:
xtreg t inv_se if cut == 1 & inter_overnight == 1, fe vce(cluster study_id)



** Random effects model

* Upward pass-through:
xtreg t inv_se if hike == 1 & inter_overnight == 1, mle vce(cluster study_id)

* Downward pass-through:
xtreg t inv_se if cut == 1 & inter_overnight == 1, mle vce(cluster study_id)

 



