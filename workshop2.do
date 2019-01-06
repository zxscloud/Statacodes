
use "Workshop 2-STATA.dta",clear
rename _all,lower

net install st0026_2.pkg

pscore offfarmwork age education gender householdsize east central northwest,pscore(ps98) blockid(blockf1) comsup level(0.001)

**NN match
attnd pchhincome offfarmwork,pscore(ps98) comsup bootstrap reps(50)

**KBM
attk pchhincome offfarmwork,pscore(ps98) comsup

**RM
attr pchhincome offfarmwork,pscore(ps98) comsup radius(0.001)

********
****NNM
psmatch2 offfarmwork age education gender householdsize east central northwest,neighbor(7) bwidth(0.4) logit outcome(pchhincome) qui

***KBM
psmatch2 offfarmwork age education gender householdsize east central northwest,kernel bwidth(0.4) logit outcome(pchhincome) qui

***RM
psmatch2 offfarmwork age education gender householdsize east central northwest,radius bwidth(0.4) logit outcome(pchhincome) qui



*****RA

teffects ra (pchhincome age education gender householdsize east central northwest) (offfarmwork),ate

*****IPW

teffects ipw (pchhincome) (offfarmwork age education gender householdsize east central northwest,logit) ,atet


****IPWRA

teffects ipwra (pchhincome age education gender householdsize east central northwest) ///
	 (offfarmwork age education gender householdsize east central northwest,logit ),ate

teffects ,coeflegend

nlcom _b[ATE:r1vs0.offfarmwork]/_b[POmean:0.offfarmwork]





teffects ipwra (pchhincome age education gender householdsize east central northwest) ///
	 (offfarmwork age education gender householdsize east central northwest,logit ),atet
	 
teffects ,coeflegend

nlcom _b[ATET:r1vs0.offfarmwork]/_b[POmean:0.offfarmwork]
***AIPW

teffects aipw (pchhincome age education gender householdsize east central northwest) ///
	 (offfarmwork age education gender householdsize east central northwest,logit )



*****Parametric approaches

***TEM

treatreg pchhincome age education gender householdsize east central northwest, ///
	treat(offfarmwork = age education gender householdsize east central northwest distancetowork) robust

predict yc11,yctrt
predict yc10,ycntrt

sum yc11 yc10

gen ATE=yc11-yc10

ttest ATE==0



****ESR

movestay pchhincome age education gender householdsize east central northwest, ///
	select(smartphoneuse = age education gender householdsize east central northwest socialnetwork) robust

//ATT
predict yc11hi,yc1_1
predict yc21hi,yc2_1
sum yc11hi yc21hi
gen atthi=yc11hi-tc21hi
ttest atthi==0

//ATU

predict yc12hi,yc1_2
predict yc22hi,yc2_2
sum yc12hi yc22hi
gen atuhi=yc12hi-tc22hi
ttest atuhi==0



