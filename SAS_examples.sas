**************************************************************;
***** PROGRAM FOR: "An Introduction to G-Methods"
***** AUTHORS: ASHLEY I NAIMI (ashley.naimi@pitt.edu)
*****          EDWARD H KENNEDY 
*****          STEPHEN R COLE
***** VERSION: 3 
***** PURPOSE: Illustrate all three g-methods
**************************************************************;
*****the data;
data a;
	input a0 z1 a1 y n;
	datalines;
	0 0 0  87.288 209271
	0 0 1 112.107 93779
	0 1 0 119.654 60657
	0 1 1 144.842 136293
	1 0 0 105.282 134781
	1 0 1 130.184 60789
	1 1 0 137.720 93903
	1 1 1 162.832 210527
	;
run;

*SECTION 1;
*STANDARD REGRESSION MODELS IN TABLE 1;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION;
title "Crude Model";
freq n;
model y = a0 a1 a0*a1;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION;
title "Z-Adjusted Model";
freq n;
model y = a0 a1 a0*a1 z1;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION;
title "A0 Only Model";
freq n;
model y = a0;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION;
title "Z1-Adjusted A0 Only Model";
freq n;
model y = a0 z1;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION;
title "A1 Only Model";
freq n;
model y = a1;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION;
title "Z1-Adjusted A1 Only Model";
freq n;
model y = a1 z1;
run;quit;run;

*MARGINAL STRUCTURAL MODEL;
ods select none;
proc logistic data=a desc noprint; 
	model a0=; 
	freq n; 
	output out=b0 p=p_num_a0; 
proc logistic data=a desc noprint; 
	model a1=; 
	freq n; 
	output out=b1 p=p_num_a1; 
proc logistic data=a desc noprint; 
	model a1= z1; 
	freq n; 
	output out=b2 p=p_den_a1; 
data c; 
	merge a b0 b1 b2;
	*STABILIZED WEIGHTS;
	*NB p_num_a0 = p_den_a0, so we use former for both;
	if a0=1 and a1=1 then sw=(p_num_a0/p_num_a0)*(p_num_a1/p_den_a1); 
	if a0=1 and a1=0 then sw=(p_num_a0/p_num_a0)*((1-p_num_a1)/(1-p_den_a1)); 
	if a0=0 and a1=1 then sw=((1-p_num_a0)/(1-p_num_a0))*(p_num_a1/p_den_a1); 
	else if a0=0 and a1=0 then sw=((1-p_num_a0)/(1-p_num_a0))*((1-p_num_a1)/(1-p_den_a1));
	pseudoN = sw*n;
	drop _level_;
	label d0=" " d1=" " d2=" ";
run;quit;run;
proc means data=c min max mean sum maxdec=2;
ods select all;
title "IPW Distribution";
freq n;var sw pseudoN;run;
proc print data=c;
title "Pseudo Data: weights and probabilities";
sum pseudoN;run;
proc genmod data=c;*NB: NOT CONSIDERING STANDARD ERROR ESTIMATION WHICH WOULD REQUIRE USE OF REPEATED STATEMENT;
title "Marginal Structural Model";
freq n;weight sw;
model y = a0 a1 a0*a1;
ods output ParameterEstimates=psi_msm(where=(parameter="a0"|parameter="a1"|parameter="a0*a1") keep=parameter estimate);
run;quit;run;
title;

data _null_;
	set psi_msm;
	if parameter="a0" then call symput('psi0_msm',estimate);
	if parameter="a1" then call symput('psi1_msm',estimate);
	if parameter="a0*a1" then call symput('psi2_msm',estimate);
run;
%put &psi0_msm;%put &psi1_msm;%put &psi2_msm;

*G-FORMULA;
*Transform to one observation per individual;
data aa;set a;do i = 1 to n;output;end;drop n i;run; 

*STEP 2: FIT A MODEL FOR EACH COMPONENT OF THE LIKELIHOOD;
ods select none;
proc genmod data=aa; model y = a1 z1 a0 a1*a0;
proc logistic data=aa desc; model a1 = z1; 
proc logistic data=aa desc; model z1 = a0; 
run;quit;run;

*STEP 3: SELECT RANDOM SAMPLE OF OBSERVED DATA AND KEEP ONLY BASELINE COVARIATES (a0);
proc surveyselect data=aa out=mc0(keep= a0 NumberHits) method=urs seed=78419 sampsize=1000000;run;
data mc0;set mc0;do i = 1 to NumberHits;output;end;keep a0;run;

*STEP 4: G-FORMULA MACRO;
%macro gform(run=, a0=, a1=, z1=); 
data mc&run;
	set mc0;
	call streaminit(987);
	z1 = rand("bernoulli",1/(1+exp(-(-0.4309 + 0.8735*&a0))));
	a1 = rand("bernoulli",1/(1+exp(-(-0.8002 + 1.6084*&z1))));
	y = rand("normal",87.2466 + 24.9699*&a1 +  32.5502*&z1 + 17.9893*&a0 + 0.0541*&a1*&a0,1);
run;
%mend;
%gform(run=1,a0=a0,a1=a1,z1=z1);*natural course;
%gform(run=2,a0=1,a1=1,z1=z1);*fully exposed;
%gform(run=3,a0=1,a1=0,z1=z1);*time zero only;
%gform(run=4,a0=0,a1=1,z1=z1);*time one only;
%gform(run=5,a0=0,a1=0,z1=z1);*fully unexposed;

ods select all;
proc means data=aa min max mean  ;title "Observed Data";var a0 z1 a1 y;run;
proc means data=mc1 min max mean ;title "Natural Course";var a0 z1 a1 y;run;
proc means data=mc2 min max mean ;title "Fully Exposed";var a0 z1 a1 y;run;
proc means data=mc3 min max mean ;title "Time Zero Only";var a0 z1 a1 y;run;
proc means data=mc4 min max mean ;title "Time One Only";var a0 z1 a1 y;run;
proc means data=mc5 min max mean ;title "Fully Unexposed";var a0 z1 a1 y;run;
title;

*G-ESTIMATION OF A SNMM;
* INSTRUMENTAL VARIABLE ESTIMATOR;
proc logistic data=a desc;
freq n;
model a1 = z1 a0 a0*z1;
output out=aa2 pred=p1;run;
data aa2;set aa2;resid1=a1-p1;
proc syslin data=aa2 2sls;
	weight n;
	where a0=1;
	endogenous a1;
	instruments resid1;
	model y = a1;
	ods output ParameterEstimates=twosls1(where=(variable="a1"));
run;quit;run;
proc syslin data=aa2 2sls;
	weight n;
	where a0=0;
	endogenous a1;
	instruments resid1;
	model y = a1;
	ods output ParameterEstimates=twosls2(where=(variable="a1"));
run;quit;run;
data twosls1;set twosls1;rename estimate=psi1;keep estimate;
data twosls2;set twosls2;rename estimate=psi02;keep estimate;
data twosls;
	merge twosls1 twosls2;
	psi2=psi1-psi02;
	_LEVEL_=1;
	drop psi02;
proc print data=twosls;run;
data aa2;merge aa2 twosls;by _LEVEL_;
y_tilde = y - a1*psi1 - a0*a1*psi2;
run;quit;run;
proc logistic data=aa2 desc;
freq n;
model a0 = ;
output out=aa2 pred=p0;run;
data aa2;set aa2;resid0=a0-p0;
proc print data=aa2 (obs=20);run;
proc syslin data=aa2 2sls;
	weight n;
	endogenous a0;
	instruments resid0;
	model y_tilde = a0;
	ods output ParameterEstimates=twosls0(where=(variable="a0"));
run;quit;run;
data twosls0;set twosls0;_level_=1;rename estimate=psi0;keep estimate _level_;
data GEst_IV;
	merge twosls0 twosls;by _level_;
	drop _level_;
run;quit;run;
proc print data=GEst_IV noobs;
title "SNMM Effect Estimates: IV Estimator";
run;quit;run;

* OLS ESTIMATOR;
ods select none;
*fit exposure models;
proc logistic data=a desc;
model a1 = z1 a0 a0*z1;
freq n;
output out=aa2 pred=p1;run;
proc logistic data=aa2 desc;
model a0 = ;
freq n;
output out=aa3 pred=p0;
data aa4;
	set aa3;
	psi1 = (a1-p1);
	psi2 = a0*(a1-p1);
proc reg data=aa4;
*estimate parameters for A1 and A1*A0 interaction using no intercept OLS;
title "G-Estimation of a SNMM: A1 and A1*A0 Effects";
ods select ParameterEstimates;
model y = psi1 psi2  / noint;
freq n;
ods output ParameterEstimates=psi_ge(where=(variable="psi1"|variable="psi2") keep=variable estimate);
run;
proc transpose data=psi_ge out=psi00(keep=psi1 psi2);id variable;run;
proc print data=psi00;run;
data _null_;
	set psi_ge;
	if variable="psi1" then call symput('psi1_snm',estimate);
	if variable="psi2" then call symput('psi2_snm',estimate);
run;
%put &psi1_snm;%put &psi2_snm;
run;quit;run;
data aa4;
	set aa4;
	*subtract A1 and A1*A0 effect from observed outcome to get potential outcome if A1=0;
	y0 = (y - &psi1_snm*a1 - &psi2_snm*a1*a0);
	psi0 = (a0-p0);
proc reg data=aa4;
*estimate parameters for A0 effect;
title "G-Estimation of a SNMM: A0 Effect";
ods output ParameterEstimates=psi01(keep=estimate rename=(estimate=psi0));
model y0 = psi0 / noint ;
freq n;
proc print data=psi01;run;
data GEst_OLS;merge psi01 psi00;
proc print data=GEst_OLS round;
title "SNMM Effect Estimates: OLS Estimator";
run;quit;run;






