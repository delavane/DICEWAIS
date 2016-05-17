$title Modeling WAIS collapse threat in DICE as a stochastic program with endogenous uncertainty
* Edited from the DICE-2013R model, version DICE2013Rv2_102213_vanilla_v24b.gms
* Delavane Diaz, delavane@stanford.edu
* Last revised Nov 24, 2015
$eolcom #

$if not set run $set run test
$if not set iteration $set iteration 1

* sensitivity cases
$if not set prtp $set prtp 0.015         # Set rate of social time preference
$if not set slrDF $set slrDF 0.0082      # Set 1m SLR damages
$if not set ECS $set ECS 2.9             # Set ECS
$if not set WAISrate $set WAISrate 10    # Set WAISrate
$if not set Ttrigger $set Ttrigger 4     # Set Ttrigger
$if not set tax $set tax 0               # Set Tax
$if not set model $set model collapse    # Set case: neglect_collapse, collapse
$if not set damages $set damages DICE    # Set SLR damage function per DICE or CIAM
$if not set HR $set HR Markov            # Set hazard rate formulation as Markov or Bayes (note Bayes creates discontinuities)

* DICE2013, modified to run in 10 year time steps (vs 5 yr) with 2010 as t=1
$if not set yr $set yr 30 # Set total time horizon
sets sequence   /0/
ct              total climate time periods /1*%yr%/
t(ct)           time periods /1*%yr%/
collapset       collapse time periods /1*10/
tfirst(t)       first time period
;
alias(t,tt);
tfirst(t) = yes$(t.val eq 1);

parameters
tstep           Years per Period /10/
** Preferences
elasmu          Elasticity of marginal utility of consumption (1.45 default computed below)
prstp           Initial rate of social time preference per year /%prtp%/
** Population and technology
gama            Capital elasticity in production function/.300/
pop0            Initial world population (millions) /6838/
popadj          Growth rate to calibrate to 2050 pop projection /0.134/
popasym         Asymptotic population (millions) /10500/
dk              Depreciation rate on capital (per year) /.100/
q0              Initial world gross output (trill 2005 USD) /63.69/
k0              Initial capital value (trill 2005 USD) /135/
a0              Initial level of total factor productivity /3.80/
ga0             Initial growth rate for TFP per period (0.079 per 5 yr --> 1.079^2-1=0.164241) /0.164/
dela            Decline rate of TFP per year /0.006/
l(ct)           Level of population and labor
al(t)           Level of total factor productivity
sigma(t)        CO2-equivalent-emissions output ratio
rr(ct)          Average utility social discount rate
ga(t)           Growth rate of productivity from
gl(t)           Growth rate of labor
gfacpop(t)      Growth factor population
** Emissions parameters
gsigma1         Initial growth of sigma (per year) /-0.01/
dsig            Decline rate of decarbonization (per YEAR - typo in code? period) /-0.001/
eland0          Carbon emissions from land 2010 (GtCO2 per year) /3.3/
deland          Decline rate of land emissions (per period) /.2/
e0              Industrial emissions 2010 (GtCO2 per year) /33.61/
miu0            Initial emissions control rate for base case 2010 /0 /  # DICE2013 0.039
forcoth(t)      Exogenous forcing for other greenhouse gases
etree(t)        Emissions from deforestation
** Carbon cycle
* Initial Conditions
mat0            Initial Concentration in atmosphere 2010 (GtC)/830.4/
mu0             Initial Concentration in upper strata 2010 (GtC) /1527./
ml0             Initial Concentration in lower strata 2010 (GtC) /10010./
mateq           Equilibrium concentration atmosphere (GtC) /588/
mueq            Equilibrium concentration in upper strata (GtC) /1350/
mleq            Equilibrium concentration in lower strata (GtC) /10000/
* Flow paramaters
* Use different parameters to reflect the 10 year time step -- percent per period
b12             Carbon cycle transition matrix - fraction of conc that goes into upper ocean (0.088 per 5 yr) /0.176/
b23             Carbon cycle transition matrix (0.00250 per 5 yr) /0.005/
* These are for declaration and are defined later
b11             Carbon cycle transition matrix
b21              Carbon cycle transition matrix
b22             Carbon cycle transition matrix
b32             Carbon cycle transition matrix
b33             Carbon cycle transition matrix
sig0             Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)
** Climate model parameters
*Use different parameters to reflect the 10 year time step
t2xco2          Equilibrium temp impact (oC per doubling CO2) /2.9/
fex0            2010 forcings of non-CO2 GHG (Wm-2) /0.25/
fex1            2100 forcings of non-CO2 GHG (Wm-2) /0.70/
tocean0         Initial lower stratum temp change (C from 1900) /.0068/
tatm0           Initial atmospheric temp change (C from 1900) /0.80/
c10             Initial climate equation coefficient for upper level (0.098 per 5 yr --> keep this as 5 yr and adjust below) /0.098/
c1beta          Regression slope coefficient(SoA~Equil TSC) /0.01243/
c1              Climate equation coefficient for upper level (0.098 per 5 yr) /0.196/ # Speed of adjustment
c3              Transfer coefficient upper to lower stratum (0.088 in DICE2013 - time independent) /0.088/
c4              Transfer coefficient for lower level (0.0250 per 5 yr) /0.05/
fco22x          Forcings of equilibrium CO2 doubling (Wm-2) /3.8/
*slr0           Initial sea level rise in 2005 per DICE2010 (m) /0.11/
slr0            Initial sea level rise since 2000 (m) in 2010 per Church and White 2011 /0.04/
Tfreeze         freezing temperature of seawater (deg C) per Shaffer 2014 /-1.8/
Ta0             mean Antarctic temperature reduced to sea level for 1961-1990 (deg C) per Shaffer 2014 /-18/
paf             polar amplification factor for Southern Hemisphere per Masson-Delmotte 2005 /1.2/
Tadj0           temperature above pre-industrial for 1961-1990 (deg C) per NCDC Global Surface Temperature Anomalies /0.29/
To0             compute pre-industrial ocean subsurface temperature adjacent to the AIS (deg C)
lam             Climate model parameter
waisvol         WAIS volume in SLR-equivalent (m) /5/
waisrate        WAIS collapserate (m per year)
beta            Decadal hazard rate coefficient as a Markovian function
WAISthresholdT  WAIS threshold
** Climate damage parameters per DICE 2010 which splits SLR and Temp damages
a1              Damage coefficient on temperature /0.00008162/
a2              Damage coefficient on temperature squared/0.00204626/
a3              Exponent on temperature damages /2/
b1              Damage coefficient on SLR/0.00518162/
b2              Damage coefficient on SLR squared/0.00305776/
b3              Exponent on SLR damages /2/
** Abatement cost
expcost2        Exponent of control cost function /2.8/
pback           Cost of backstop 2005$ per tCO2 2010 /344/
gback           Initial cost decline backstop cost per period /.025/
limmiu(t)       Upper limit on control rate (allow 20% negative emissions after 2100)
tnopol          Period before which no emissions controls base /45/
cprice0         Initial base carbon price (2005$ per tCO2) /1.0/
gcprice         Growth rate of base carbon price per year /.02/
gsig(t)         Change in sigma (cumulative improvement of energy efficiency)
cost1(t)        Adjusted cost for backstop
pbacktime(t)    Backstop price
optlrsav        Optimal long-run savings rate used for transversality
cpricebase(t)   Carbon price in base case
** Availability of fossil fuels
fosslim         Maximum cumulative extraction fossil fuels (GtC) /6000/
;

* Read in uncertain parameters and make any adjustments
parameters inputdata, reltarget, tax;
$gdxin DICEWAISinputs
$load inputdata
t2xco2 = inputdata('%iteration%','ECS');
waisrate = inputdata('%iteration%','WAISrate')/1000; # Convert WAIS collapserate from mm/yr to m/yr (e.g., 0.01 m per year Bamber 2013)
beta = inputdata('%iteration%','hazardbeta');
tax = %tax%;

*Transient TSC Correction ("Speed of Adjustment Parameter")
c1 = 2*(c10 + c1beta*(t2xco2-2.9)); # 2x --> to reflect 10 yr periods
* Maintain same linear and quadratic ratio in SLR DF
b1=%slrDF%*0.00518162/0.00823938;
b2=%slrDF%*0.00305776/0.00823938;
* adjust consumption elasticity when we alter prstp: r=elasmu*g+prstp with r=0.04623, g=0.02154
elasmu = (0.04623-%prtp%)/0.02154;

$if %scenario% == tax a1=0; a2=0; b1=0; b2=0;

* Parameters for long-run consistency of carbon cycle
b11 = 1 - b12;
b21 = b12*MATEQ/MUEQ;
b22 = 1 - b21 - b23;
b32 = b23*mueq/mleq;
b33 = 1 - b32 ;
* Further definitions of parameters
sig0 = e0/(q0*(1-miu0));
lam = fco22x/t2xco2;
l("1") = pop0;
loop(t, l(t+1)=l(t););
loop(t, l(t+1)=l(t)*(popasym/L(t))**popadj ;);
loop(ct$(not t(ct)), l(ct)=l(ct-1); ); # hold labor constant after end of time horizon
ga(t)=ga0*exp(-dela*tstep*((t.val-1)));
al("1") = a0; loop(t, al(t+1)=al(t)/((1-ga(t))););
gsig("1")=gsigma1; loop(t,gsig(t+1)=gsig(t)*((1+dsig)**tstep) ;);
sigma("1")=sig0; loop(t,sigma(t+1)=(sigma(t)*exp(gsig(t)*tstep)););
pbacktime(t)=pback*(1-gback)**(t.val-1);
cost1(t) = pbacktime(t)*sigma(t)/expcost2/1000;
etree(t) = eland0*(1-deland)**(t.val-1);
rr(ct) = 1/((1+prstp)**(tstep*(ct.val-1)));
forcoth(t) = fex0+ (1/18)*(fex1-fex0)*(t.val-1)$(t.val lt 19)+ (fex1-fex0)$(t.val ge 19);
optlrsav = (dk + .004)/(dk + .004*elasmu + prstp)*gama;
*Base Case Carbon Price
cpricebase(t)= cprice0*(1+gcprice)**(tstep*(t.val-1));
limmiu(t)=1; # except cases below:
limmiu('1')=0.1;
limmiu(t)$(t.val>11)=1.2;



** This section corresponds to the model case specified earlier (e.g., neglect_collapse, collapse, collapse%certain%, ev_collapse
$goto %model%


$label neglect_collapse
set s           states (year of collapse) /0/
as(s,t)         active states;
parameter
o(s,t)          state offset pointer
collapse(s,t)   indicator for collapse trigger (0 or 1)
;
collapse(s,t) = 0;
as(s,t) = yes;
o(s,t) = 0;
$goto model


$label collapse
set
s states (year of collapse) /set.sequence,set.collapset/ # collapse horizon can be less than the total horizon
as(s,t)         active states;
parameter
o(s,t)          state offset pointer
collapse(s,t)   indicator for collapse trigger (0 or 1)
$if set SP pr(s) state probability # note probabilities will be fixed exogenously
;
collapse(s,t) = 0;
loop(s$(not sameas(s,'0')),
 o(s,t)$(t.val<=s.val) = -s.val;
 as(s,t) = yes$(o(s,t) = 0);
* If a collapse occurs, increased slr rate does not apply until start of NEXT time period
 collapse(s,t)$(t.val > s.val) = 1;
);
o('0',t) = 0;
as('0',t) = yes;
$goto model


** Main model code continues here
$label model
$macro sw(s,t) s+o(s,t),t

VARIABLES
MIU(s,t)        Emission control rate GHGs
FORC(s,t)       Increase in radiative forcing (watts per m2 from 1900)
TATM(s,t)       Increase temperature of atmosphere (degrees C from 1900)
TOCEAN(s,t)     Increase temperatureof lower oceans (degrees C from 1900)
MAT(s,t)        Carbon concentration increase in atmosphere (GtC from 1750)
MU(s,t)         Carbon concentration increase in shallow oceans (GtC from 1750)
ML(s,t)         Carbon concentration increase in lower oceans (GtC from 1750)
E(s,t)          Total CO2 emissions (GtCO2 per year)
EIND(s,t)       Industrial emissions (GtCO2 per year)
C(s,t)          Consumption (trillions 2005 US dollars per year)
K(s,t)          Capital stock (trillions 2005 US dollars)
CPC(s,t)        Per capita consumption (thousands 2005 USD per year)
I(s,t)          Investment (trillions 2005 USD per year)
SR(s,t)         Gross savings rate as fraction of gross world product
RI(s,t)         Real interest rate (per annum)
Y(s,t)          Gross world product net of abatement and damages (trillions 2005 USD per year)
YGROSS(s,t)     Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
DAMAGES(s,t)    Damages (trillions 2005 USD per year)
DAMFRAC(s,t)    Damages as fraction of gross output
ABATECOST(s,t)  Cost of emissions reductions (trillions 2005 USD per year)
MCABATE(s,t)    Marginal cost of abatement (2005$ per ton CO2)
CCA(s,t)        Cumulative industrial carbon emissions (GTC)
PERIODU(s,t)    One period utility function
CPRICE(s,t)     Carbon price (2005$ per ton of CO2)
CEMUTOTPER(s,t) Period utility
UTILITY(s)      Welfare function
EU              Expected Utility
* Define states and probabilities
PR(s)            state probability
HR(t)            collapse hazard rate
PI(t)            probability of no collapse at start of time t
* Coastal impacts
WAIS(s,t)       SLR component from WAIS collapse in terms of m per period
SLR(s,t)        total sea level rise in year t under state scenario s (m)
COASTDAMAGES(s,t) Coastal damages ($T - note calibrated as $B so convert)
;
NONNEGATIVE VARIABLES MIU, TATM, MAT, MU, ML, Y, YGROSS, C, K, I, SLR, DAMAGES, COASTDAMAGES, PR, HR, PI;

EQUATIONS
*Emissions and Damages
EEQ(s,t)        Emissions equation
EINDEQ(s,t)     Industrial emissions
CCACCA(s,t)     Cumulative carbon emissions
FORCE(s,t)      Radiative forcing equation
DAMFRACEQ(s,t)  Equation for damage fraction
DAMEQ(s,t)      Damage equation
ABATEEQ(s,t)    Cost of emissions reductions equation
MCABATEEQ(s,t)  Equation for MC abatement
CARBPRICEEQ(s,t) Carbon price equation from abatement
*Climate and carbon cycle
MMAT(s,t)       Atmospheric concentration equation
MMU(s,t)        Shallow ocean concentration
MML(s,t)        Lower ocean concentration
TATMEQ(s,t)     Temperature-climate equation for atmosphere
TOCEANEQ(s,t)   Temperature-climate equation for lower oceans
*Economic variables
YGROSSEQ(s,t)   Output gross equation
YY(s,t)         Output net equation
CC(s,t)         Consumption equation
CPCE(s,t)       Per capita consumption definition
SEQ(s,t)        Savings rate equation
KK(s,t)         Capital balance equation
RIEQ(s,t)       Interest rate equation
* Utility
CEMUTOTPEREQ(s,t) Period utility
PERIODUEQ(s,t) Instantaneous utility function equation
UTIL(s)         Objective function
WELFARE         Objective function
* Define states and probabilities
HRDEF(t)         Defines HR
PIDEF(t)         Defines PI
PRDEF(s,t)       Defines PR
PR0(s)           Ensure no collapse scenario accounted for
* Coastal impacts
WAISDEF(s,t)    Define SLR component from WAIS collapse in terms of m per period
SLRDEF(s,t)     Define total sea level rise in year t under state s (m)
COASTDAMAGESDEF(s,t) Define coastal damages ($T)
;

$if set SP $goto DICE

** Characterize WAIS collapse hazard
set tp(t)       Time periods with tipping points;
tp(t) = yes$(t.val <= smax(s,s.val));

* Hazard rate depends on the temperature in s=0 world
* beta has been calibrated as a decadal hazard rate of warming since 2000
HRDEF(tp(t))$(t.val>1).. HR(t) =e=
$if %HR%=='Markov' beta*(TATM('0',t)-0.6)**2
;
* Hazard rate parameters beta and alpha are calibrated offline in DICEWAIS_uncertain_parameters.R (previously collapsecalibrationv3.gms)
HR.FX('1')=0;

* Probability that collapse has not occurred at start of time period t (this is only enforced for periods that have tipping potential)
PIDEF(tp(t)).. PI(t) =e=
prod(tt$(ord(tt) lt ord(t)), 1-HR(tt) )
;
* Probability that a tipping point occurs in time period t:
PRDEF(s,t)$sameas(s,t).. PR(s) =e= HR(t) * PI(t);
PR0(s)$(s.val=0).. PR(s) =e= 1 - sum(tp(t), HR(t) * PI(t));


** Equations of the DICE 2013 model
$label DICE
*Emissions and Damages
eeq(as(s,t)).. E(s,t) =E= EIND(s,t) + etree(t); # annual emissions
eindeq(as(s,t)).. EIND(s,t)=E= sigma(t) * YGROSS(s,t) * (1-(MIU(s,t)));
ccacca(as(s,t+1))..CCA(s,t+1) =E= CCA(sw(s,t))+ EIND(sw(s,t))*tstep/3.666;
force(as(s,t)).. FORC(s,t)=E= fco22x * ((log((MAT(s,t)/588.000))/log(2))) + forcoth(t);
damfraceq(as(s,t)).. DAMFRAC(s,t) =E= a1*TATM(s,t) + a2*TATM(s,t)**a3;
dameq(as(s,t)).. DAMAGES(s,t) =E= YGROSS(s,t) * DAMFRAC(s,t);
abateeq(as(s,t)).. ABATECOST(s,t) =E= YGROSS(s,t) * cost1(t) * (MIU(s,t)**expcost2);
mcabateeq(as(s,t)).. MCABATE(s,t) =E= pbacktime(t) * MIU(s,t)**(expcost2-1);
carbpriceeq(as(s,t)).. CPRICE(s,t) =E= pbacktime(t) * MIU(s,t)**(expcost2-1);

*Climate and carbon cycle
*All of these equations now use different parameters to reflect the 10 year time step
mmat(as(s,t+1)).. MAT(s,t+1) =E= MAT(sw(s,t))*b11 + MU(sw(s,t))*b21 + (E(sw(s,t))*(tstep/3.666));
mml(as(s,t+1)).. ML(s,t+1)=E= ML(sw(s,t))*b33 + MU(sw(s,t))*b23;
mmu(as(s,t+1)).. MU(s,t+1)=E= MAT(sw(s,t))*b12 + MU(sw(s,t))*b22 + ML(sw(s,t))*b32;
tatmeq(as(s,t+1))..TATM(s,t+1) =E= TATM(sw(s,t)) + c1 * ((FORC(s,t+1)-(fco22x/t2xco2)*TATM(sw(s,t)))-(c3*(TATM(sw(s,t))-TOCEAN(sw(s,t)))));
toceaneq(as(s,t+1)).. TOCEAN(s,t+1) =E= TOCEAN(sw(s,t)) + c4*(TATM(sw(s,t))-TOCEAN(sw(s,t)));

*Economic variables
ygrosseq(as(s,t))..YGROSS(s,t) =E= (al(t)*(l(t)/1000)**(1-gama))*(K(s,t)**gama);
yy(as(s,t)).. Y(s,t) =E= YGROSS(s,t) - ABATECOST(s,t) - DAMAGES(s,t) - COASTDAMAGES(s,t);
cc(as(s,t)).. C(s,t) =E= Y(s,t) - I(s,t);
cpce(as(s,t)).. CPC(s,t) =E= 1000 * C(s,t) / l(t);
seq(as(s,t)).. I(s,t) =E= SR(s,t) * Y(s,t);
kk(as(s,t+1)).. K(s,t+1) =L= (1-dk)**tstep * K(sw(s,t)) + tstep * I(sw(s,t));
rieq(as(s,t+1)).. RI(s,t) =E= (1+prstp) * (CPC(s,t+1)/CPC(sw(s,t)))**(elasmu/tstep) - 1;
*Utility
periodueq(as(s,t)).. PERIODU(s,t) =E= ((C(s,t)*1000/l(t))**(1-elasmu)-1)/(1-elasmu);
cemutotpereq(as(s,t)).. CEMUTOTPER(s,t) =E= PERIODU(s,t) * l(t) * rr(t);
util(s).. UTILITY(s)=E= tstep * [sum(t, CEMUTOTPER(sw(s,t)))];
WELFARE.. EU =e= 1e-3*sum(s, PR(s) * UTILITY(s) );

** additional equations related to SLR
* Compute SLR component from WAIS discharge (m per period)
WAISDEF(as(s,t+1)).. WAIS(s,t) =E=
$if set montecarlo              inputdata('%iteration%','WAISbaserate')/1000*tstep +
$if not set montecarlo          0.2833/1000*tstep +              # mm/yr per Sheperd et al 2012
                                   collapse(sw(s,t))*            # this turns collapse on/off
                                         waisrate*tstep          # m/period, constant discharge
;

* "Baseline" SLR per DICE2010 - thermal expansion, GSIC, GIS component (not WAIS) in decadal rates as function of temp; add WAIS
SLRDEF(as(s,t+1)).. SLR(s,t+1) =E= (0.00779+0.0314*0.26+0.00176*7.3)*TATM(sw(s,t)) + SLR(sw(s,t)) + WAIS(sw(s,t)); # rate terms in m/period

** SLR Damages
* DICE2010 SLR damage function
$if %damages%==DICE     COASTDAMAGESDEF(as(s,t)).. COASTDAMAGES(s,t) =E= YGROSS(s,t)*(1-1/(1+b1*SLR(s,t)+b2*SLR(s,t)**b3));
* CIAM least cost SLR damage function (Note we convert CIAM damage function from $B to $T to be consistent with DICE)
$if not %damages%==DICE COASTDAMAGESDEF(as(s,t)).. COASTDAMAGES(s,t) =E= YGROSS(s,t)*[0.0009624531*SLR(s,t)];



*Limits
CCA.up(s,t) = fosslim;
MIU.up(s,t) = limmiu(t);
SR.up(s,t) = 1;
** Upper and lower bounds for stability
K.LO(s,t) = 1;
MAT.LO(s,t) = 10;
MU.LO(s,t)= 100;
ML.LO(s,t)= 1000;
C.LO(s,t) = 2;
TOCEAN.UP(s,t) = 20;
TOCEAN.LO(s,t) = -1;
TATM.UP(s,t) = 10;
TATM.LO(s,t) = 0.6;
CPC.LO(s,t) = .01;

* Initial conditions
CCA.FX(s,tfirst) = 90;
K.FX(s,tfirst) = k0;
MAT.FX(s,tfirst) = mat0;
MU.FX(s,tfirst) = mu0;
ML.FX(s,tfirst) = ml0;
TATM.FX(s,tfirst) = tatm0;
TOCEAN.FX(s,tfirst) = tocean0;

* DICE-WAIS additions
PI.L(t)$(ord(t)>1) = 0.1;
PI.FX(tfirst) = 1; PI.FX('2') = 1;
PR.L(s)$(ord(s)>1) = 0.1;
MIU.L(s,t) = min(t.val/10,1);
TATM.L(s,t) = tatm0+t.val/10;
MIU.FX(as(s,tfirst)) = 0; # No mitigation has happened in first period
SLR.FX(s,tfirst) = slr0;

CEMUTOTPER.L(as(s,t))=1000;
SR.FX(as(s,t))$(ord(t)>1) = optlrsav;

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = off;
option limrow = 0;
option limcol = 0;
model DICESLR /all/;

$if %scenario% == tax CCA.up(s,t) = inf; MIU.l(s,t)=1; cprice.fx(s,t)$((ord(t)>1) and (ord(t)<10)) = min(pbacktime(t), tax*(1.04**(tstep*(t.val-2))) );


solve DICESLR maximizing EU using nlp;
solve DICESLR maximizing EU using nlp;
