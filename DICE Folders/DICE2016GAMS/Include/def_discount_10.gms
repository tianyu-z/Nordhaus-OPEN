* Discount program


elasmu   = .001;
prstp    = .01;
rr(t) = 1/((1+prstp)**(tstep*(t.val-1)));
optlrsav = (dk + .004)/(dk + .004*elasmu + prstp)*gama;

solve CO2 maximizing UTILITY using nlp ;
solve CO2 maximizing UTILITY using nlp ;

* Calculate social cost of carbon and other variables

scc(t) = -1000*eeq.m(t)/(.00001+cc.m(t));
atfrac(t) = ((mat.l(t)-588)/(ccatot.l(t)+.000001  ));
atfrac2010(t) = ((mat.l(t)-mat0)/(.00001+ccatot.l(t)-ccatot.l('1')  ));
ppm(t)    = mat.l(t)/2.13;
put /
put /"DICE-2106R3 disc = 1%";
put /"ifopt =" ifopt;
put "elastmu =" elasmu;
put "prstp =" prstp;
put // "Period";
Loop (T, put T.val);
put / "Year" ;
Loop (T, put (2010+(TSTEP*T.val) ));
put / "Industrial Emissions GTCO2 per year" ;
Loop (T, put EIND.l(T));
put / "Atmospheric concentration C (ppm)" ;
Loop (T, put (MAT.l(T)/2.13));
put / "Atmospheric Temperature " ;
Loop (T, put TATM.l(T));
put / "Output Net Net) " ;
Loop (T, put Y.l(T));
put / "Climate Damages fraction output" ;
Loop (T, put DAMFRAC.l(T));
put / "Consumption Per Capita " ;
Loop (T, put CPC.l(T));
put / "Carbon Price (per t CO2)" ;
Loop (T, put cprice.l(T));
put / "Emissions Control Rate" ;
Loop (T, put MIU.l(T));
put / "Social cost of carbon" ;
Loop (T, put scc(T));
put / "Interest Rate " ;
Loop (T, put RI.l(T));
put / "Population" ;
Loop (T, put L(T));
put / "TFP" ;
Loop (T, put AL(T));
put / "Output gross,gross" ;
Loop (T, put YGROSS.L(t));
put / "Change tfp" ;
Loop (T, put ga(t));
put / "Capital" ;
Loop (T, put k.l(t));
 put / "s" ;
Loop (T, put s.l(t));
  put / "I" ;
Loop (T, put I.l(t));
   put / "Y gross net" ;
Loop (T, put ynet.l(t));
   put / "damages" ;
Loop (T, put damages.l(t));
  put / "damfrac" ;
Loop (T, put damfrac.l(t));
 put / "abatement" ;
Loop (T, put abatecost.l(t));
 put / "sigma" ;
Loop (T, put sigma(t));
 put / "Forcings" ;
Loop (T, put forc.l(t));
put / "Other Forcings" ;
Loop (T, put forcoth(t));
put / "Period utilty" ;
Loop (T, put periodu.l(t));
put / "Consumption" ;
Loop (T, put C.l(t));
put / "Objective" ;
put utility.l;
put / "Land emissions" ;
Loop (T, put etree(t));
put / "Cumulative ind emissions" ;
Loop (T, put cca.l(t));
put / "Cumulative total emissions" ;
Loop (T, put ccatot.l(t));
put / "Atmospheric concentrations Gt" ;
Loop (T, put mat.l(t));
put / "Atmospheric concentrations ppm" ;
Loop (T, put ppm(t));
put / "Total Emissions GTCO2 per year" ;
Loop (T, put E.l(T));
put / "Atmospheric concentrations upper" ;
Loop (T, put mu.l(t));
put / "Atmospheric concentrations lower" ;
Loop (T, put ml.l(t));
put / "Atmospheric fraction since 1850" ;
Loop (T, put atfrac(t));
put / "Atmospheric fraction since 2010" ;
Loop (T, put atfrac2010(t));
put / "Average forward looking 100 year temp" ;
Loop (T, put avtemp100.l(t));
put / "Average forward looking 200 year temp" ;
Loop (T, put avtemp200.l(t));