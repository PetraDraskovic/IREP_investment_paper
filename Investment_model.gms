$macro c(name) equation name;\
name

set t /t1*t168/;
set s /s1*s4/;
set i /i1*i4/;
set j /j1*j3/;
set w /w1*w8/;

*CROPEX day ahead prices (import price profiles (first sheet) form Input_data_investment.xlsx)
table cropexDA(t,w)
*$call =xls2gms r=input!a1:i169 i=C:\Users\milja\OneDrive\Desktop\cijene_rep_tjedni.xlsx o=C:\Users\milja\OneDrive\Desktop\inc_file\prices.inc
$include C:\Users\milja\OneDrive\Desktop\inc_file\prices.inc
;


*wind curves (import wind profiles (second sheet) form Input_data_investment.xlsx)
table wind_real(t,w,s)
*$call =xls2gms r=input!a1:ag169 i=C:\Users\milja\OneDrive\Desktop\solari_rep_tjedni.xlsx o=C:\Users\milja\OneDrive\Desktop\inc_file\wr.inc
$include C:\Users\milja\OneDrive\Desktop\inc_file\wr.inc
;

scalars
eff /0.92/
M /10000/
M_dual /10000/
coeff /0.4/
P_ex /20/
lamdba_voda /0.397/
lamdba_h /5/
eff_voda /10/
Gamma /72/
soe0 /0/
Ibsse /400000/
Ibssp /400000/
Ipv /1000000/
Iel /1250000/
disc_rate /0.03/
Tbss /15/
Tpv/20/
Tel/20/
;



parameter prob(s)/
s1       0.25
s2       0.25
s3       0.25
s4       0.25
/;

parameter prob_w(w)/
w1       1
w2       11
w3       2
w4       8
w5       5
w6       6
w7       18
w8       1
/;

parameter r(i)/
i1       0
i2       0.23
i3       0.947
i4       1
/;

parameter f(i)/
i1       0.823
i2       0.658
i3       0.046
i4       0
/;


positive variable ch(t,w,s), dis(t,w,s), soe(t,w,s) ,
                  soc_p(t,w,s,j), windy_real(t,w,s), omega(w,s), z_dual(t,w,s), mu1(t,w,s), mu2(t,w,s),
                  el_crtano(t,w,s), el_h(t,w,s), x_h(t,w,s), p_max, soemax, C_bss, C_pv, P_pv, Pel, C_el;
binary variable x, bat(t,w,s), x1(t,w,s),x2(t,w,s), b_h(t,w,s);
free variable cost, deviation(t,w,s), market_pos(t,w),  y_dual(t,w,s),  el_pom1(t,w,s),  el_pom2(t,w,s),  el_pom3(t,w,s);

c(eq1).. cost =e= sum(w, prob_w(w)* (sum(s,prob(s)*(sum(t,cropexDA(t,w)*deviation(t,w,s)
          - z_dual(t,w,s) + cropexDA(t,w)* coeff*y_dual(t,w,s)
          + (x_h(t,w,s)*lamdba_h)- ((lamdba_voda * x_h(t,w,s) * eff_voda)/1000))
          - omega(w,s)*Gamma))
          + sum(t,cropexDA(t,w)*market_pos(t,w))))- C_bss - C_pv - C_el;


*annual costs constraints
c(eq50).. C_bss =e= (soemax*Ibsse + p_max*Ibssp)*((disc_rate*((disc_rate+1)**Tbss))/(((disc_rate+1)**Tbss)-1));
c(eq51).. C_pv =e=  P_pv*Ipv*((disc_rate*((disc_rate+1)**Tpv))/(((disc_rate+1)**Tpv)-1));
c(eq55).. C_el =e=  Pel*Iel*((disc_rate*((disc_rate+1)**Tel))/(((disc_rate+1)**Tel)-1));


*market position and deviation constraints
c(eq3(t,w,s)).. market_pos(t,w) + deviation(t,w,s) =e= dis(t,w,s) - ch(t,w,s) - el_h(t,w,s) + windy_real(t,w,s);
c(eq19(t,w,s)).. windy_real(t,w,s) =l= wind_real(t,w,s)*P_pv;



*grid exchange constraints
c(eq22(t,w)).. market_pos(t,w)  =l= P_ex ;
c(eq23(t,w)).. market_pos(t,w)  =g= -P_ex ;

c(eq22a(t,w,s)).. market_pos(t,w) + deviation(t,w,s) =l= P_ex ;
c(eq23a(t,w,s)).. market_pos(t,w) + deviation(t,w,s) =g= -P_ex ;

*battery constraints
c(eq6(t,w,s)).. soe(t,w,s) =l= soemax;
c(eq16(t,w,s)).. soe(t,w,s) =e= sum(j, soc_p(t,w,s,j));
c(eq17(t,w,s,j)).. soc_p(t,w,s,j) =l= sum(i$(ord(i)=ord(j)),R(i+1)-R(i))*(soemax) ;
c(eq18(t,w,s))$(ord(t) gt 1).. ch(t,w,s) *eff =l= F('i1')*soemax - sum(i$(ord(i) ge 2), sum(j$(ord(j)=ord(i)-1), soc_p(t-1,w,s,j))*(F(i-1)-F(i))/(R(i)-R(i-1)));
c(eq14(t,w,s))$(ord(t) gt 1).. soe(t,w,s) =e= soe(t-1,w,s) + eff*ch(t,w,s) -  dis(t,w,s)/eff;
c(eq15(t,w,s)).. soe('t1',w,s) =e= soe0 + eff*ch('t1',w,s) - dis('t1',w,s)/eff;
*c(eq7(t,w,s)).. -bat(t,w,s)*M =l= ch(t,w,s);
c(eq8(t,w,s)).. ch(t,w,s) =l= bat(t,w,s)*M  ;
c(eq9(t,w,s)).. ch(t,w,s) =l= p_max ;
c(eq11(t,w,s)).. dis(t,w,s) =l= (1-bat(t,w,s))*M  ;
c(eq12(t,w,s)).. dis(t,w,s) =l= p_max ;


*robust constraints
c(eq25(t,w,s)).. omega(w,s) + z_dual(t,w,s) =g= 2 * cropexDA(t,w) * coeff * y_dual(t,w,s);
c(eq26(t,w,s)).. y_dual(t,w,s) =g= deviation(t,w,s);
c(eq27(t,w,s)).. y_dual(t,w,s) =g= -deviation(t,w,s);
c(eq28(t,w,s))..  1-mu1(t,w,s)-mu2(t,w,s) =e= 0;
c(eq29(t,w,s))..  deviation(t,w,s)-y_dual(t,w,s) =g= -(1-x1(t,w,s))*M_dual;
c(eq30(t,w,s))..  mu1(t,w,s) =l= x1(t,w,s)*M_dual;
c(eq31(t,w,s))..  -deviation(t,w,s)-y_dual(t,w,s) =g= -(1-x2(t,w,s))*M_dual;
c(eq32(t,w,s))..  mu2(t,w,s) =l= x2(t,w,s)*M_dual;

*electrolyzer constraints
c(eq39(t,w,s))..  el_pom1(t,w,s) =l= el_h(t,w,s);
c(eq53(t,w,s))..  -M * b_h(t,w,s) =l= el_pom1(t,w,s);
c(eq60(t,w,s))..                      el_pom1(t,w,s) =l= M * b_h(t,w,s);
c(eq61(t,w,s))..  -M * (1 - b_h(t,w,s)) =l= el_pom1(t,w,s) - 0.1* Pel;
c(eq62(t,w,s))..                            el_pom1(t,w,s) - 0.1* Pel =l= M * (1 - b_h(t,w,s));

c(eq63(t,w,s))..  el_h(t,w,s) =l= el_pom2(t,w,s);
c(eq64(t,w,s))..  -M * b_h(t,w,s) =l= el_pom2(t,w,s);
c(eq65(t,w,s))..                      el_pom2(t,w,s) =l= M * b_h(t,w,s);
c(eq66(t,w,s))..  -M * (1 - b_h(t,w,s)) =l= el_pom2(t,w,s) - Pel;
c(eq67(t,w,s))..                            el_pom2(t,w,s) - Pel =l= M * (1 - b_h(t,w,s));



c(eq41(t,w,s)).. el_crtano(t,w,s) =e= (x_h(t,w,s)*39.4)/1000;
c(eq42(t,w,s)).. el_crtano(t,w,s) =e= 0.689*(el_h(t,w,s)) + el_pom3(t,w,s);

c(eq57(t,w,s)).. -M * b_h(t,w,s) =l= el_pom3(t,w,s);
c(eq68(t,w,s))..                     el_pom3(t,w,s) =l= M * b_h(t,w,s);
c(eq69(t,w,s)).. -M * (1 - b_h(t,w,s)) =l= el_pom3(t,w,s) - 0.011 * Pel;
c(eq70(t,w,s))..                           el_pom3(t,w,s) - 0.011 * Pel =l= M * (1 - b_h(t,w,s));




model project /all/;
option optCR=0.005;
option threads=2;
*option reslim=2000;
option mip=gurobi;
option rmip=gurobi;

file opt gurobi option file /gurobi.opt/;
put opt;
put '';
put 'method 2';
put 'IntFeasTol 1e-9';
putclose;

project.optfile=1;

solve project using MIP maximizing cost;
execute_unload "robust_stohast_bat_el.gdx"   ch, dis, el_h , soe,  x_h,  market_pos, soemax, p_max, deviation, Pel, P_pv, C_bss, C_pv, C_el, windy_real, cost, x_h ;

