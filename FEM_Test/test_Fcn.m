%[text] # Fcn - Symbolic Function
%%
%[text] ## Constructor
disp(Fcn("D2", "x+y").var);

disp(Fcn("D2R", "l+m;").var);

disp(Fcn("D2R1", "s").var);

disp(Fcn("D2T", "x+y").var);
disp(Fcn("D2T", "x+y").parm);

disp(Fcn("D2TR", "l+m").var);
disp(Fcn("D2TR", "l+m").parm);

disp(Fcn("D2L", "x+y").var);
disp(Fcn("D2L", "x+y").parm);

disp(Fcn("D2LR", "s").var);
disp(Fcn("D2LR", "s").parm);

disp(Fcn("D2", "x+y").fun);
disp(Fcn("D2", "[x;y]").fun);
disp(Fcn("D2T", "[x,y;y,x]").fun);

disp(Fcn("D2T", "[c1*x;c2*y]", "[c1;c2]").fun);
disp(Fcn("D2T", "[c1*x;c2*y]", "[c1;c2]").coef);
%%
%[text] ## Get functions
fcn = Fcn("D2T", "[c1*x;c2*y]", "[c1;c2]");
disp(fcn.nVar);
disp(fcn.nFun);
disp(fcn.sParm);
disp(fcn.nCoef);

fcn = Fcn("D2T", "[x,y;y,x]");
disp(fcn.nFun);
%%
%[text] ### getDomn
fcns = [Fcn("D2", "1"), Fcn("D2T", "1")];
disp(fcns.getDomn);

% % Invalid: incompatible domains raise an error.
% fcns = [Fcn("D2T", "1"), Fcn("D2TR", "1")];
% disp(fcns.getDomn);
%%
%[text] ### getFun
disp(Fcn("D2", "x+y").getFun);
disp(Fcn("D2T", "x+y+x1+x2").getFun);
disp(Fcn("D2T", "c1*x+c2*y+x1+x2", "[c1;c2]").getFun);
disp(Fcn("D2", "c1*x+c2*y+x1+x2").getFun("parm", {"[x1,x2;y1,y2]"}, "coef", {"[c1;c2]"}));
%%
fun = Fcn("D2T", "c1*x+c2*y+x1+x2", "[c1;c2]").getFun;
disp(fun(str2sym("[P;P]"), str2sym("[X1,X2,X3;Y1,Y2,Y3]"), str2sym("[C1;C2]")));
disp(fun(str2sym("[P1,P2,P3;P1,P2,P3]"), str2sym("[X1,X2,X3;Y1,Y2,Y3]"), str2sym("[C1;C2]")));

fun = Fcn("D2", "1").getFun;
disp(fun([1, 2, 3; 1, 2, 3]));
%%
%[text] ## `Overload intrinsic operators`
fcn1 = Fcn("D2", "x+y");
fcn2 = Fcn("D2", "[sin(x);cos(y)]");

fcn0 = fcn1 + fcn2; disp(fcn0.fun);
fcn0 = fcn1 - fcn2; disp(fcn0.fun);
fcn0 = -fcn1; disp(fcn0.fun);
fcn0 = fcn1 * fcn2; disp(fcn0.fun);
fcn0 = fcn1 * 2; disp(fcn0.fun);
fcn0 = fcn1 \ fcn2; disp(fcn0.fun);
fcn0 = fcn1 \ 2; disp(fcn0.fun);

fcn1 = Fcn("D2", "[x+y;x-y]");
fcn2 = Fcn("D2", "[sin(x);cos(y)]");

fcn0 = fcn1 .* fcn2; disp(fcn0.fun);
fcn0 = fcn1 .\ fcn2; disp(fcn0.fun);
fcn0 = fcn1 .^ 2; disp(fcn0.fun);
fcn0 = fcn2.'; disp(fcn0.fun);

% % Invalid: incompatible domains and coefficients raise an error.
% Fcn("D2", "1") .* Fcn("D2R", "1")
% Fcn("D2", "1", "c1") .* Fcn("D2", "1", "c2")
%%
%[text] ### Domn merge
fcn = Fcn("VOID", "1") .* Fcn("D2", "1"); disp(fcn.domn);
fcn = Fcn("D2", "1") .* Fcn("D2T", "1"); disp(fcn.domn);
fcn = Fcn("D2", "1") .* Fcn("D2L", "1"); disp(fcn.domn);
fcn = Fcn("D2R", "1") .* Fcn("D2TR", "1"); disp(fcn.domn);
fcn = Fcn("D2R1", "1") .* Fcn("D2LR", "1"); disp(fcn.domn);
fcn = Fcn("D2L", "1") .* Fcn("D2T", "1"); disp(fcn.domn);
%%
%[text] ### `Coefficient concatenation`
fcn1 = Fcn("D2", "1", "[c1;c2]");
fcn2 = Fcn("D2", "1", "[c3;c4]");

fcn0 = fcn1 + fcn2; disp(fcn0.coef);
fcn0 = fcn1 - fcn2; disp(fcn0.coef);
%%
%[text] ## Overloaded intrinsic functions
fcn1 = Fcn("D2", "[x+y;x-y]");
fcn2 = Fcn("D2", "[sin(x);cos(y)]");
fcn3 = Fcn("D2", "(1+x)^2/(1+2*x+x^2)");

disp(sum(fcn1).fun);
disp(abs(fcn1).fun);
disp(dot(fcn1, fcn2).fun);
disp(diff(fcn1, "y", 1).fun);
disp(simplify(fcn3).fun);
disp(subs(fcn1, "x", "1").fun);
%%
%[text] ## `Mathematical operations`
%%
%[text] ### eval
disp(Fcn("D2T", "x+y").eval([1; 2]));
disp(Fcn("D2T", "x+y").eval("[x1;y1]"));
%%
%[text] ### comb
disp(comb([Fcn("D2", "x"), Fcn("D2", "y")], [1; 2]).fun);

disp(comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], [1; 2]).fun);
disp(comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], [1; 2]).coef);

disp(comb([Fcn("D2", "x"), Fcn("D2", "y")], "[c1;c2]").fun);
disp(comb([Fcn("D2", "x"), Fcn("D2", "y")], "[c1;c2]").coef);

% % Invalid: functions already carry coefficients, and weights are of symbolic type.
% disp(comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], "[c1;c2]").fun);
%%
%[text] ### tfm
%[text] #### D2TR $\\rightarrow$ D2T
fcn = Fcn("D2TR", "[l;m]").tfm("D2T");
disp(fcn.domn);
disp(fcn.fun);
disp(fcn.eval("[x1;y1]"));
disp(fcn.eval("[x2;y2]"));
disp(fcn.eval("[x3;y3]"));
%%
%[text] #### D2T $\\rightarrow$ D2TR
fcn = Fcn("D2T", "[x;y]").tfm("D2TR");
disp(fcn.domn);
disp(fcn.fun);
disp(fcn.eval("[0;0]"));
disp(fcn.eval("[1;0]"));
disp(fcn.eval("[0;1]"));
%%
%[text] #### D2L $\\rightarrow$ D2LR
fcn = Fcn("D2L", "[x;y]").tfm("D2LR");
disp(fcn.domn);
disp(fcn.fun);
disp(fcn.eval("0"));
disp(fcn.eval("0.5"));
disp(fcn.eval("1"));
%%
%[text] #### array operation
fcns = [Fcn("D2T", "[x;y]"), Fcn("D2T", "[y;x]")];
disp(fcns.tfm("D2TR").fun);
%%
%[text] ### dif
%[text] #### vector function
fcn = Fcn("D2T", ["sin(x)+cos(y)"; "x-y"]);
%[text] $\[\\partial\_x v\_x \\, ; v\_y\]$
disp(fcn.dif([1, 0; 0, 0]).fun);
%[text] $\[\\partial\_y v\_x \\, ; v\_y\]$
disp(fcn.dif([0, 0; 1, 0]).fun);
%[text] $\[v\_x \\, ; \\partial\_x v\_y\]$
disp(fcn.dif([0, 1; 0, 0]).fun);
%[text] $\[v\_x \\, ; \\partial\_y v\_y\]$
disp(fcn.dif([0, 0; 0, 1]).fun);
%[text] $\[v\_x \\, ; 0\]$
disp(fcn.dif([0, nan; 0, nan]).fun);
%[text] $\[0 \\, ; v\_y\]$
disp(fcn.dif([nan, 0; nan, 0]).fun);
%%
%[text] #### matrix function
fcn = Fcn("D2T", "[x,1;1,y]");
%[text] #### 
%[text] $\[ \\partial\_x v\_{xx} \\, , v\_{xy} \] \\\\\n\[ \\partial\_x v\_{yx} \\, , v\_{yy} \]$
disp(fcn.dif([1, 0; 1, 0; 0, 0; 0, 0]').fun);
%[text] #### 
%[text] $\[ v\_{xx} \\, , \\partial\_y v\_{xy} \] \\\\\n\[ v\_{yx} \\, , \\partial\_y v\_{yy} \]$
disp(fcn.dif([0, 0; 0, 0; 0, 1; 0, 1]').fun);
%%
%[text] #### array operation
fcns = [Fcn("D2T", "sin(x)+cos(y)"), Fcn("D2T", "x-y")];
disp([fcns.dif([1; 0]).fun]);
%%
%[text] #### dimension upgrade
%[text] scalar $\\rightarrow$ vector
fcn = Fcn("D2T", "x-y");
grad = cat(3, [1; 0], [0; 1]);
disp(fcn.dif(grad).fun);
%[text] vector $\\rightarrow$ matrix
fcn = Fcn("D2T", "[x-y;x^2+y^2]");
grad = cat(3, [1, 1; 0, 0], [0, 0; 1, 1]);
disp(fcn.dif(grad).fun);
%%
%[text] ### int
%[text] $\\int\_K v \\, dx \\, dy$
disp(Fcn("D2", "1").int("D2T").fun);
disp(Fcn("D2T", "x").int("D2T").fun);
%%
%[text] $\\int\_{\\hat{K}} \\hat{v} \\, d \\hat{x} \\, d \\hat{y}$
disp(Fcn("D2R", "1").int("D2TR").fun);
disp(Fcn("D2TR", "l").int("D2TR").fun);
%%
%[text] $\\int\_e v \\, ds$
disp(Fcn("D2", "1").int("D2L").fun);
disp(Fcn("D2L", "x").int("D2L").fun);
%%
%[text] $\\int\_e v \\, ds$
disp(Fcn("D2T", "1").int("D2L", 2).fun);
disp(Fcn("D2T", "x").int("D2L", 2).fun);
%%
%[text] $\\int\_{\\hat{e}} \\hat{v} \\, d \\hat{s}$
disp(Fcn("D2R1", "1").int("D2LR").fun);
disp(Fcn("D2LR", "s").int("D2LR").fun);
%%
%[text] $\\int\_{\\hat{e}} \\hat{v} \\, d \\hat{s}$
disp(Fcn("D2R", "1").int("D2LR", 2).fun);
disp(Fcn("D2TR", "s").int("D2LR", 2).fun);
%%
%[text] ## `Parameter and coefficient management`
fcn = Fcn("D2T", "c1*x+c2*y+x1+x2", "[c1;c2]");
%[text] ### clrParm
fcn0 = fcn.clrParm;
disp(fcn0.domn); disp(fcn0.fun); disp(fcn0.parm); disp(fcn0.coef);
%[text] ### clrCoef
fcn0 = fcn.clrCoef;
disp(fcn0.domn); disp(fcn0.fun); disp(fcn0.parm); disp(fcn0.coef);
%[text] ### subParm
fcn0 = fcn.subParm("[X1,X2,X3;Y1,Y2,Y3]");
disp(fcn0.domn); disp(fcn0.fun); disp(fcn0.parm); disp(fcn0.coef);
%[text] ### subCoef
fcn0 = fcn.subCoef("[C1;C2]");
disp(fcn0.domn); disp(fcn0.fun); disp(fcn0.parm); disp(fcn0.coef);
%%
%[text] ## Static functions
%[text] ### `cst`
disp(Fcn.cst(1).domn);
disp(Fcn.cst(1).fun);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
