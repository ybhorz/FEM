%[text] ## Fcn
%[text] ### Constructor
Fcn("D2", "[x+y;x-y]").var
Fcn("D2", "[x+y;x-y]").parm
Fcn("D2", "[x+y;x-y]").fun
%%
Fcn("D2R", "[l+m;l-m]").var
Fcn("D2R", "[l+m;l-m]").parm
Fcn("D2R", "[l+m;l-m]").fun
%%
Fcn("D2R1", "[s;s]").var
Fcn("D2R1", "[s;s]").parm
Fcn("D2R1", "[s;s]").fun
%%
Fcn("D2T", "[x+y;x-y]").var
Fcn("D2T", "[x+y;x-y]").parm
Fcn("D2T", "[x+y;x-y]").fun
%%
Fcn("D2TR", "[l+m;l-m]").var
Fcn("D2TR", "[l+m;l-m]").parm
Fcn("D2TR", "[l+m;l-m]").fun
%%
Fcn("D2L", "[x+y;x-y]").var
Fcn("D2L", "[x+y;x-y]").parm
Fcn("D2L", "[x+y;x-y]").fun
%%
Fcn("D2LR", "[s;s]").var
Fcn("D2LR", "[s;s]").parm
Fcn("D2LR", "[s;s]").fun
%%
Fcn("D2T", "[x+y;x-y]", "[c1;c2;c3]").coef
%%
Fcn("D2T", "[x,1;1,y]").fun
%%
%[text] ### Get functions
fcn = Fcn("D2T", "[x+y;x-y]", "[c1;c2;c3]");
fcn.nVar
fcn.nFun
fcn.sParm
fcn.nCoef
%%
Fcn("D2T", "[x,1;1,y]").nFun
%%
%[text] getDomn
fcns = [Fcn("D2", "1"), Fcn("D2T", "1")];
fcns.getDomn

% Check.
% fcns = [Fcn("D2T", "1"), Fcn("D2TR", "1")];
% fcns.getDomn
%%
%[text] getFun
Fcn("D2", "x+y").getFun
Fcn("D2T", "x+y+x1+x2").getFun
Fcn("D2T", "c1*x+c2*y+x1+x2", "[c1;c2]").getFun
Fcn("D2", "x+y+x1+x2").getFun("coef", {"[x1,x2;y1,y2]"})
%%
fun = Fcn("D2T", "x+x2").getFun;
fun = fun([1, 2, 3; 4, 5, 6], [1, 2, 3; 4, 5, 6])

fun = Fcn("D2T", "1").getFun;
fun = fun([1, 2, 3; 4, 5, 6], [1, 2, 3; 4, 5, 6])

fun = Fcn("D2T", "x2").getFun;
fun = fun(str2sym("[1,2,3;4,5,6]"), str2sym("[x1,x2,x3;y1,y2,y3]"))
%%
%[text] ### Overload intrinsic operators
fcn1 = Fcn("D2T", "x+y");
fcn2 = Fcn("D2T", "[sin(x);cos(y)]");
fcn0 = fcn1 + fcn2; fcn0.fun
fcn0 = fcn1 - fcn2; fcn0.fun
fcn0 = -fcn2; fcn0.fun
fcn0 = fcn1 * fcn2; fcn0.fun
fcn0 = fcn1 * 2; fcn0.fun
fcn0 = fcn1 \ fcn2; fcn0.fun
fcn0 = fcn1 \ 2; fcn0.fun
fcn0 = fcn1 .* fcn2; fcn0.fun
fcn0 = fcn1 .\ fcn2; fcn0.fun
fcn0 = fcn1 .^ 2; fcn0.fun
fcn0 = fcn2.'; fcn0.fun

% Check
% Fcn("D2") .* Fcn("D2R")
% Fcn("D2", "x", "c1") .* Fcn("D2", "x", "c2")
%%
%[text] merge
fcn = Fcn("VOID", "1") .* Fcn("D2", "1"); fcn.domn
fcn = Fcn("D2", "1") .* Fcn("VOID", "1"); fcn.domn
fcn = Fcn("D2", "1") .* Fcn("D2T", "1"); fcn.domn
fcn = Fcn("D2T", "1") .* Fcn("D2", "1"); fcn.domn
%%
%[text] plus
fcn1 = Fcn("D2T", "x+y", "[c1;c2]");
fcn2 = Fcn("D2T", "[sin(x);cos(y)]", "[b1;b2]");
fcn0 = fcn1 + fcn2; fcn0.coef
%%
%[text] ### `Overload intrinsic functions`
fcn1 = Fcn("D2T", "[x+y;x-y]");
fcn2 = Fcn("D2T", "[sin(x);cos(y)]");

sum(fcn1).fun
abs(fcn1).fun
dot(fcn1, fcn2).fun
diff(fcn1, "x", 1).fun
%%
%[text] ### Mathematical operations
%[text] eval
Fcn("D2T", "x+y").eval([1; 2])
Fcn("D2T", "x+y").eval("[x1;y1]")
%%
%[text] comb
comb([Fcn("D2", "x"), Fcn("D2", "y")], [1;2]).fun
comb([Fcn("D2", "x"), Fcn("D2", "y")], [1;2]).coef

comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], [1;2]).fun
comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], [1;2]).coef


comb([Fcn("D2", "x"), Fcn("D2", "y")], "[c1;c2]").fun
comb([Fcn("D2", "x"), Fcn("D2", "y")], "[c1;c2]").coef

% Check
% comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], "[c1;c2]").fun
% comb([Fcn("D2", "x", "[c1;c2]"), Fcn("D2", "y", "[b1;b2]")], "[c1;c2]").coef
%%
%[text] tfm
fcn = Fcn("D2TR", "[l;m]").tfm("D2T");
fcn.domn
fcn.var
fcn.fun
fcn.parm
fcn.eval("[x1;y1]")
fcn.eval("[x2;y2]")
fcn.eval("[x3;y3]")
%%
fcn = Fcn("D2T", "[x;y]").tfm("D2TR");
fcn.domn
fcn.var
fcn.fun
fcn.parm
fcn.eval("[0;0]")
fcn.eval("[1;0]")
fcn.eval("[0;1]")
%%
fcn = Fcn("D2L", "[x;y]").tfm("D2LR");
fcn.domn
fcn.var
fcn.fun
fcn.parm
fcn.eval("0")
fcn.eval("0.5")
fcn.eval("1")
%%
fcns = [Fcn("D2T", "[x;y]"), Fcn("D2T", "[y;x]")];
fcns.tfm("D2TR").fun
%%
%[text] dif
fcn = Fcn("D2T", ["sin(x)+cos(y)"; "x^3+y^4"]);
fcn.dif([1, 0; 0, 0]).fun
fcn.dif([0, 0; 1, 0]).fun
fcn.dif([0, 1; 0, 0]).fun
fcn.dif([0, 0; 0, 1]).fun
fcn.dif([0, nan; 0, nan]).fun
%%
fcn = Fcn("D2T", "[x,1;1,y]");
fcn.dif([1,0;1,0;0,0;0,0]').fun
fcn.dif([0,0;0,0;0,1;0,1]').fun
%%
fcns = [Fcn("D2T", "sin(x)+cos(y)"), Fcn("D2T", "x^3+y^4")];
fcns.dif([1; 0]).fun
%%
fcn = Fcn("D2T", "x*y");
grad = cat(3, [1;0], [0;1]);
fcn.dif(grad).fun
%%
fcn = Fcn("D2T", "[1+x+y;1+x^2+y^2]");
grad = cat(3, [1,1;0,0], [0,0;1,1]);
fcn.dif(grad).fun
%%
%[text] int
Fcn("D2", "1").int("D2T")
Fcn("D2", "x").int("D2T")
Fcn("D2T", "x").int("D2T")
%%
Fcn("D2R", "1").int("D2TR")
Fcn("D2R", "l").int("D2TR")
Fcn("D2TR", "l").int("D2TR")
%%
Fcn("D2", "1").int("D2L")
Fcn("D2", "x").int("D2L")
Fcn("D2L", "x").int("D2L")
%%
Fcn("D2T", "1").int("D2L", 2)
Fcn("D2T", "x").int("D2L", 2)
%%
Fcn("D2R1", "1").int("D2LR")
Fcn("D2R1", "s").int("D2LR")
Fcn("D2LR", "s").int("D2LR")
%%
Fcn("D2R", "1").int("D2LR", 2)
Fcn("D2R", "s").int("D2LR", 2)
Fcn("D2TR", "s").int("D2LR", 2)
%%
%[text] ### Parameter and coefficient management
%[text] subParm
fcn = Fcn("D2T", "x+y+x1+x2").subParm("[a,b,c;d,e,f]");
fcn.domn
fcn.fun
fcn.parm
%%
fcns = [Fcn("D2T", "x+y+x1+x2"), Fcn("D2T", "x+y+y1+y2")];
fcns.subParm("[a,b,c;d,e,f]").fun
%%
%[text] subCoef
fcn = Fcn("D2T", "c1*x+c2*y", "[c1;c2]").subCoef("[a;b]");
fcn.domn
fcn.fun
fcn.coef
%%
%[text] ### Static functions
%[text] cst
Fcn.cst(1)
%%
%[text] ## Tfm
%[text] Constructor
Tfm("D2T")
Tfm("D2L")
%%
%[text] JDet
Tfm("D2T").JDet.domn
Tfm("D2T").JDet.var
Tfm("D2T").JDet.fun
Tfm("D2T").JDet.parm
%%
%[text] JNorm
Tfm("D2L").JNorm.domn
Tfm("D2L").JNorm.var
Tfm("D2L").JNorm.fun
Tfm("D2L").JNorm.parm
%%
%[text] `refTfm`
Tfm("D2T").refTfm.domn
Tfm("D2T").refTfm.var
Tfm("D2T").refTfm.fun
Tfm("D2T").refTfm.parm

Tfm("D2T").refTfm.eval("[x1;y1]")
Tfm("D2T").refTfm.eval("[x2;y2]")
Tfm("D2T").refTfm.eval("[x3;y3]")
%%
%[text] orgTfm
Tfm("D2T").orgTfm.domn
Tfm("D2T").orgTfm.var
Tfm("D2T").orgTfm.fun
Tfm("D2T").orgTfm.parm

Tfm("D2T").orgTfm.eval("[0;0]")
Tfm("D2T").orgTfm.eval("[1;0]")
Tfm("D2T").orgTfm.eval("[0;1]")
%%
Tfm("D2L").orgTfm.domn
Tfm("D2L").orgTfm.var
Tfm("D2L").orgTfm.fun
Tfm("D2L").orgTfm.parm

Tfm("D2L").orgTfm.eval("0")
Tfm("D2L").orgTfm.eval("0.5")
Tfm("D2L").orgTfm.eval("1")

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
