%[text] domn = D2T
parm = [1, 1; 4, 0; 2, 3]';
%%
fcn = Fcn("D2", "x+1");
fun = fcn.getFun;
intVal = fcn.int("D2T").getFun; intVal([0; 0], parm)
GInt("D2T", 1).eval(fun, parm)
%%
fcn = Fcn("D2", "x^2+x+1");
fun = fcn.getFun;
intVal = fcn.int("D2T").getFun; intVal([0; 0], parm)
GInt("D2T", 2).eval(fun, parm)
%%
fcn = Fcn("D2", "x^3+x^2+x+1");
fun = fcn.getFun;
intVal = fcn.int("D2T").getFun; intVal([0; 0], parm)
GInt("D2T", 3).eval(fun, parm)
%%
fcn = Fcn("D2", "x^4+x^3+x^2+x+1");
fun = fcn.getFun;
intVal = fcn.int("D2T").getFun; intVal([0; 0], parm)
GInt("D2T", 4).eval(fun, parm)
%%
%[text] domn = D2L
parm = [1, 1; 2, 3]';
%%
fcn = Fcn("D2", "x+1");
fun = fcn.getFun;
intVal = fcn.int("D2L").getFun; intVal(0, parm)
GInt("D2L", 1).eval(fun, parm)
%%
fcn = Fcn("D2", "x^2+x+1");
fun = fcn.getFun;
intVal = fcn.int("D2L").getFun; intVal(0, parm)
GInt("D2L", 2).eval(fun, parm)
%%
fcn = Fcn("D2", "x^3+x^2+x+1");
fun = fcn.getFun;
intVal = fcn.int("D2L").getFun; intVal(0, parm)
GInt("D2L", 3).eval(fun, parm)
%%
fcn = Fcn("D2", "x^4+x^3+x^2+x+1");
fun = fcn.getFun;
intVal = fcn.int("D2L").getFun; intVal(0, parm)
GInt("D2L", 4).eval(fun, parm)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
