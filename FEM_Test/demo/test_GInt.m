%[text] # GInt - Gauss integration
%%
%[text] ## Integration on D2T
parm = [-1, 0; 1, 0; 0, 1]';
%%
%[text] ### Degree 1
fcn = Fcn("D2", "x+1");
disp(GInt("D2T", 1).eval(fcn.getFun, parm));

intVal = fcn.int("D2T").getFun;
disp(intVal([0; 0], parm));
%%
%[text] ### Degree 2
fcn = Fcn("D2", "x^2+x+1");
disp(GInt("D2T", 2).eval(fcn.getFun, parm));

intVal = fcn.int("D2T").getFun;
disp(intVal([0; 0], parm));
%%
%[text] ### Degree 3
fcn = Fcn("D2", "x^3+x^2+x+1");
disp(GInt("D2T", 3).eval(fcn.getFun, parm));

intVal = fcn.int("D2T").getFun;
disp(intVal([0; 0], parm));
%%
%[text] ### Degree 4
fcn = Fcn("D2", "x^4+x^3+x^2+x+1");
disp(GInt("D2T", 4).eval(fcn.getFun, parm));

intVal = fcn.int("D2T").getFun;
disp(intVal([0; 0], parm));
%%
%[text] ## Integration on D2TR
%[text] ### Degree 1
fcn = Fcn("D2R", "l+1");
disp(GInt("D2T", 1).eval(fcn.getFun));

intVal = fcn.int("D2TR").getFun;
disp(intVal([0; 0]));
%%
%[text] ### Degree 2
fcn = Fcn("D2R", "l^2+l+1");
disp(GInt("D2T", 2).eval(fcn.getFun));

intVal = fcn.int("D2TR").getFun;
disp(intVal([0; 0]));
%%
%[text] ### Degree 3
fcn = Fcn("D2R", "l^3+l^2+l+1");
disp(GInt("D2T", 3).eval(fcn.getFun));

intVal = fcn.int("D2TR").getFun;
disp(intVal([0; 0]));
%%
%[text] ### Degree 4
fcn = Fcn("D2R", "l^4+l^3+l^2+l+1");
disp(GInt("D2T", 4).eval(fcn.getFun));

intVal = fcn.int("D2TR").getFun;
disp(intVal([0; 0]));
%%
%[text] ## Integration on D2L
parm = [0, 0; 1, 1]';
%[text] ### Degree 1
fcn = Fcn("D2", "x+1");
disp(GInt("D2L", 1).eval(fcn.getFun, parm));

intVal = fcn.int("D2L").getFun;
disp(intVal(0, parm));
%%
%[text] ### Degree 3
fcn = Fcn("D2", "x^3+x^2+x+1");
disp(GInt("D2L", 3).eval(fcn.getFun, parm));

intVal = fcn.int("D2L").getFun;
disp(intVal(0, parm));
%%
%[text] ### Degree 5
fcn = Fcn("D2", "x^5+x^4+x^3+x^2+x+1");
disp(GInt("D2L", 5).eval(fcn.getFun, parm));

intVal = fcn.int("D2L").getFun;
disp(intVal(0, parm));
%%
%[text] ### Degree 7
fcn = Fcn("D2", "x^7+x^6+x^5+x^4+x^3+x^2+x+1");
disp(GInt("D2L", 7).eval(fcn.getFun, parm));

intVal = fcn.int("D2L").getFun;
disp(intVal(0, parm));
%%
%[text] ## Integration on D2LR
%[text] ### Degree 1
fcn = Fcn("D2R1", "s+1");
disp(GInt("D2L", 1).eval(fcn.getFun));

intVal = fcn.int("D2LR").getFun;
disp(intVal(0));
%%
%[text] ### Degree 3
fcn = Fcn("D2R1", "s^3+s^2+s+1");
disp(GInt("D2L", 3).eval(fcn.getFun));

intVal = fcn.int("D2LR").getFun;
disp(intVal(0));
%%
%[text] ### Degree 5
fcn = Fcn("D2R1", "s^5+s^4+s^3+s^2+s+1");
disp(GInt("D2L", 5).eval(fcn.getFun));

intVal = fcn.int("D2LR").getFun;
disp(intVal(0));
%%
%[text] ### Degree 7
fcn = Fcn("D2R1", "s^7+s^6+s^5+s^4+s^3+s^2+s+1");
disp(GInt("D2L", 7).eval(fcn.getFun));

intVal = fcn.int("D2LR").getFun;
disp(intVal(0));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
