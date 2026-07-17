%[text] Lagrange element P1
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
base = P1_FE.base;

base(1).fun
base(1).eval("[x1;y1]")
base(1).eval("[x2;y2]")
base(1).eval("[x3;y3]")
%%
%[text] Lagrange element P2
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0])]);
base = P2_FE.base;


base(1).fun
base(1).eval("[x1;y1]")
base(1).eval("[x2;y2]")
base(1).eval("[(x1+x2)/2;(y1+y2)/2]")

base(4).fun
base(4).eval("[x1;y1]")
base(4).eval("[x2;y2]")
base(4).eval("[(x1+x2)/2;(y1+y2)/2]")
%%
%[text] CR element
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0]));
base = CR_FE.base;

base(1).fun
base(1).eval("[x1/2+x2/2;y1/2+y2/2]")
base(1).eval("[x2/2+x3/2;y2/2+y3/2]")
base(1).eval("[x3/2+x1/2;y3/2+y1/2]")
%%
%[text] RT0 element
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'",...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true));
base = RT0_FE.base;

base(1).fun
int(dot(base(1), MshEnt("D2T").UNV(1)), "D2L", 1).fun
int(dot(base(1), MshEnt("D2T").UNV(2)), "D2L", 2).fun
int(dot(base(1), MshEnt("D2T").UNV(3)), "D2L", 3).fun
%%
%[text] DG P1 element
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
base = DG1_FE.base;

base(1).fun
base(1).eval("[x1;y1]")
base(1).eval("[x2;y2]")
base(1).eval("[x3;y3]")
%%
%[text] SDG P1 element
%[text] - scalar \
SDG_1S_FE = FE("D2T", "[1,x,y]", ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, [Fcn("D2R1", "1"), Fcn("D2R1", "s")], [0; 0], "EntIdx", 1), ...
    MoDoF("D2", MshEnt("D2T").msh, 2, Fcn.cst(1), [0; 0])]);
base = SDG_1S_FE.base;

base(1).fun
int(base(1).clrParm.tfm("D2LR") .* Tfm("D2L").JNorm, "D2LR").fun
int(base(1).clrParm.tfm("D2LR") .* Fcn("D2R1", "s") .* Tfm("D2L").JNorm, "D2LR").fun
int(base(1), "D2T").fun

base(2).fun
int(base(2).clrParm.tfm("D2LR") .* Tfm("D2L").JNorm, "D2LR").fun
int(base(2).clrParm.tfm("D2LR") .* Fcn("D2R1", "s") .* Tfm("D2L").JNorm, "D2LR").fun
int(base(2), "D2T").fun

base(3).fun
int(base(3).clrParm.tfm("D2LR") .* Tfm("D2L").JNorm, "D2LR").fun
int(base(3).clrParm.tfm("D2LR") .* Fcn("D2R1", "s") .* Tfm("D2L").JNorm, "D2LR").fun
int(base(3), "D2T").fun
%%
SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0; 0], "EntIdx", 1), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "EntIdx", 3, "share", false)]);
base = SDG_1S_FE.base;

base(1).fun
base(1).eval("[x1;y1]")
base(1).eval("[x2;y2]")
base(1).eval("[x3;y3]")
%%
%[text] - vector \
SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, [Fcn("D2R1", "1"), Fcn("D2R1", "s")], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, [Fcn("D2R1", "1"), Fcn("D2R1", "s")], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", 1, "orien", true, "share", false)]);
base = SDG_1V_FE.base;

base(5).fun
int(dot(base(5), MshEnt("D2T").UNV(1)).clrParm.tfm("D2LR") .* Tfm("D2L").JNorm, "D2LR").fun
int(dot(base(5), MshEnt("D2T").UNV(1)).clrParm.tfm("D2LR") .* Fcn("D2R1", "s") .* Tfm("D2L").JNorm, "D2LR").fun

base(6).fun
int(dot(base(6), MshEnt("D2T").UNV(1)).clrParm.tfm("D2LR") .* Tfm("D2L").JNorm, "D2LR").fun
int(dot(base(6), MshEnt("D2T").UNV(1)).clrParm.tfm("D2LR") .* Fcn("D2R1", "s") .* Tfm("D2L").JNorm, "D2LR").fun
%%
SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", 1, "orien", true, "share", false)]);
base = SDG_1V_FE.base;

base(1).fun
dot(base(1), MshEnt("D2T").UNV(2)).eval("[x2;y2]")
dot(base(1), MshEnt("D2T").UNV(2)).eval("[x3;y3]")
%%
%[text] - matrix \
SDG_1M_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 2]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", 1, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", [2, 3], "share", false, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UTV], "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);
base = SDG_1M_FE.base;

% load("SDG_1M_FE.mat");
% base = SDG_1M_FE.base;

base(1).fun
dot(base(1) * MshEnt("D2T").UNV(1), MshEnt("D2T").UNV(1)).eval("[x1;y1]")
dot(base(1) * MshEnt("D2T").UNV(1), MshEnt("D2T").UNV(1)).eval("[x2;y2]")
dot(base(1) * MshEnt("D2T").UNV(1), MshEnt("D2T").UTV(1)).eval("[x1;y1]")
dot(base(1) * MshEnt("D2T").UNV(1), MshEnt("D2T").UTV(1)).eval("[x2;y2]")

base(7).fun
dot(base(7) * MshEnt("D2T").UNV(1), MshEnt("D2T").UNV(1)).eval("[x1;y1]")
dot(base(7) * MshEnt("D2T").UNV(1), MshEnt("D2T").UNV(1)).eval("[x2;y2]")
dot(base(7) * MshEnt("D2T").UNV(1), MshEnt("D2T").UTV(1)).eval("[x1;y1]")
dot(base(7) * MshEnt("D2T").UNV(1), MshEnt("D2T").UTV(1)).eval("[x2;y2]")
%%
%[text] SDG P0 element
%[text] - scalar \
SDG_0S_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0], "EntIdx", 1));
base = SDG_0S_FE.base;

base.fun
%%
%[text] - vector \
SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true));
base = SDG_0V_FE.base;

base(1).fun
dot(base(1), MshEnt("D2T").UNV(2)).simplify.fun
%%
%[text] - matrix \
SDG_0M_FE = FE("D2T", FE.repFS("1", [2, 2]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", 1, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UTV], "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);
base = SDG_0M_FE.base;

base(1).fun
dot(base(1) * MshEnt("D2T").UNV(1), MshEnt("D2T").UNV(1)).simplify.fun
dot(base(1) * MshEnt("D2T").UNV(1), MshEnt("D2T").UTV(1)).simplify.fun

base(2).fun
dot(base(2) * MshEnt("D2T").UNV(1), MshEnt("D2T").UNV(1)).simplify.fun
dot(base(2) * MshEnt("D2T").UNV(1), MshEnt("D2T").UTV(1)).simplify.fun
%%
%[text] Hybrid P1 element
P1H_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
P1H_FE.base.fun
%%
%[text] FE.repFS
FE.repFS("[1,x,y]", [3, 1])
FE.repFS("[1,x,y]", [2, 3])

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
