clear;
node = Node([0, 0; 1, 0; 0, 1; 1, 1]', [-1, -2, -4, -3]);
elem = Elem([2, 3, 1; 3, 2, 4]', [5, -3, 1; -5, 4, -2]', [5, 7]);
edge = Edge([1, 2; 3, 4; 1, 3; 2, 4; 2, 3]', [1, 0; -2, 0; -1, 0; 2, 0; 1, -2]', [1, 3, 4, 2, 0]);
msh = Msh("D2T", node, elem, edge);
msh.figure;
%%
%[text] Lagrange element P1
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
P1_FES = FES(msh, P1_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.getFun;
FEfcn = P1_FES.proj(fcn); FEfun = FEfcn.getFun;
node = msh.node.coord(:, msh.elem.node(:, 2));
fun(node)
FEfun(node, FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] Lagrange element P2
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0])]);
P2_FES = FES(msh, P2_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.getFun;
FEfcn = P2_FES.proj(fcn); FEfun = FEfcn.getFun;
node = [1, 0; 1, 1; 0, 1; 1, 1/2; 1/2, 1; 1/2, 1/2]';
fun(node)
FEfun(node, FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] CR element
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0]));
CR_FES = FES(msh, CR_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.getFun;
FEfcn = CR_FES.proj(fcn); FEfun = FEfcn.getFun;
node = [1, 1/2; 1/2, 1; 1/2, 1/2]';
fun(node)
FEfun(node, FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] RT0 element
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'",...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true));
RT0_FES = FES(msh, RT0_FE);

fcn = Fcn("D2", "[1+2*x+3*y;1]");
FEfcn = RT0_FES.proj(fcn);
intFun = dot(fcn, MshEnt("D2L").UNV).getFun;
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;
intFEfcn = dot(FEfcn.subParm(ElPmSym), MshEnt("D2L").UNV);
intFEfun = intFEfcn.getFun("parm", {ElPmSym, EgPmSym});

ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
GInt("D2L", 1).eval(@(x) intFun(x, EgParm), EgParm)
GInt("D2L", 1).eval(@(x) intFEfun(x, ElParm, EgParm, ElCoef), EgParm)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
GInt("D2L", 1).eval(@(x) intFun(x, EgParm), EgParm)
GInt("D2L", 1).eval(@(x) intFEfun(x, ElParm, EgParm, ElCoef), EgParm)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
GInt("D2L", 1).eval(@(x) intFun(x, EgParm), EgParm)
GInt("D2L", 1).eval(@(x) intFEfun(x, ElParm, EgParm, ElCoef), EgParm)
%%
%[text] DG P1 element
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
DG1_FES = FES(msh, DG1_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.getFun;
FEfcn = DG1_FES.proj(fcn); FEfun = FEfcn.getFun;
node = msh.node.coord(:, msh.elem.node(:, 2));
fun(node)
FEfun(node, FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] SDG P1 element
%[text] - scalar \
SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0; 0], "EntIdx", 1), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "EntIdx", 3, "share", false)]);
SDG1_FES = FES(msh, SDG_1S_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.getFun;
FEfcn = SDG1_FES.proj(fcn); FEfun = FEfcn.getFun;
node = msh.node.coord(:, msh.elem.node(:, 2));
fun(node)
FEfun(node, FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] - vector \
SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", 1, "orien", true, "share", false)]);
SDG_1V_FES = FES(msh, SDG_1V_FE);

fcn = Fcn("D2", "[1+2*x+3*y;1]");
FEfcn = SDG_1V_FES.proj(fcn);
NFun = dot(fcn, MshEnt("D2L").UNV).getFun;
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;
NFEfcn = dot(FEfcn.subParm(ElPmSym), MshEnt("D2L").UNV);
NFEfun = NFEfcn.getFun("parm", {ElPmSym, EgPmSym});

ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NFun(EgParm, EgParm)
NFEfun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NFun(EgParm, EgParm)
NFEfun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NFun(EgParm, EgParm)
NFEfun(EgParm, ElParm, EgParm, ElCoef)
%%
%[text] - matrix \
% SDG_1M_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 2]), ...
%     [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", 1, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", [2, 3], "share", false, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UTV], "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);

load("SDG_1M_FE.mat");
SDG_1M_FES = FES(msh, SDG_1M_FE);

fcn = Fcn("D2", "[1+2*x+3*y,1;1,x+2*y+3]");
FEfcn = SDG_1M_FES.proj(fcn);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;

NNFun = dot(fcn * MshEnt("D2L").UNV, MshEnt("D2L").UNV).getFun;
NNFEfcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UNV);
NNFEfun = NNFEfcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NNFun(EgParm, EgParm)
NNFEfun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NNFun(EgParm, EgParm)
NNFEfun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NNFun(EgParm, EgParm)
NNFEfun(EgParm, ElParm, EgParm, ElCoef)

NTFun = dot(fcn * MshEnt("D2L").UNV, MshEnt("D2L").UTV).getFun;
NTFEfcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UTV);
NTFEfun = NTFEfcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NTFun(EgParm, EgParm)
NTFEfun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NTFun(EgParm, EgParm)
NTFEfun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NTFun(EgParm, EgParm)
NTFEfun(EgParm, ElParm, EgParm, ElCoef)
%%
%[text] SDG P0 element
%[text] - scalar \
SDG_0S_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0], "EntIdx", 1));
SDG_0S_FES = FES(msh, SDG_0S_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.getFun;
FEfcn = SDG_0S_FES.proj(fcn); FEfun = FEfcn.getFun;
node = [0.5, 0.5]';
fun(node)
FEfun(node, FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] - vector \
SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true));
SDG_0V_FES = FES(msh, SDG_0V_FE);

fcn = Fcn("D2", "[1+2*x+3*y;1]");
FEfcn = SDG_1V_FES.proj(fcn);
NFun = dot(fcn, MshEnt("D2L").UNV).getFun;
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;
NFEfcn = dot(FEfcn.subParm(ElPmSym), MshEnt("D2L").UNV);
NFEfun = NFEfcn.getFun("parm", {ElPmSym, EgPmSym});

ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5)); node = mean(EgParm, 2);
NFun(node, EgParm)
NFEfun(node, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4)); node = mean(EgParm, 2);
NFun(node, EgParm)
NFEfun(node, ElParm, EgParm, ElCoef)
%%
%[text] - matrix \
SDG_0M_FE = FE("D2T", FE.repFS("1", [2, 2]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", 1, "orien", true, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UTV], "orien", true, "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);
SDG_0M_FES = FES(msh, SDG_0M_FE);

fcn = Fcn("D2", "[1+2*x+3*y,1;1,x+2*y+3]");
FEfcn = SDG_1M_FES.proj(fcn);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;

NNFun = dot(fcn * MshEnt("D2L").UNV, MshEnt("D2L").UNV).getFun;
NNFEfcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UNV);
NNFEfun = NNFEfcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5)); node = mean(EgParm, 2);
NNFun(node, EgParm)
NNFEfun(node, ElParm, EgParm, ElCoef)

NTFun = dot(fcn * MshEnt("D2L").UNV, MshEnt("D2L").UTV).getFun;
NTFEfcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UTV);
NTFEfun = NTFEfcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2); ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5)); node = mean(EgParm, 2);
NTFun(node, EgParm)
NTFEfun(node, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4)); node = mean(EgParm, 2);
NTFun(node, EgParm)
NTFEfun(node, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2)); node = mean(EgParm, 2);
NTFun(node, EgParm)
NTFEfun(node, ElParm, EgParm, ElCoef)
%%
%[text] Hybrid P1 element
P1H_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
P1H_FES = FES(msh, P1H_FE);

fcn = Fcn("D2", "1+2*x+3*y"); fun = fcn.tfm("D2LR").getFun;
FEfcn = P1H_FES.proj(fcn); FEfun = FEfcn.getFun;
ElParm = FEfcn.ElParm(:, :, 5); ElCoef = FEfcn.ElCoef(:, 5);
node = [0, 1];
fun(node, ElParm)
FEfun(node, ElParm, ElCoef)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
