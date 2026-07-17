clear;
node = Node([0, 0; 1, 0; 0, 1; 1, 1]', [-1, -2, -4, -3]);
elem = Elem([2, 3, 1; 3, 2, 4]', [5, -3, 1; -5, 4, -2]', [5, 7]);
edge = Edge([1, 2; 3, 4; 1, 3; 2, 4; 2, 3]', [1, 0; -2, 0; -1, 0; 2, 0; 1, -2]', [1, 3, 4, 2, 0]);
msh = Msh("D2T", node, elem, edge);
msh.figure;

Bdtype = [2, 3, 0, -3, -4];
BdNode = find(ismember(node.type, Bdtype));
BdEdge = find(ismember(edge.type, Bdtype));
%%
%[text] Lagrange element P1
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
P1_BC = BC(Fcn("D2", "x+y"), "node", BdNode);
P1_FES = FES(msh, P1_FE, P1_BC);

P1_FES.GlDoFs.msh.type
P1_FES.GlDoFs.EntIdx
P1_FES.Lc2Gl
P1_FES.BC.DoFIdx
P1_FES.BC.DoFVal

DoFVal = zeros(P1_FES.nGlDoF,1); DoFVal(3) = 1;
FEfcn = FEF(P1_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun(msh.node.coord(:, msh.elem.node(:, 2)), FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] Lagrange element P2
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0])]);
P2_BC = BC(Fcn("D2", "x+y"), "node", BdNode, "edge", BdEdge);
P2_FES = FES(msh, P2_FE, P2_BC);

P2_FES.GlDoFs(1).msh.type
P2_FES.GlDoFs(1).EntIdx

P2_FES.GlDoFs(2).msh.type
P2_FES.GlDoFs(2).EntIdx

P2_FES.Lc2Gl

P2_FES.BC.DoFIdx
P2_FES.BC.DoFVal

DoFVal = zeros(P2_FES.nGlDoF,1); DoFVal(3) = 1;
FEfcn = FEF(P2_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun([0, 1; 0.5, 0.5; 1, 0]', FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))

DoFVal = zeros(P2_FES.nGlDoF,1); DoFVal(9) = 1;
FEfcn = FEF(P2_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun([0, 1; 0.5, 0.5; 1, 0]', FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] CR element
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0]));
CR_BC = BC(Fcn("D2", "x+y"), "edge", BdEdge);
CR_FES = FES(msh, CR_FE, CR_BC);

CR_FES.GlDoFs.msh.type
CR_FES.GlDoFs.EntIdx
CR_FES.Lc2Gl
CR_FES.BC.DoFIdx
CR_FES.BC.DoFVal

DoFVal = zeros(CR_FES.nGlDoF,1); DoFVal(5) = 1;
FEfcn = FEF(CR_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun([0.5, 0.5; 1, 0.5; 0.5, 1]', FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] RT0 element
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'",...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true));
RT0_BC = BC(dot(Fcn("D2", "[y;x]"), MshEnt("D2L").UNV), "edge", BdEdge);
RT0_FES = FES(msh, RT0_FE, RT0_BC);

RT0_FES.GlDoFs.msh.type
RT0_FES.GlDoFs.EntIdx
RT0_FES.Lc2Gl
RT0_FES.BC.DoFIdx
RT0_FES.BC.DoFVal

DoFVal = zeros(RT0_FES.nGlDoF,1); DoFVal(5) = 1;
FEfcn = FEF(RT0_FES, DoFVal);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;
intFcn = dot(FEfcn.subParm(ElPmSym), MshEnt("D2L").UNV);
intFun = intFcn.getFun("parm", {ElPmSym, EgPmSym});

ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
GInt("D2L", 1).eval(@(x) intFun(x, ElParm, EgParm, ElCoef), EgParm)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
GInt("D2L", 1).eval(@(x) intFun(x, ElParm, EgParm, ElCoef), EgParm)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
GInt("D2L", 1).eval(@(x) intFun(x, ElParm, EgParm, ElCoef), EgParm)
%%
%[text] DG P1 element
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
DG1_BC = BC(Fcn("D2", "x+y"), "node", BdNode);
DG1_FES = FES(msh, DG1_FE, DG1_BC);

DG1_FES.GlDoFs.msh.type
DG1_FES.GlDoFs.EntIdx
DG1_FES.Lc2Gl
DG1_FES.BC.DoFIdx
DG1_FES.BC.DoFVal

DoFVal = zeros(DG1_FES.nGlDoF,1); DoFVal(4) = 1;
FEfcn = FEF(DG1_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun(msh.node.coord(:, msh.elem.node(:, 2)), FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] SDG P1 element
%[text] - scalar \
SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0; 0], "EntIdx", 1), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "EntIdx", 3, "share", false)]);
SDG_1S_BC = BC(Fcn("D2", "x+y"), "node", BdNode, "edge", BdEdge);
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);

SDG_1S_FES.GlDoFs(1).msh.type
SDG_1S_FES.GlDoFs(1).EntIdx

SDG_1S_FES.GlDoFs(2).msh.type
SDG_1S_FES.GlDoFs(2).EntIdx

SDG_1S_FES.Lc2Gl
SDG_1S_FES.BC.DoFIdx
SDG_1S_FES.BC.DoFVal

DoFVal = zeros(SDG_1S_FES.nGlDoF,1); DoFVal(1) = 1;
FEfcn = FEF(SDG_1S_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun(msh.node.coord(:, msh.elem.node(:, 2)), FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] - vector \
SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", 1, "orien", true, "share", false)]);
SDG_1V_BC = BC(dot(Fcn("D2", "[y;x]"), MshEnt("D2L").UNV), "edge", BdEdge);
SDG_1V_FES = FES(msh, SDG_1V_FE, SDG_1V_BC);

SDG_1V_FES.GlDoFs(1).msh.type
SDG_1V_FES.GlDoFs(1).EntIdx

SDG_1V_FES.GlDoFs(2).msh.type
SDG_1V_FES.GlDoFs(2).EntIdx

SDG_1V_FES.Lc2Gl
SDG_1V_FES.BC.DoFIdx
SDG_1V_FES.BC.DoFVal

DoFVal = zeros(SDG_1V_FES.nGlDoF,1); DoFVal(4) = 1;
FEfcn = FEF(SDG_1V_FES, DoFVal);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;
NFcn = dot(FEfcn.subParm(ElPmSym), MshEnt("D2L").UNV);
NFun = NFcn.getFun("parm", {ElPmSym, EgPmSym});

ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NFun(EgParm, ElParm, EgParm, ElCoef)
%%
%[text] - matrix \
% SDG_1M_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 2]), ...
%     [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", 1, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", [2, 3], "share", false, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UTV], "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);

load("SDG_1M_FE.mat");

real = Fcn("D2", "[y,x;x,y]");
BdFcn = [dot(real * MshEnt("D2L").UNV, MshEnt("D2L").UNV), ...
    dot(real * MshEnt("D2L").UNV, MshEnt("D2L").UNV), ...
    dot(real * MshEnt("D2L").UNV, MshEnt("D2L").UTV)];
SDG_1M_BC = BC(BdFcn, "edge", BdEdge);
SDG_1M_FES = FES(msh, SDG_1M_FE, SDG_1M_BC);

SDG_1M_FES.GlDoFs(1).msh.type
SDG_1M_FES.GlDoFs(1).EntIdx

SDG_1M_FES.GlDoFs(2).msh.type
SDG_1M_FES.GlDoFs(2).EntIdx

SDG_1M_FES.GlDoFs(3).msh.type
SDG_1M_FES.GlDoFs(3).EntIdx

SDG_1M_FES.Lc2Gl
SDG_1M_FES.BC.DoFIdx
SDG_1M_FES.BC.DoFVal

DoFVal = zeros(SDG_1M_FES.nGlDoF,1); DoFVal(2) = 1;
FEfcn = FEF(SDG_1M_FES, DoFVal);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;

NNFcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UNV);
NNFun = NNFcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NNFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NNFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NNFun(EgParm, ElParm, EgParm, ElCoef)

NTFcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UTV);
NTFun = NTFcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NTFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NTFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NTFun(EgParm, ElParm, EgParm, ElCoef)
%%
%[text] SDG P0 element
%[text] - scalar \
SDG_0S_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0], "EntIdx", 1));
SDG_0S_BC = BC(Fcn("D2", "x+y"), "node", BdNode, "edge", BdEdge);
SDG_0S_FES = FES(msh, SDG_0S_FE, SDG_0S_BC);

SDG_0S_FES.GlDoFs(1).msh.type
SDG_0S_FES.GlDoFs(1).EntIdx

SDG_0S_FES.Lc2Gl
SDG_0S_FES.BC.DoFIdx
SDG_0S_FES.BC.DoFVal

DoFVal = zeros(SDG_0S_FES.nGlDoF,1); DoFVal(1) = 1;
FEfcn = FEF(SDG_0S_FES, DoFVal); FEfun = FEfcn.getFun;
FEfun(msh.node.coord(:, msh.elem.node(:, 2)), FEfcn.ElParm(:, :, 2), FEfcn.ElCoef(:, 2))
%%
%[text] - vector \
SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0,0;0,0], "coef", MshEnt("D2L").UNV, "EntIdx", [2, 3], "orien", true));
SDG_0V_BC = BC(dot(Fcn("D2", "[y;x]"), MshEnt("D2L").UNV), "edge", BdEdge);
SDG_0V_FES = FES(msh, SDG_0V_FE, SDG_0V_BC);

SDG_0V_FES.GlDoFs(1).msh.type
SDG_0V_FES.GlDoFs(1).EntIdx

SDG_0V_FES.Lc2Gl
SDG_0V_FES.BC.DoFIdx
SDG_0V_FES.BC.DoFVal

DoFVal = zeros(SDG_0V_FES.nGlDoF,1); DoFVal(4) = 1;
FEfcn = FEF(SDG_0V_FES, DoFVal);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;
NFcn = dot(FEfcn.subParm(ElPmSym), MshEnt("D2L").UNV);
NFun = NFcn.getFun("parm", {ElPmSym, EgPmSym});

ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NFun(EgParm, ElParm, EgParm, ElCoef)
%%
%[text] - matrix \
SDG_0M_FE = FE("D2T", FE.repFS("1", [2, 2]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UNV], "EntIdx", 1, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, zeros(2,4), "coef", [MshEnt("D2L").UNV, MshEnt("D2L").UTV], "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);
real = Fcn("D2", "[y,x;x,y]");
BdFcn = [dot(real * MshEnt("D2L").UNV, MshEnt("D2L").UNV), ...
    dot(real * MshEnt("D2L").UNV, MshEnt("D2L").UTV)];
SDG_0M_BC = BC(BdFcn, "edge", BdEdge);
SDG_0M_FES = FES(msh, SDG_0M_FE, SDG_0M_BC);

SDG_0M_FES.GlDoFs(1).msh.type
SDG_0M_FES.GlDoFs(1).EntIdx

SDG_0M_FES.GlDoFs(2).msh.type
SDG_0M_FES.GlDoFs(2).EntIdx

SDG_0M_FES.Lc2Gl
SDG_0M_FES.BC.DoFIdx
SDG_0M_FES.BC.DoFVal

DoFVal = zeros(SDG_0M_FES.nGlDoF,1); DoFVal(1) = 1;
FEfcn = FEF(SDG_0M_FES, DoFVal);
ElPmSym = sym("ElParm", FEfcn.sElParm); EgPmSym = MshEnt("D2L").parm;

NNFcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UNV);
NNFun = NNFcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NNFun(EgParm, ElParm, EgParm, ElCoef)

NTFcn = dot(FEfcn.subParm(ElPmSym) * MshEnt("D2L").UNV, MshEnt("D2L").UTV);
NTFun = NTFcn.getFun("parm", {ElPmSym, EgPmSym});
ElParm = FEfcn.ElParm(:, :, 2);
ElCoef = FEfcn.ElCoef(:, 2);
EgParm = msh.node.coord(:, msh.edge.node(:, 5));
NTFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 4));
NTFun(EgParm, ElParm, EgParm, ElCoef)
EgParm = msh.node.coord(:, msh.edge.node(:, 2));
NTFun(EgParm, ElParm, EgParm, ElCoef)
%%
%[text] Hybrid P1 element
P1H_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
P1H_BC = BC(Fcn("D2", "x+y"), "edge", BdEdge);
P1H_FES = FES(msh, P1H_FE, P1H_BC);

P1H_FES.GlDoFs.msh.type
P1H_FES.GlDoFs.EntIdx
P1H_FES.Lc2Gl
P1H_FES.BC.DoFIdx
P1H_FES.BC.DoFVal

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
