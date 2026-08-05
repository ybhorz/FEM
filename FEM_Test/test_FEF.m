%[text] # FEF - finite element function
%%
node = Node([0, 1; 1, 1; 0, 2; 1, 2]', [-1, -2, -4, -3]);
elem = Elem([2, 3, 1; 3, 2, 4]', [5, -3, 1; -5, 4, -2]', [5, 7]);
edge = Edge([1, 2; 3, 4; 1, 3; 2, 4; 2, 3]', [1, 0; -2, 0; -1, 0; 2, 0; 1, -2]', [1, 3, 4, 2, 0]);
msh = Msh("D2T", node, elem, edge);
msh.figure;
%%
%[text] ## $P\_0$ element function
P0_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0]));
P0_FES = FES(msh, P0_FE);

FEFcn = FEF(P0_FES, [0, 1]);
plotFEF(FEFcn);
%%
%[text] ## Continuous $P\_1$ element function
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
P1_FES = FES(msh, P1_FE);

FEFcn = FEF(P1_FES, [0, 1, 3, 2]);
plotFEF(FEFcn);
%%
%[text] ## Discontinuous $P\_1$ element function
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
DG1_FES = FES(msh, DG1_FE);

FEFcn = FEF(DG1_FES, [1, 2, 0, 3, 2, 4]);
plotFEF(FEFcn);
%%
%[text] ## CR element function
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0]));
CR_FES = FES(msh, CR_FE);

FEFcn = FEF(CR_FES, [0, 1, 2, 3, 4]);
plotFEF(FEFcn);
%%
%[text] ## $RT\_0$ element function
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'", ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true));
RT0_FES = FES(msh, RT0_FE);

FEFcn = FEF(RT0_FES, [0, 1, 2, 3, 4]);
plotFEF(FEFcn);
%%
%[text] ## $BDM\_1$ element function
BDM1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true)]);
BDM1_FES = FES(msh, BDM1_FE);

FEFcn = FEF(BDM1_FES, [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]);
plotFEF(FEFcn);
%%
%[text] ## Vector $P\_1$ element function
VP1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]);
VP1_FES = FES(msh, VP1_FE);

FEFcn = FEF(VP1_FES, [0, 1, 3, 2, 0, 3, 1, 2]);
plotFEF(FEFcn);
%%
%[text] ## FEF.Method
%[text] ### subt
FEFcn1 = FEF(P1_FES, [0, 1, 3, 2]);
FEFcn2 = FEF(P1_FES, [0, 3, 1, 2]);
plotFEF(subt(FEFcn1, FEFcn2));
%%
%[text] ## multi
FESs = [P1_FES, CR_FES];

DoFVal = [0, 1, 3, 2, 0, 1, 2, 3, 4];

[FEFcn1, FEFcn2] = FEF.multi(FESs, DoFVal);
plotFEF(FEFcn1); plotFEF(FEFcn2);

FEFs = FEF.multi(FESs, DoFVal);
plotFEF(FEFs(1)); plotFEF(FEFs(2));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
