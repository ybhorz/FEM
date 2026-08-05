%[text] # BC - boundary condition
%%
node = Node([0, 1; 1, 1; 0, 2; 1, 2]', [-1, -2, -4, -3]);
elem = Elem([2, 3, 1; 3, 2, 4]', [5, -3, 1; -5, 4, -2]', [5, 7]);
edge = Edge([1, 2; 3, 4; 1, 3; 2, 4; 2, 3]', [1, 0; -2, 0; -1, 0; 2, 0; 1, -2]', [1, 3, 4, 2, 0]);
msh = Msh("D2T", node, elem, edge);
msh.figure;

Bdtype = [1, 2, -1, -2, -3];
BdNode = find(ismember(node.type, Bdtype));
BdEdge = find(ismember(edge.type, Bdtype));
%%
%[text] ## Continuous $P\_1$ element
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
g = Fcn("D2", "x^2+y^2"); P1_BC = BC(g, "node", BdNode);
P1_FES = FES(msh, P1_FE, P1_BC);

gFun = g.getFun;
for iBCDoF = 1:P1_FES.BC.nDoF
    iDoF = P1_FES.BC.DoFIdx(iBCDoF);
    [~, iEnt, ~] = ind2sub(P1_FES.GlDoFs, iDoF);
    iNode = P1_FES.GlDoFs.EntIdx(iEnt);
    pnt = msh.node.coord(:, iNode);
    fprintf("DoF index: %d, DoF value: %.2f, node: %d, g's value: %.2f\n", iDoF, P1_FES.BC.DoFVal(iBCDoF), iNode, gFun(pnt));
end
%%
%[text] ## Continuous $P\_2$ element
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0])]);
g = Fcn("D2", "x^3+y^3"); P2_BC =BC(g, "node", BdNode, "edge", BdEdge);
P2_FES = FES(msh, P2_FE, P2_BC);

gFun = g.getFun;
for iBCDoF = 1:P2_FES.BC.nDoF
    iDoF = P2_FES.BC.DoFIdx(iBCDoF);
    [iGlDoF, iEnt, ~] = ind2sub(P2_FES.GlDoFs, iDoF);
    switch iGlDoF
        case 1
            iNode = P2_FES.GlDoFs(1).EntIdx(iEnt);
            pnt = msh.node.coord(:, iNode);
            fprintf("DoF index: %d, DoF value: %.2f, node: %d, g's value: %.2f\n", iDoF, P2_FES.BC.DoFVal(iBCDoF), iNode, gFun(pnt));
        case 2
            iEdge = P2_FES.GlDoFs(2).EntIdx(iEnt);
            pnt = 1/2 * sum(msh.node.coord(:, msh.edge.node(:, iEdge)), 2);
            fprintf("DoF index: %d, DoF value: %.2f, midpoint of edge: %d, g's value: %.2f\n", iDoF, P2_FES.BC.DoFVal(iBCDoF), iEdge, gFun(pnt));
    end
end
%%
%[text] ## Discontinuous $P\_1$ element
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
g = Fcn("D2", "x^2+y^2"); DG1_BC = BC(g, "node", BdNode);
DG1_FES = FES(msh, DG1_FE, DG1_BC);

gFun = g.getFun;
for iBCDoF = 1:DG1_FES.BC.nDoF
    iDoF = DG1_FES.BC.DoFIdx(iBCDoF);
    [~, iEnt, ~] = ind2sub(DG1_FES.GlDoFs, iDoF);
    iNode = DG1_FES.GlDoFs.EntIdx(iEnt);
    pnt = msh.node.coord(:, iNode);
    fprintf("DoF index: %d, DoF value: %.2f, node: %d, g's value: %.2f\n", iDoF, DG1_FES.BC.DoFVal(iBCDoF), iNode, gFun(pnt));
end
%%
%[text] ## CR element
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0]));
g = Fcn("D2", "x^2+y^2"); CR_BC =BC(g, "edge", BdEdge);
CR_FES = FES(msh, CR_FE, CR_BC);

gFun = g.getFun;
for iBCDoF = 1:CR_FES.BC.nDoF
    iDoF = CR_FES.BC.DoFIdx(iBCDoF);
    [~, iEnt, ~] = ind2sub(CR_FES.GlDoFs, iDoF);
    iEdge = CR_FES.GlDoFs.EntIdx(iEnt);
    pnt = 1/2 * sum(msh.node.coord(:, msh.edge.node(:, iEdge)), 2);
    fprintf("DoF index: %d, DoF value: %.2f, midpoint of edge: %d, g's value: %.2f\n", iDoF, CR_FES.BC.DoFVal(iBCDoF), iEdge, gFun(pnt));
end
%%
%[text] ## $RT\_0$ element
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'", ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, "GInt", GInt("D2L", 4)));
g = dot(Fcn("D2", "[y^2;x^2]"), MshEnt("D2L").UNV); RT0_BC = BC(g, "edge", BdEdge);
RT0_FES = FES(msh, RT0_FE, RT0_BC);

for iBCDoF = 1:RT0_FES.BC.nDoF
    iDoF = RT0_FES.BC.DoFIdx(iBCDoF);
    [~, iEnt, ~] = ind2sub(RT0_FES.GlDoFs, iDoF);
    iEdge = RT0_FES.GlDoFs.EntIdx(iEnt);
    intFun = g.int("D2L").getFun; EgPram = msh.node.coord(:, msh.edge.node(:, iEdge));
    fprintf("DoF index: %d, DoF value: %.2f, edge: %d, g's integration: %.2f\n", iDoF, RT0_FES.BC.DoFVal(iBCDoF), iEdge, intFun([0;0], EgPram));
end
%%
%[text] ## $BDM\_1$ element
BDM1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, "GInt", GInt("D2L", 4))]);
g = dot(Fcn("D2", "[y^2;x^2]"), MshEnt("D2L").UNV); BDM1_BC = BC(g, "edge", BdEdge);
BDM1_FES = FES(msh, BDM1_FE, BDM1_BC);

for iBCDoF = 1:BDM1_FES.BC.nDoF
    iDoF = BDM1_FES.BC.DoFIdx(iBCDoF);
    [iGlDoF, iEnt, ~] = ind2sub(BDM1_FES.GlDoFs, iDoF);
    switch iGlDoF
        case 1
            iEdge = BDM1_FES.GlDoFs(1).EntIdx(iEnt);
            intFun = g.int("D2L").getFun; EgPram = msh.node.coord(:, msh.edge.node(:, iEdge));
            fprintf("DoF index: %d, DoF value: %.2f, edge: %d, g's integration 1: %.2f\n", iDoF, BDM1_FES.BC.DoFVal(iBCDoF), iEdge, intFun([0;0], EgPram));
        case 2
            iEdge = BDM1_FES.GlDoFs(2).EntIdx(iEnt);
            intFcn = g .* Fcn("D2L", "(x-x1+y-y1)/(x2-x1+y2-y1)-1/2");
            intFun = intFcn.int("D2L").getFun;
            EgPram = msh.node.coord(:, msh.edge.node(:, iEdge));
            fprintf("DoF index: %d, DoF value: %.2f, edge: %d, g's integration 2: %.2f\n", iDoF, BDM1_FES.BC.DoFVal(iBCDoF), iEdge, intFun([0;0], EgPram));
    end
end
%%
%[text] ## Vector $P\_1$ element
VP1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]);
g = [Fcn("D2", "x^2+y^2"), Fcn("D2", "x^2-y^2")]; VP1_BC = BC(g, "node", BdNode);
VP1_FES = FES(msh, VP1_FE, VP1_BC);

g1Fun = g(1).getFun; g2Fun = g(2).getFun;
for iBCDoF = 1:VP1_FES.BC.nDoF
    iDoF = VP1_FES.BC.DoFIdx(iBCDoF);
    [iGlDoF, iEnt, ~] = ind2sub(VP1_FES.GlDoFs, iDoF);
    switch iGlDoF
        case 1
            iNode = VP1_FES.GlDoFs(1).EntIdx(iEnt);
            pnt = msh.node.coord(:, iNode);
            fprintf("DoF index: %d, DoF value: %.2f, node: %d, g's value: %.2f\n", iDoF, VP1_FES.BC.DoFVal(iBCDoF), iNode, g1Fun(pnt));
        case 2
            iNode = VP1_FES.GlDoFs(2).EntIdx(iEnt);
            pnt = msh.node.coord(:, iNode);
            fprintf("DoF index: %d, DoF value: %.2f, node: %d, g's value: %.2f\n", iDoF, VP1_FES.BC.DoFVal(iBCDoF), iNode, g2Fun(pnt));
    end
end
%%
%[text] ## Trace P1 element
TP1_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
g = Fcn("D2", "x^2+y^2"); TP1_BC = BC(g, "edge", BdEdge);
TP1_FES = FES(msh, TP1_FE, TP1_BC);

gFun = g.getFun;
for iBCDoF = 1:TP1_FES.BC.nDoF
    iDoF = TP1_FES.BC.DoFIdx(iBCDoF);
    [~, iEnt, iSamp] = ind2sub(TP1_FES.GlDoFs, iDoF);
    iEdge = TP1_FES.GlDoFs.EntIdx(iEnt);
    switch iSamp
        case 1
            pnt = msh.node.coord(:, msh.edge.node(1, iEdge));
            fprintf("DoF index: %d, DoF value: %.2f, vertex 1 of edge: %d, g's value: %.2f\n", iDoF, TP1_FES.BC.DoFVal(iBCDoF), iEdge, gFun(pnt));
        case 2
            pnt = msh.node.coord(:, msh.edge.node(2, iEdge));
            fprintf("DoF index: %d, DoF value: %.2f, vertex 2 of edge: %d, g's value: %.2f\n", iDoF, TP1_FES.BC.DoFVal(iBCDoF), iEdge, gFun(pnt));
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
