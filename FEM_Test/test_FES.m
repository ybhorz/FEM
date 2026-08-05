%[text] # FES - finite element space
%%
node = Node([0, 1; 1, 1; 0, 2; 1, 2]', [-1, -2, -4, -3]);
elem = Elem([2, 3, 1; 3, 2, 4]', [5, -3, 1; -5, 4, -2]', [5, 7]);
edge = Edge([1, 2; 3, 4; 1, 3; 2, 4; 2, 3]', [1, 0; -2, 0; -1, 0; 2, 0; 1, -2]', [1, 3, 4, 2, 0]);
msh = Msh("D2T", node, elem, edge);
msh.figure;
%%
%[text] ## $P\_0$ element space
%[text] $P\_0$ element
%[text] - Domain: triangle
%[text] - Function space: $P\_0$
%[text] - Nodal DoF: $v|\_c$ at centroid \
%[text] $P\_0$ element space
%[text] $U\_h = \\{ v \\in L^2 (\\Omega) : v|\_{K} \\in P\_0 (K), \\forall K \\in T\_h \\}$
P0_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0]));
P0_FES = FES(msh, P0_FE);

disp(P0_FES.GlDoFs.EntDim);
disp(P0_FES.GlDoFs.EntIdx);
disp(P0_FES.Lc2Gl);
%[text] $P\_0$ base function
iDoF = 1;
BsDoF = zeros(P0_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(P0_FES, BsDoF); BsFun = BsFcn.getFun;

for iElem = 1:msh.nElem
    cPnt = 1/3 * sum(msh.node.coord(:, msh.elem.node(:, iElem)), 2);
    fprintf("centroid of element: %d, base's value: %.2f\n", iElem, ...
        BsFun(cPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
end
%[text] $P\_0$ projection
fcn = Fcn("D2", "x+y"); fun = fcn.getFun;
PFcn = P0_FES.proj(fcn); PFun = PFcn.getFun;

for iElem = 1:msh.nElem
    cPnt = 1/3 * sum(msh.node.coord(:, msh.elem.node(:, iElem)), 2);
    fprintf("centroid of element: %d, function's value: %.2f, projection's value: %.2f\n", iElem, ...
        fun(cPnt), PFun(cPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
end
%%
%[text] ## Continuous $P\_1$ element space
%[text] $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Continuous $P\_1$ element space
%[text] - $U\_h = \\{ v \\in H^1 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h \\}$ \
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
P1_FES = FES(msh, P1_FE);

disp(P1_FES.GlDoFs.EntDim);
disp(P1_FES.GlDoFs.EntIdx);
disp(P1_FES.Lc2Gl);
%[text] $P\_1$ base function
iDoF = 1;
BsDoF = zeros(P1_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(P1_FES, BsDoF); BsFun = BsFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, base's value: %.2f\n", iElem, iNode, ...
            BsFun(vPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] $P\_1$ projection
fcn = Fcn("D2", "x^2+y^2"); fun = fcn.getFun;
PFcn = P1_FES.proj(fcn); PFun = PFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, function's value: %.2f, projection's value: %.2f\n", iElem, iNode, ...
            fun(vPnt), PFun(vPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## Continuous $P\_2$ element space
%[text] $P\_2$ element
%[text] - Domain: triangle
%[text] - Function space: $P\_2$
%[text] - Nodal DoF: $v|\_c$ at element vertex and edge midpoint \
%[text] Continuous $P\_2$ element space
%[text] - $U\_h = \\{ v \\in H^1 (\\Omega) : v|\_{K} \\in P\_2 (K), \\forall K \\in T\_h\\}$ \
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0])]);
P2_FES = FES(msh, P2_FE);

disp(P2_FES.GlDoFs(1).EntDim);
disp(P2_FES.GlDoFs(1).EntIdx);
disp(P2_FES.GlDoFs(2).EntDim)
disp(P2_FES.GlDoFs(2).EntIdx);
disp(P2_FES.Lc2Gl);
%[text] $P\_2$ base function
iDoF = 1;
BsDoF = zeros(P2_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(P2_FES, BsDoF); BsFun = BsFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, base's value: %.2f\n", iElem, iNode, ...
            BsFun(vPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        mPnt = 1/2 * sum(msh.node.coord(:, msh.edge.node(:, iEdge)), 2);
        fprintf("element: %d, midpoint of edge: %d, base's value: %.2f\n", iElem, iEdge, ...
            BsFun(mPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] $P\_2$ projection
fcn = Fcn("D2", "x^3+y^3"); fun = fcn.getFun;
PFcn = P2_FES.proj(fcn); PFun = PFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, function's value: %.2f, projection's value: %.2f\n", iElem, iNode, ...
            fun(vPnt), PFun(vPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        mPnt = 1/2 * sum(msh.node.coord(:, msh.edge.node(:, iEdge)), 2);
        fprintf("element: %d, midpoint of edge: %d, function's value: %.2f, projection's value: %.2f\n", iElem, iEdge, ...
            fun(mPnt), PFun(mPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## Discontinuous $P\_1$ element space
%[text] $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Discontinuous $P\_1$ element space
%[text] - $U\_h = \\{ v \\in L^2 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h \\}$ \
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
DG1_FES = FES(msh, DG1_FE);

disp(DG1_FES.GlDoFs.EntDim);
disp(DG1_FES.GlDoFs.EntIdx);
disp(DG1_FES.Lc2Gl);
%[text] $P\_1$ base function
iDoF = 1;
BsDoF = zeros(DG1_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(DG1_FES, BsDoF); BsFun = BsFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, base value: %.2f\n", iElem, iNode, ...
            BsFun(vPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] $P\_1$ projection
fcn = Fcn("D2", "x^2+y^2"); fun = fcn.getFun;
PFcn = DG1_FES.proj(fcn); PFun = PFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, function's value: %.2f, projection's value: %.2f\n", iElem, iNode, ...
            fun(vPnt), PFun(vPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## CR element space
%[text] CR element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at edge midpoint \
%[text] CR element space
%[text] - $U\_h = \\{ v \\in L^2 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h; \\ \[v\]\_{c(e)} = 0, \\forall e \\in E^o\_h \\}$ \
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0]));
CR_FES = FES(msh, CR_FE);

disp(CR_FES.GlDoFs.EntDim);
disp(CR_FES.GlDoFs.EntIdx);
disp(CR_FES.Lc2Gl);
%[text] CR base function
iDoF = 1;
BsDoF = zeros(CR_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(CR_FES, BsDoF); BsFun = BsFcn.getFun;

for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        mPnt = 1/2 * sum(msh.node.coord(:, msh.edge.node(:, iEdge)), 2);
        fprintf("element: %d, midpoint of edge: %d, base's value: %.2f\n", iElem, iEdge, ...
            BsFun(mPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] CR projection
fcn = Fcn("D2", "x^2+y^2"); fun = fcn.getFun;
PFcn = CR_FES.proj(fcn); PFun = PFcn.getFun;

for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        mPnt = 1/2 * sum(msh.node.coord(:, msh.edge.node(:, iEdge)), 2);
        fprintf("element: %d, midpoint of edge: %d, function's value: %.2f, projection's value: %.2f\n", iElem, iEdge, ...
            fun(mPnt), PFun(mPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## $RT\_0$ element space
%[text] $RT\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_0 + x P\_0$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ on edge \
%[text] $RT\_0$ element space
%[text] - $U\_h  = \\{ v \\in H(div; \\Omega) : v|\_K \\in RT\_0 (K), \\forall K \\in T\_h \\}$ \
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'", ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, "GInt", GInt("D2L", 4)));
RT0_FES = FES(msh, RT0_FE);

disp(RT0_FES.GlDoFs.EntDim);
disp(RT0_FES.GlDoFs.EntIdx);
disp(RT0_FES.Lc2Gl);
%[text] $RT\_0$ base function
iDoF = 1;
BsDoF = zeros(RT0_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(RT0_FES, BsDoF);

for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        intFun = dot(BsFcn, MshEnt("D2T").UNV(iElEg)).int("D2L", iElEg).getFun;
        fprintf("element: %d, edge: %d, base's integration: %.2f\n", iElem, iEdge, ...
            intFun([0;0], BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] $RT\_0$ projection
fcn = Fcn("D2", "[x^2;y^2]");
PFcn = RT0_FES.proj(fcn);

for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        intFun = dot(fcn, MshEnt("D2T").UNV(iElEg)).int("D2L", iElEg).getFun;
        intPFun = dot(PFcn, MshEnt("D2T").UNV(iElEg)).int("D2L", iElEg).getFun;
        fprintf("element: %d, edge: %d, function's integration: %.2f, projection's integration: %.2f\n", iElem, iEdge, ...
            intFun([0;0], PFcn.ElParm(:, :, iElem)), ...
            intPFun([0;0], PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## $BDM\_1$ element space
%[text] $BDM\_1$ element
%[text] - Domain: triangle
%[text] - Function space: $P^2\_1$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ and $\\int\_e (v \\cdot n) (s - 1/2) \\, ds$ on edge \
%[text] $BDM\_1$ element space
%[text] - $U\_h  = \\{ v \\in H(div; \\Omega) : v|\_K \\in BDM\_1 (K), \\forall K \\in T\_h \\}$ \
BDM1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true, "GInt", GInt("D2L", 4))]);
BDM1_FES = FES(msh, BDM1_FE);

disp(BDM1_FES.GlDoFs(1).EntDim);
disp(BDM1_FES.GlDoFs(1).EntIdx);
disp(BDM1_FES.GlDoFs(2).EntDim);
disp(BDM1_FES.GlDoFs(2).EntIdx);
disp(BDM1_FES.Lc2Gl);
%[text] $BDM\_1$ base function
iDoF = 1;
BsDoF = zeros(BDM1_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(BDM1_FES, BsDoF);

for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        intFun = dot(BsFcn, MshEnt("D2T").UNV(iElEg)).int("D2L", iElEg).getFun;
        fprintf("element: %d, edge: %d, base's integration 1: %.2f\n", iElem, iEdge, ...
            intFun([0;0], BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end

tsts = [Fcn("D2T", "(x-x1+y-y1)/(x2-x1+y2-y1)-1/2"); ...
    Fcn("D2T", "(x-x2+y-y2)/(x3-x2+y3-y2)-1/2"); ...
    Fcn("D2T", "(x-x3+y-y3)/(x1-x3+y1-y3)-1/2")];
for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        intFcn = dot(BsFcn, MshEnt("D2T").UNV(iElEg)) .* tsts(iElEg);
        intFun = intFcn.int("D2L", iElEg).getFun;
        fprintf("element: %d, edge: %d, base's integration 2: %.2f\n", iElem, iEdge, ...
            intFun([0;0], BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] $BDM\_1$ projection
fcn = Fcn("D2", "[x^2;y^2]");
PFcn = BDM1_FES.proj(fcn);

for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        intFun = dot(fcn, MshEnt("D2T").UNV(iElEg)).int("D2L", iElEg).getFun;
        intPFun = dot(PFcn, MshEnt("D2T").UNV(iElEg)).int("D2L", iElEg).getFun;
        fprintf("element: %d, edge: %d, function's integration 1: %.2f, projection's integration 1: %.2f\n", iElem, iEdge, ...
            intFun([0;0], PFcn.ElParm(:, :, iElem)), ...
            intPFun([0;0], PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end

tsts = [Fcn("D2T", "(x-x1+y-y1)/(x2-x1+y2-y1)-1/2"); ...
    Fcn("D2T", "(x-x2+y-y2)/(x3-x2+y3-y2)-1/2"); ...
    Fcn("D2T", "(x-x3+y-y3)/(x1-x3+y1-y3)-1/2")];
for iElem = 1:msh.nElem
    for iElEg = 1:msh.elem.nEdge
        iEdge = abs(msh.elem.edge(iElEg, iElem));
        intFcn = dot(fcn, MshEnt("D2T").UNV(iElEg)) .* tsts(iElEg);
        intFun = intFcn.int("D2L", iElEg).getFun;
        intPFcn = dot(PFcn, MshEnt("D2T").UNV(iElEg)) .* tsts(iElEg);
        intPFun = intPFcn.int("D2L", iElEg).getFun;
        fprintf("element: %d, edge: %d, function's integration 2: %.2f, projection's integration 2: %.2f\n", iElem, iEdge, ...
            intFun([0;0], PFcn.ElParm(:, :, iElem)), ...
            intPFun([0;0], PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## Vector $P\_1$ element
%[text] Vector $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1^2$
%[text] - Nodal DoF: $v\_x \\, |\_c$ and $v\_y \\, |\_c$ at element vertex \
%[text] Vector $P\_1$ element space
%[text] - $U\_h = \\{ v \\in H^1 (\\Omega) : v|\_{K} \\in P\_1 (K)^2, \\forall K \\in T\_h \\}$ \
VP1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]);
VP1_FES = FES(msh, VP1_FE);

disp(VP1_FES.GlDoFs(1).EntDim);
disp(VP1_FES.GlDoFs(1).EntIdx);
disp(VP1_FES.GlDoFs(2).EntDim);
disp(VP1_FES.GlDoFs(2).EntIdx);
disp(VP1_FES.Lc2Gl);
%[text] $P\_1$ base function
iDoF = 1;
BsDoF = zeros(VP1_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(VP1_FES, BsDoF); BsFun = BsFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, base's value: [%.2f;%.2f]\n", iElem, iNode, ...
            BsFun(vPnt, BsFcn.ElParm(:, :, iElem), BsFcn.ElCoef(:, iElem)));
    end
end
%[text] Vector $P\_1$ projection
fcn = Fcn("D2", "[x^2;y^2]"); fun = fcn.getFun;
PFcn = VP1_FES.proj(fcn); PFun = PFcn.getFun;

for iElem = 1:msh.nElem
    for iElNd = 1:msh.elem.nNode
        iNode = msh.elem.node(iElNd, iElem);
        vPnt = msh.node.coord(:, iNode);
        fprintf("element: %d, node: %d, function's value: [%.2f;%.2f], projection's value: [%.2f;%.2f]\n", iElem, iNode, ...
            fun(vPnt), PFun(vPnt, PFcn.ElParm(:, :, iElem), PFcn.ElCoef(:, iElem)));
    end
end
%%
%[text] ## Trace P1 element space
%[text] Trace P1 element
%[text] - Element: reference edge \[0, 1\]
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $\\hat{\\mu} |\_c$ at edge vertex \
%[text] Trace $P\_1$ element space
%[text] - $\\Lambda\_h = \\{ \\mu \\in L^2 (E\_h) : \\mu |\_e \\in P\_1 (e), \\forall e \\in E\_h \\}$ \
TP1_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
TP1_FES = FES(msh, TP1_FE);

disp(TP1_FES.GlDoFs.EntDim);
disp(TP1_FES.GlDoFs.EntIdx);
disp(TP1_FES.Lc2Gl);
%[text] Trace $P\_1$ base function
iDoF = 1;
BsDoF = zeros(TP1_FES.nGlDoF, 1); BsDoF(iDoF) = 1;
BsFcn = FEF(TP1_FES, BsDoF); BsFun = BsFcn.getFun;

for iElEg = 1:msh.nEdge
    for iElNd = 1:msh.edge.nNode
        iNode = msh.edge.node(iElNd, iElEg);
        rPnt = iElNd - 1;
        fprintf("edge: %d, node: %d, base's value: %.2f\n", iElEg, iNode, ...
                BsFun(rPnt, BsFcn.ElParm(:, :, iElEg), BsFcn.ElCoef(:, iElEg)));
    end
end
%[text] Trace $P\_1$ projection
fcn = Fcn("D2", "x^2+y^2"); fun = fcn.getFun;
PFcn = TP1_FES.proj(fcn); PFun = PFcn.getFun;

for iElEg = 1:msh.nEdge
    for iElNd = 1:msh.edge.nNode
        iNode = msh.edge.node(iElNd, iElEg);
        rPnt = iElNd - 1;
        vPnt = msh.node.coord(:, iNode);
        fprintf("edge: %d, node: %d, function's value: %.2f, function's value: %.2f\n", iElEg, iNode, ...
                fun(vPnt), BsFun(rPnt, PFcn.ElParm(:, :, iElEg), PFcn.ElCoef(:, iElEg)));
    end
end
%%
%[text] ## FES.Method
%[text] ### getMsh, cumDoF
FESs = [VP1_FES, RT0_FES];

disp(FESs.getMsh.type); 

disp(FESs.cumDoF);
disp(FESs.cumDoF(1));
disp(FESs.cumDoF(2));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
