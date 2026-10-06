%[text] Brinkman equation
%[text] $\\begin{cases}\n  - \\nu \\Delta u + \\alpha u + \\nabla p = f &\\text{in } \\Omega \\\\\n  \\nabla \\cdot u = 0 &\\text{in } \\Omega \\\\\n  u = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  \\sigma - \\nabla u = 0 &\\text{in } \\Omega \\\\\n  -\\nu \\nabla \\cdot \\sigma + \\alpha u + \\nabla p = f &\\text{in } \\Omega \\\\\n  \\nabla \\cdot u = 0 &\\text{in } \\Omega \\\\\n  u = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
v = 1; a = 10;
u = Fcn("D2", "[sin(pi*x)^2*sin(2*pi*y); - sin(2*pi*x)*sin(pi*y)^2]");
p = Fcn("D2", "10*(x-1/2)*(y-1/2)");
d0_u = [0, 0; 0, 0]'; div_u = [1, 0; 0, 1]; dx_u = [1, 0; 1, 0]'; dy_u = [0, 1; 0, 1]'; grad_u = cat(3, dx_u, dy_u);
s = dif(u, grad_u);
d0_s = [0, 0; 0, 0; 0, 0; 0, 0]'; div_s = [1, 0; 1, 0; 0, 1; 0 ,1]';
d0_p = [0, 0]'; dx_p = [1, 0]'; dy_p = [0, 1]'; grad_p = cat(3, dx_p, dy_p);
f = - sum(dif(s, div_s), 2) * v + u * a + dif(p, grad_p);
d0_g = 0; d0_l = [0, 0];

UNV = MshEnt("D2L").UNV; UTV = MshEnt("D2L").UTV;
len = MshEnt("D2L").len;
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured triangulation + Alfeld splitting
%[text] Sub-element: triangle $\[v\_1, v\_2, c\]$, edge 1 is a primal edge, edges 2, 3 are dual edges.
msh = mshSplit(mshD2TS([0, 1, 0, 1], 16));
PrEdge = find(ismember(msh.edge.type, [0, 1, 2, 3, 4]));
PrOEdge = find(ismember(msh.edge.type, 0));
DlEdge = find(ismember(msh.edge.type, 1i));
BdEdge = find(ismember(msh.edge.type, [1, 2, 3, 4]));
%%
%[text] Scalar $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $q|\_c$ at element vertex \
%[text] Scalar $SDG\_1$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_1 (K), \\forall K \\in T\_h; \\ \[q\]\_e = 0, \\forall e \\in F^{pr, o}\_h \\}$ (pressure fixed at one DoF) \
D2TElem = MshEnt("D2T").msh;
SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", D2TElem, 1, [0, 1], d0_p, "EntIdx", 1), ...
    NdDoF("D2", D2TElem, 0, [], d0_p, "EntIdx", 3, "share", false)]);
SDG_1S_BC = BC(p, "DoF", 1);
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);
%[text] Vector $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_1$
%[text] - Nodal DoF: $v \\cdot n |\_c$ at edge vertex \
%[text] Vector $SDG\_1$ element space
%[text] - $U\_h = \\{ v \\in L^2 (\\Omega, R^2) : v|\_K \\in P\_1 (K)^2, \\forall K \\in T\_h; \\ \[v \\cdot n\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", D2TElem, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true), ...
    NdDoF("D2", D2TElem, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", 1, "orien", true, "share", false)]);
SDG_1V_FES = FES(msh, SDG_1V_FE);
%[text] Matrix $DG\_1$ element (the continuity of the SDG matrix element is imposed weakly below)
%[text] - Element: triangle
%[text] - Function space: $P^{2 \\times 2}\_1$
%[text] - Nodal DoF: each component at element vertex (not shared) \
%[text] Matrix $DG\_1$ element space
%[text] - $\\Sigma\_h = \\{ \\tau \\in L^2 (\\Omega, R^{2 \\times 2}) : \\tau|\_K \\in P\_1 (K)^{2 \\times 2}, \\forall K \\in T\_h \\}$ \
DG_1M_DoFs = NdDoF.empty;
for iComp = 1:4
    d0_si = nan(2, 4); d0_si(:, iComp) = 0;
    DG_1M_DoFs(end + 1) = NdDoF("D2", D2TElem, 0, [], d0_si, "share", false);
end
DG_1M_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 2]), DG_1M_DoFs, "map", "affine");
DG_1M_FES = FES(msh, DG_1M_FE);
%[text] Multipliers: trace spaces on edges
%[text] - Element: edge
%[text] - Function space: $P\_1$ (tangential velocity $\\gamma\_h t$ on dual edges), $P^2\_1$ (velocity $\\lambda\_h$ on interior primal edges)
%[text] - Nodal DoF: (each component) at edge vertex
%[text] - Edges where a multiplier is not used are set to zero ($\\lambda\_h = 0$ on $\\partial \\Omega$ imposes $u = 0$) \
TP_1S_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_g, "share", false));
TP_1S_FES = FES(msh, TP_1S_FE, BC(Fcn("D2", "0"), "edge", PrEdge));
TP_1V_FE = FE("D2LR", FE.repFS("[1,s]", [2, 1]), ...
    [NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_l, "coef", Fcn.cst([1; 0]), "share", false), ...
    NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_l, "coef", Fcn.cst([0; 1]), "share", false)]);
TP_1V_FES = FES(msh, TP_1V_FE, BC(Fcn("D2", "0"), "edge", [DlEdge, BdEdge]));
%%
%[text] Hybridized SDG scheme (Zhao, Chung and Lam, CMAME 364, 2020): $\\sigma\_h$ is discontinuous; $\[\\sigma\_h n\] = 0$ on interior primal edges and $\[\\sigma\_h n \\cdot t\] = 0$ on dual edges are imposed weakly by $\\lambda\_h$ and $\\gamma\_h$. Find $\\sigma\_h \\in \\Sigma\_h$, $u\_h \\in U\_h$, $\\gamma\_h$, $\\lambda\_h$, and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}\n  \\sum\_{K \\in T\_h} \\int\_K \\sigma\_h \\tau\_h \\, dx \\, dy + \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla \\cdot \\tau\_h \\, dx \\, dy - \\sum\_{e \\in F\_h^{dl}} \\int\_e (u\_h \\cdot n \[ \\tau\_h n \\cdot n \] + \\gamma\_h \[ \\tau\_h n \\cdot t \]) \\, ds - \\sum\_{e \\in F\_h^{pr, o}} \\int\_e \\lambda\_h \[ \\tau\_h n \] \\, ds = 0 &\\forall \\tau\_h \\in \\Sigma\_h \\\\\n  \\alpha \\sum\_{K \\in T\_h} \\int\_K u\_h v\_h \\, dx \\, dy - \\nu \\sum\_{K \\in T\_h} \\int\_K \\nabla \\cdot \\sigma\_h v\_h \\, dx \\, dy + \\nu \\sum\_{e \\in F\_h^{dl}} \\int\_e \[ \\sigma\_h n \\cdot n \] v\_h \\cdot n \\, ds - \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla \\cdot v\_h \\, dx \\, dy + \\sum\_{e \\in F\_h^{pr}} \\int\_e p\_h \[ v\_h \\cdot n \] \\, ds = \\sum\_{K \\in T\_h} \\int\_K f v\_h \\, dx \\, dy &\\forall v\_h \\in U\_h \\\\\n  \\nu \\sum\_{e \\in F\_h^{dl}} \\int\_e \[ \\sigma\_h n \\cdot t \] \\eta\_h \\, ds = 0, \\quad \\nu \\sum\_{e \\in F\_h^{pr, o}} \\int\_e \[ \\sigma\_h n \] \\mu\_h \\, ds = 0 &\\forall \\eta\_h, \\mu\_h \\\\\n  - \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla q\_h \\, dx \\, dy + \\sum\_{e \\in F\_h^{dl}} \\int\_e u\_h \\cdot n \[ q\_h \] \\, ds = 0 &\\forall q\_h \\in P\_h\n\\end{cases}$
Sh = DG_1M_FES; ord_Sh = 1;
Uh = SDG_1V_FES; ord_Uh = 1;
Gh = TP_1S_FES; ord_Gh = 1;
Lh = TP_1V_FES; ord_Lh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;

trls = [Sh, Uh, Gh, Lh, Ph]; tsts = [Sh, Uh, Gh, Lh, Ph];
is = 1; iu = 2; ig = 3; il = 4; ip = 5; it = 1; iv = 2; ih = 3; im = 4; iq = 5;

d0_t = d0_s; div_t = div_s; d0_v = d0_u; div_v = div_u; grad_v = grad_u; d0_h = d0_g; d0_m = d0_l; d0_q = d0_p; grad_q = grad_p;
Auv = DLF(msh, 2, Fcn.cst(a), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D2T", ord_Uh * 2));
Ast = DLF(msh, 2, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D2T", ord_Sh * 2));
But = [DLF(msh, 2, Fcn.cst(1), d0_u, div_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D2T", ord_Uh + ord_Sh - 1)), ...
    DLF.interface(msh, 1, [-UNV, UNV, UNV], d0_u, d0_t, "EntIdx", DlEdge, "iTrl", iu, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, coef(1)) .* dot(tst * coef(2), coef(3)), "GInt", GInt("D2L", ord_Uh + ord_Sh))];
Bgt = DLF.interface(msh, 1, [-UNV, UTV], d0_g, d0_t, "EntIdx", DlEdge, "iTrl", ig, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) trl .* dot(tst * coef(1), coef(2)), "GInt", GInt("D2L", ord_Gh + ord_Sh));
Blt = DLF.interface(msh, 1, -UNV, d0_l, d0_t, "EntIdx", PrOEdge, "iTrl", il, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D2L", ord_Lh + ord_Sh));
Bsv = [DLF(msh, 2, Fcn.cst(-v), div_s, d0_v, "iTrl", is, "iTst", iv, "form", @(coef, trl, tst) coef .* dot(sum(trl, 2), tst), "GInt", GInt("D2T", ord_Sh + ord_Uh - 1)), ...
    DLF.interface(msh, 1, [UNV * v, UNV, UNV], d0_s, d0_v, "EntIdx", DlEdge, "iTrl", is, "iTst", iv, "trlOpr", "jump", "form", @(coef, trl, tst) dot(tst, coef(1)) .* dot(trl * coef(2), coef(3)), "GInt", GInt("D2L", ord_Sh + ord_Uh))];
Bsh = DLF.interface(msh, 1, [UNV * v, UTV], d0_s, d0_h, "EntIdx", DlEdge, "iTrl", is, "iTst", ih, "trlOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef(1), coef(2)) .* tst, "GInt", GInt("D2L", ord_Sh + ord_Gh));
Bsm = DLF.interface(msh, 1, UNV * v, d0_s, d0_m, "EntIdx", PrOEdge, "iTrl", is, "iTst", im, "trlOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D2L", ord_Sh + ord_Lh));
Bpv = [DLF(msh, 2, Fcn.cst(-1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1)), ...
    DLF.interface(msh, 1, UNV, d0_p, d0_v, "EntIdx", PrEdge, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D2L", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 2, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh + ord_Ph - 1)), ...
    DLF.interface(msh, 1, UNV, d0_u, d0_q, "EntIdx", DlEdge, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", ord_Uh + ord_Ph))];
Fv = SLF(msh, 2, f, d0_v, "iTst", iv, "GInt", GInt("D2T", 2 + ord_Uh));
%%
%[text] Solution.
%[text] Local elimination: $\\sigma\_h$ element by element, then $(u\_h, \\gamma\_h)$ macro element by macro element (original triangles); the global system is for $(\\lambda\_h, p\_h)$.
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Ast, But, Bgt, Blt, Bsv, Bsh, Bsm, Bpv, Buq], Fv);
[sh, uh, ~, ~, ph] = FEF.multi(trls, condSolve(Stiff, Load, trls, {is, [iu, ig]}));
%%
%[text] Error.
%[text] $\\| v \\|\_{H^1, h}^2 = \\sum\_{K \\in T\_h} \\| \\nabla v \\|\_K^2 + \\sum\_{e \\in F\_h^{pr}} h^{-1}\_e \\| \[ v \] \\|\_e^2 + \\sum\_{e \\in F\_h^{dl}} h^{-1}\_e \\| \[ v \\cdot t \] \\|\_e^2$
L2Norm_Sh = Norm(msh, 2, d0_s, "GInt", GInt("D2T", (1 + ord_Sh) * 2));
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (1 + ord_Uh) * 2));
H1Norm_Uh = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", PrEdge, "coef", len \ 1, "form", @(coef, fcn, pow) coef(1) .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", DlEdge, "coef", [len \ 1, UTV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(fcn, coef(2))) .^ pow, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (1 + ord_Ph) * 2));
% Pressure is determined up to a constant (fixed at one DoF): error is measured with mean removed,
% |p - ph - c|_L2^2 = |p - ph|_L2^2 - |Omega| c^2, c = mean(p - ph), |Omega| = 1.
MeanNorm_Ph = Norm(msh, 2, d0_p, "pow", 1, "form", @(coef, fcn, pow) fcn, "GInt", GInt("D2T", (1 + ord_Ph) * 2));
eMean_Ph = eNorm(msh, p, ph, MeanNorm_Ph);
fprintf('|s-sh|_L2: %e, |u-uh|_L2: %e, |u-uh|_H1: %e, |p-ph|_L2: %e\n', eNorm(msh, s, sh, L2Norm_Sh), eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, H1Norm_Uh), sqrt(eNorm(msh, p, ph, L2Norm_Ph) ^ 2 - eMean_Ph ^ 2));
%%
%[text] Plot.
plotFEF(sh);
plotFcn(msh, s);
%%
plotFEF(uh);
plotFcn(msh, u);
%%
plotFEF(ph);
plotFcn(msh, p);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
