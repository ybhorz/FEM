%[text] Stokes equation
%[text] $\\begin{cases}\n  - \\nu \\Delta u + \\nabla p = f &\\text{in } \\Omega \\\\\n  \\nabla \\cdot u = 0 &\\text{in } \\Omega \\\\\n  u = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  \\sigma - \\nabla u = 0 &\\text{in } \\Omega \\\\\n  -\\nu \\nabla \\cdot \\sigma + \\nabla p = f &\\text{in } \\Omega \\\\\n  \\nabla \\cdot u = 0 &\\text{in } \\Omega \\\\\n  u = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
v = 1;
u = Fcn("D2", "[sin(pi*x)^2*sin(2*pi*y); - sin(2*pi*x)*sin(pi*y)^2]");
p = Fcn("D2", "10*(x-1/2)*(y-1/2)");
d0_u = [0, 0; 0, 0]'; div_u = [1, 0; 0, 1]; dx_u = [1, 0; 1, 0]'; dy_u = [0, 1; 0, 1]'; grad_u = cat(3, dx_u, dy_u);
s = dif(u, grad_u);
d0_s = [0, 0; 0, 0; 0, 0; 0, 0]'; div_s = [1, 0; 1, 0; 0, 1; 0 ,1]';
d0_p = [0, 0]'; dx_p = [1, 0]'; dy_p = [0, 1]'; grad_p = cat(3, dx_p, dy_p);
f = - sum(dif(s, div_s), 2) * v + dif(p, grad_p);

UNV = MshEnt("D2L").UNV;
len = MshEnt("D2L").len;
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured triangulation + Alfeld splitting
msh = mshSplit(mshD2TS([0, 1, 0, 1], 16));
PrEdge = find(ismember(msh.edge.type, [0, 1, 2, 3, 4]));
DlEdge = find(ismember(msh.edge.type, 1i));
BdEdge = find(ismember(msh.edge.type, [1, 2, 3, 4]));
BdNode = find(ismember(msh.node.type, [1, 2, 3, 4, -1, -2, -3, -4]));
DlNode = find(ismember(msh.node.type, 1i));
%%
%[text] Scalar $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $q|\_c$ at one element vertex and two edge midpoint \
%[text] Scalar $SDG\_1$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_1 (K), \\forall K \\in T\_h; \\ \[q\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
D2TElem = MshEnt("D2T").msh;
SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", D2TElem, 0, [], d0_p, "EntIdx", 3), ...
    NdDoF("D2", D2TElem, 1, 1/2, d0_p, "EntIdx", [2, 3])]);
SDG_1S_BC = BC(p, "node", nan, "edge", DlEdge(1));
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);
%[text] Vector $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Vector $SDG\_1$ element space
%[text] - $U\_{h, 0} = \\{ v \\in L^2 (\\Omega, R^2) : v|\_K \\in P\_1 (K)^2, \\forall K \\in T\_h; \\ \[v\]\_e = 0, \\forall e \\in F^{pr, o}\_h; \\ v|\_\\Gamma = 0 \\}$ \
d0_ux = [0, nan; 0, nan]; d0_uy = [nan, 0; nan, 0];
SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", D2TElem, 1, [0, 1], d0_ux, "EntIdx", 1), NdDoF("D2", D2TElem, 1, [0, 1], d0_uy, "EntIdx", 1), ...
    NdDoF("D2", D2TElem, 0, [], d0_ux, "EntIdx", 3, "share", false), NdDoF("D2", D2TElem, 0, [], d0_uy, "EntIdx", 3, "share", false)]);
SDG_1V_BC = BC(Fcn("D2", "0"), "node", BdNode, "edge", BdEdge);
SDG_1V_FES = FES(msh, SDG_1V_FE, SDG_1V_BC);
%[text] Matrix $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P^{2\\times 2}\_1$
%[text] - Nodal DoF: $v n |\_c$ at edge vertex \
%[text] Matrix $SDG\_1$ element space
%[text] - $\\Sigma\_h = \\{ \\tau \\in L^2 (\\Omega, R^{2 \\times 2}) : \\tau|\_K \\in P\_1 (K)^{2 \\times 2}, \\forall K \\in T\_h; \\ \[\\tau \\cdot n\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
d0_sx = [0, 0; nan, nan; 0, 0; nan, nan]'; d0_sy = [nan, nan; 0, 0; nan, nan; 0, 0]';
SDG_1M_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 2]), ...
    [NdDoF("D2", D2TElem, 1, [0, 1], d0_sx, "coef", UNV, "EntIdx", [2, 3], "form", @(coef, fcn) sum(fcn * coef(1)), "orien", true), ...
    NdDoF("D2", D2TElem, 1, [0, 1], d0_sy, "coef", UNV, "EntIdx", [2, 3], "form", @(coef, fcn) sum(fcn * coef(1)), "orien", true), ...
    NdDoF("D2", D2TElem, 1, [0, 1], d0_sx, "coef", UNV, "EntIdx", 1, "form", @(coef, fcn) sum(fcn * coef(1)), "orien", true, "share", false), ...
    NdDoF("D2", D2TElem, 1, [0, 1], d0_sy, "coef", UNV, "EntIdx", 1, "form", @(coef, fcn) sum(fcn * coef(1)), "orien", true, "share", false)]);
SDG_1M_FES = FES(msh, SDG_1M_FE);
%%
%[text] Scalar $SDG\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P\_0$
%[text] - Nodal DoF: $q|\_c$ at one element vertex \
%[text] Scalar $SDG\_0$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_0 (K), \\forall K \\in T\_h; \\ \[q\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
D2TElem = MshEnt("D2T").msh;
SDG_0S_FE = FE("D2T", "1", NdDoF("D2", D2TElem, 0, [], d0_p, "EntIdx", 3));
SDG_0S_BC = BC(p, "node", DlNode(1));
SDG_0S_FES = FES(msh, SDG_0S_FE, SDG_0S_BC);
%[text] Vector $SDG\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_0$
%[text] - Nodal DoF: $v|\_c$ at midpoint of one edge \
%[text] Vector $SDG\_0$ element space
%[text] - $U\_{h, 0} = \\{ v \\in L^2 (\\Omega, R^2) : v|\_K \\in P\_0 (K)^2, \\forall K \\in T\_h; \\ \[v\]\_e = 0, \\forall e \\in F^{pr, o}\_h; \\ v|\_\\Gamma = 0 \\}$ \
d0_ux = [0, nan; 0, nan]; d0_uy = [nan, 0; nan, 0];
SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
    [NdDoF("D2", D2TElem, 1, 1/2, d0_ux, "EntIdx", 1), NdDoF("D2", D2TElem, 1, 1/2, d0_uy, "EntIdx", 1)]);
SDG_0V_BC = BC(Fcn("D2", "0"), "edge", BdEdge);
SDG_0V_FES = FES(msh, SDG_0V_FE, SDG_0V_BC);
%[text] Matrix $SDG\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P^{2\\times 2}\_0$
%[text] - Nodal DoF: $v n |\_c$ at midpoint of two edge \
%[text] Matrix $SDG\_0$ element space
%[text] - $\\Sigma\_h = \\{ \\tau \\in L^2 (\\Omega, R^{2 \\times 2}) : \\tau|\_K \\in P\_0 (K)^{2 \\times 2}, \\forall K \\in T\_h; \\ \[\\tau \\cdot n\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
d0_sx = [0, 0; nan, nan; 0, 0; nan, nan]'; d0_sy = [nan, nan; 0, 0; nan, nan; 0, 0]';
SDG_0M_FE = FE("D2T", FE.repFS("1", [2, 2]), ...
    [NdDoF("D2", D2TElem, 1, 1/2, d0_sx, "coef", UNV, "EntIdx", [2, 3], "form", @(coef, fcn) sum(fcn * coef(1)), "orien", true), ...
    NdDoF("D2", D2TElem, 1, 1/2, d0_sy, "coef", UNV, "EntIdx", [2, 3], "form", @(coef, fcn) sum(fcn * coef(1)), "orien", true)]);
SDG_0M_FES = FES(msh, SDG_0M_FE);
%%
%[text] SDG scheme: find $\\sigma\_h \\in \\Sigma\_h$, $u\_h \\in U\_{h, 0$, and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}\n  \\sum\_{K \\in T\_h} \\int\_K \\sigma\_h \\tau\_h \\, dx \\, dy + \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla \\cdot \\tau\_h \\, dx \\, dy - \\sum\_{e \\in F\_h^{pr}} \\int\_e u\_h \[ \\tau\_h n\] \\, ds = 0 &\\forall \\tau\_h \\in \\Sigma\_h \\\\\n  \\nu \\sum\_{K \\in T\_h} \\int\_K \\sigma\_h \\nabla v\_h \\, dx \\, dy - \\nu \\sum\_{e \\in F\_h^{dl}} \\int\_e \\sigma\_h n \[ v\_h \] \\, ds - \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla \\cdot v\_h \\, dx \\, dy + \\sum\_{e \\in F\_h^{dl}} \\int\_e p\_h \[ v\_h \\cdot n \] \\, ds = \\sum\_{K \\in T\_h} \\int\_K f v\_h \\, dx \\, dy &\\forall v\_h \\in U\_{h, 0} \\\\\n  - \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla q\_h \\, dx \\, dy + \\sum\_{e \\in F\_h^{pr}} \\int\_e u\_h \\cdot n \[ q\_h \] \\, ds = 0 &\\forall q\_h \\in P\_h\n\\end{cases}$
Sh = SDG_1M_FES; ord_Sh = 1;
Uh = SDG_1V_FES; ord_Uh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;

% Sh = SDG_0M_FES; ord_Sh = 1;
% Uh = SDG_0V_FES; ord_Uh = 1;
% Ph = SDG_0S_FES; ord_Ph = 1;

trls = [Sh, Uh, Ph]; tsts = [Sh, Uh, Ph];
is = 1; iu = 2; ip = 3; it = 1; iv = 2; iq = 3;

d0_t = d0_s; div_t = div_s; d0_v = d0_u; div_v = div_u; grad_v = grad_u; d0_q = d0_p; grad_q = grad_p;
Ast = DLF(msh, 2, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D2T", ord_Sh * 2));
But = [DLF(msh, 2, Fcn.cst(1), d0_u, div_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D2T", ord_Uh + ord_Sh - 1)), ...
    DLF.interface(msh, 1, -UNV, d0_u, d0_t, "EntIdx", PrEdge, "iTrl", iu, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D2L", ord_Uh + ord_Sh))];
Bsv = [DLF(msh, 2, Fcn.cst(v), d0_s, grad_v, "iTrl", is, "iTst", iv, "GInt", GInt("D2T", ord_Sh + ord_Uh - 1)), ...
    DLF.interface(msh, 1, -UNV * v, d0_s, d0_v, "EntIdx", DlEdge, "iTrl", is, "iTst", iv, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D2L", ord_Sh + ord_Uh))];
Bpv = [DLF(msh, 2, Fcn.cst(-1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1)), ...
    DLF.interface(msh, 1, UNV, d0_p, d0_v, "EntIdx", DlEdge, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D2L", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 2, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh + ord_Ph - 1)), ...
    DLF.interface(msh, 1, UNV, d0_u, d0_q, "EntIdx", PrEdge, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", ord_Uh + ord_Ph))];
Fv = SLF(msh, 2, f, d0_v, "iTst", iv, "GInt", GInt("D2T", 2 + ord_Uh));
%%
%[text] Solution.
[Stiff, Load] = assemble(msh, trls, tsts, [Ast, But, Bsv, Bpv, Buq], Fv);
[sh, uh, ph] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error.
%[text] $\\| v \\|\_{H^1, h}^2 = \\sum\_{K \\in T\_h} \\| \\nabla v \\|\_K^2 + \\sum\_{e \\in F\_h^{dl}} h^{-1}\_e \\| \[v\] \\|\_e^2$
L2Norm_Sh = Norm(msh, 2, d0_s, "GInt", GInt("D2T", (1 + ord_Sh) * 2));
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (1 + ord_Uh) * 2));
H1Norm_Uh = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", DlEdge, "coef", len \ 1, "form", @(coef, fcn, pow) coef .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (1 + ord_Ph) * 2));
fprintf('|s-sh|_L2: %e, |u-uh|_L2: %e, |u-uh|_H1: %e, |p-ph|_L2: %e\n', eNorm(msh, s, sh, L2Norm_Sh), eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, H1Norm_Uh), eNorm(msh, p, ph, L2Norm_Ph));
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
