%[text] Poisson equation
%[text] $\\begin{cases}\n  - \\Delta p = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on }\\Gamma\_D \\\\\n  \\frac{\\partial p}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  u - \\nabla p = 0 & \\text{in } \\Omega \\\\ \n  -\\nabla \\cdot u = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on } \\Gamma\_D \\\\\n  u \\cdot n = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
p = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y");
d0_p = [0; 0]; dx_p = [1; 0]; dy_p = [0; 1]; grad_p = cat(3, dx_p, dy_p);
u = dif(p, grad_p);
d0_u = [0, 0; 0, 0]; div_u = [1, 0; 0, 1];
f = -sum(dif(u, div_u));
g_D = p;
UNV = MshEnt("D2L").UNV;
g_N = dot(u, UNV);
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured triangulation + Alfeld splitting
%[text] Dirichlet boundary: lower and right boundary
%[text] Neumann boundary: upper and left boundary
msh = mshSplit(mshD2TS([0, 1, 0, 1], 16));
DirBd = [1, 2, -1, -2, -3];
NeuBd = [3, 4];
DirNode = find(ismember(msh.node.type, DirBd));
DirEdge = find(ismember(msh.edge.type, DirBd));
NeuEdge = find(ismember(msh.edge.type, NeuBd));
PrEdge = find(ismember(msh.edge.type, [0, 1, 2, 3, 4]));
PrOEdge = find(ismember(msh.edge.type, 0));
DlEdge = find(ismember(msh.edge.type, 1i));
%%
%[text] Scalar $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Moment DoF: $\\int\_e q \\, ds$ and $\\int\_e q (s-1/2) \\, ds$ on one edge, and $\\int\_K q \\, dx \\, dy$ on element
%[text] - or Nodal DoF: $q|\_c$ at element vertex \
%[text] Scalar $SDG\_1$ element space
%[text] - $P\_{h, g} = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_1 (K), \\forall K \\in T\_h; \\ \[q\]\_e = 0, \\forall e \\in F^{pr, o}\_h; \\ q|\_{\\Gamma\_D} = g \\}$ \
SDG_1S_FE = FE("D2T", "[1,x,y]", [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), d0_p, "EntIdx", 1, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), d0_p, "EntIdx", 1, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 2, Fcn.cst(1), d0_p, "GInt", GInt("D2T", 2))]);
% SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], d0_p, "EntIdx", 1), ...
    % NdDoF("D2", MshEnt("D2T").msh, 0, [], d0_p, "EntIdx", 3, "share", false)]);
SDG_1S_BC = BC(g_D, "node", DirNode, "edge", DirEdge);
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);
%[text] Vector $SDG\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_1$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ and $\\int\_e v \\cdot n (s - 1/2) \\, ds$ on two edge, and $\\int\_K v \\, dx \\, dy$ on element
%[text] - or Nodal DoF: $v \\cdot n |\_c$ at edge vertex \
%[text] Vector $SDG\_1$ element space
%[text] - $U\_h = \\{ v \\in L^2 (\\Omega; R^2) : v|\_K \\in P\_1 (K)^2, \\forall K \\in T\_h; \\ \[v \\cdot n\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
SDG_1V_FE  = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s - 1/2"), d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 2, [Fcn.cst([1;0]), Fcn.cst([0;1])], d0_u, "form", @(coef, fcn, tst) sum(coef .* fcn .* tst), "GInt", GInt("D2T", 2))]);
% SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
%     [NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", 1, "orien", true, "share", false)]);
SDG_1V_FES = FES(msh, SDG_1V_FE);
%%
%[text] Scalar $SDG\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P\_0$
%[text] - Moment DoF: $\\int\_e q \\, ds$ on one edge
%[text] - or Nodal DoF: $q|\_c$ at midpoint  of  one edge \
%[text] Scalar $SDG\_0$ element space
%[text] - $P\_{h, g} = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_0 (K), \\forall K \\in T\_h; \\ \[q\]\_e = 0, \\forall e \\in F^{pr, o}\_h; \\ q|\_{\\Gamma\_D} = g \\}$ \
SDG_0S_FE = FE("D2T", "1", MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), d0_p, "EntIdx", 1, "GInt", GInt("D2L", 4)));
% SDG_0S_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, d0_p, "EntIdx", 1));
SDG_0S_BC = BC(g_D, "edge", DirEdge);
SDG_0S_FES = FES(msh, SDG_0S_FE, SDG_0S_BC);
%[text] Vector $SDG\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_0$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ on two edge
%[text] - or Nodal DoF: $v \\cdot n |\_c$ at midpoint of two edge \
%[text] Vector $SDG\_0$ element space
%[text] - $U\_h = \\{ v \\in L^2 (\\Omega; R^2) : v|\_K \\in P\_0 (K)^2, \\forall K \\in T\_h; \\ \[v \\cdot n\]\_e = 0, \\forall e \\in F^{dl}\_h \\}$ \
SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true, "GInt", GInt("D2L", 4)));
% SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true));
SDG_0V_FES = FES(msh, SDG_0V_FE);
%%
%[text] SDG scheme: find $u\_h \\in U\_h$ and $p\_h \\in P\_{h, g}$ such that
%[text] $\\begin{cases}\n  \\sum\_{K \\in T\_h} \\int\_K u\_h \\cdot v\_h \\, dx \\, dy + \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla\\cdot v\_h \\, dx \\, dy - \\sum\_{e \\in F\_h^{pr}} \\int\_e p\_h \[v\_h \\cdot n\] \\, ds  = 0, & \\forall v\_h \\in U\_h \\\\\n  \\sum\_{K \\in T\_h} \\int\_K  u\_h \\nabla q\_h \\, dx \\, dy - \\sum\_{e \\in F\_h^{dl}} \\int\_e u\_h \\cdot n \[q\_h\] \\, ds = \\sum\_{K \\in T\_h} \\int\_K f q\_h \\, dx \\, dy + \\sum\_{e \\in \\Gamma\_h^N} \\int\_e g\_N q\_h \\, ds, &\\forall q\_h \\in P\_{h, 0}\n\\end{cases}$
Uh = SDG_1V_FES; ord_Uh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;

% Uh = SDG_0V_FES; ord_Uh = 1;
% Ph = SDG_0S_FES; ord_Ph = 1;

trls = [Uh, Ph]; tsts = [Uh, Ph];
iu = 1; ip = 2; iv = 1; iq = 2;
d0_v = d0_u; div_v = div_u; d0_q = d0_p; grad_q = grad_p;
Auv = DLF(msh, 2, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D2T", ord_Uh * 2));
Bpv = [DLF(msh, 2, Fcn.cst(1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1)), ...
    DLF.interface(msh, 1, -UNV, d0_p, d0_v, "EntIdx", PrEdge, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D2L", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 2, Fcn.cst(1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh + ord_Ph - 1)), ...
    DLF.interface(msh, 1, -UNV, d0_u, d0_q, "EntIdx", DlEdge, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", ord_Uh + ord_Ph))];
Fq = [SLF(msh, 2, f, d0_q, "iTst", iq, "GInt", GInt("D2T", 1 + ord_Ph)), ...
    SLF(msh, 1, g_N, d0_q, "iTst", iq, "EntIdx", NeuEdge, "GInt", GInt("D2L", 1  + ord_Ph))];
%%
%[text] Solution
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Buq], Fq);
[uh, ph] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error
%[text] $|| v ||^2\_{div, h} = \\sum\_{K \\in T\_h} || \\nabla \\cdot v ||^2\_K + \\sum\_{e \\in F\_h^{pr, o}} h^{-1}\_e || \[v \\cdot n\]||^2\_e \\\\\n|| q ||^2\_{H^1, h} = \\sum\_{K \\in T\_h} ||\\nabla q ||^2\_K + \\sum\_{e \\in F^{dl}\_h} h^{-1}\_e || \[q\] ||^2\_e$
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (1 + ord_Uh) * 2));
divNorm_Uh = [Norm(msh, 2, div_u, "form", @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow, "GInt", GInt("D2T", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", PrOEdge, "coef", [MshEnt("D2L").len \ 1, UNV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(coef(2), fcn)) .^ pow, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (1 + ord_Ph) * 2));
H1Norm_Ph = [Norm(msh, 2, grad_p, "GInt", GInt("D2T", ord_Ph * 2)), ...
    Norm(msh, 1, d0_p, "EntIdx", DlEdge, "coef", MshEnt("D2L").len \ 1, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Ph * 2))];
fprintf('|u-uh|_L2: %e, |u-uh|_div: %e, |p-ph|_L2: %e, |p-ph|_H1: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, divNorm_Uh), eNorm(msh, p, ph, L2Norm_Ph), eNorm(msh, p, ph, H1Norm_Ph));

Pu = Uh.proj(u);
Pp = Ph.proj(p);
fprintf('|Pu-uh|_L2: %e, |Pu-uh|_div: %e, |Pp-ph|_L2: %e, |Pp-ph|_H1: %e\n', eNorm(msh, Fcn.cst(0), subt(Pu, uh), L2Norm_Uh), eNorm(msh, Fcn.cst(0), subt(Pu, uh), divNorm_Uh), eNorm(msh, Fcn.cst(0), subt(Pp, ph), L2Norm_Ph), eNorm(msh, Fcn.cst(0), subt(Pp, ph), H1Norm_Ph));
%%
%[text] Plot
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
