%[text] Poisson equation
%[text] $\\begin{cases}\n  - \\Delta p = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on } \\Gamma\_D \\\\\n  \\frac{\\partial p}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  u - \\nabla p = 0 & \\text{in } \\Omega \\\\\n  -\\nabla \\cdot u = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on } \\Gamma\_D \\\\\n  u \\cdot n = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
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
%[text] Mesh: structured triangulation
%[text] Dirichlet boundary: lower and right boundary
%[text] Neumann boundary: upper and left boundary
msh = mshD2TS([0, 1, 0, 1], 16);
DirBd = [1, 2];
NeuBd = [3, 4];
DirEdge = find(ismember(msh.edge.type, DirBd));
NeuEdge = find(ismember(msh.edge.type, NeuBd));
%%
%[text] $RT\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_0 + x P\_0$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ on edge
%[text] - or Nodal DoF: $v \\cdot n |\_c$ at edge midpoint \
%[text] $RT\_0$ element space
%[text] - $U\_{h, g}  = \\{ v \\in H(div; \\Omega) : v|\_K \\in RT\_0 (K), \\forall K \\in T\_h; \\ v \\cdot n = g \\}$ \
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'",...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), d0_u, "coef", UNV, "orien", true, "GInt", GInt("D2L", 4)));
% RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'",...
%     NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, d0_u, "coef", UNV, "orien", true));
RT0_BC = BC(g_N, "edge", NeuEdge);
RT0_FES = FES(msh, RT0_FE, RT0_BC);
%%
%[text] $BDM\_1$ element
%[text] - Domain: triangle
%[text] - Function space: $P^2\_1$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ and $\\int\_e (v \\cdot n) (s - 1/2) \\, ds$ on edge
%[text] - or Nodal DoF: $v \\cdot n |\_c$ at edge vertex \
%[text] $BDM\_1$ element space
%[text] - $U\_{h, g}  = \\{ v \\in H(div; \\Omega) : v|\_K \\in BDM\_1 (K), \\forall K \\in T\_h; \\ v \\cdot n |\_{\\Gamma\_N} = g \\}$ \
BDM1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), d0_u, "coef", UNV, "orien", true, "GInt", GInt("D2L", 4)), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), d0_u, "coef", UNV, "orien", true, "GInt", GInt("D2L", 4))]);
% BDM1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
%     NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], d0_u, "coef", UNV, "orien", true));
BDM1_BC = BC(g_N, "edge", NeuEdge);
BDM1_FES = FES(msh, BDM1_FE, BDM1_BC);
%%
%[text] $P\_0$ element
%[text] - Domain: triangle
%[text] - Function space: $P\_0$
%[text] - Moment DoF: $\\int\_K q \\, dx \\, dy$ on element
%[text] - or Nodal DoF: $q|\_c$ at centroid \
%[text] $P\_0$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_0 (K), \\forall K \\in T\_h \\}$ \
P0_FE = FE("D2T", "1", MoDoF("D2", MshEnt("D2T").msh, 2, Fcn.cst(1), d0_p, "GInt", GInt("D2T", 2)));
% P0_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], d0_p));
P0_FES = FES(msh, P0_FE);
%%
%[text] MFE scheme: find $u\_h \\in U\_{h, g}$ and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}\n  \\int\_\\Omega u\_h \\cdot v\_h \\, dx \\, dy + \\int\_\\Omega p\_h \\nabla \\cdot v\_h \\, dx \\, dy = \\int\_{\\Gamma\_D} g\_D v\_h \\cdot n \\, ds, & \\forall v\_h \\in U\_{h, 0} \\\\\n  -\\int\_\\Omega \\nabla \\cdot u\_h q\_h \\, dx \\, dy = \\int\_\\Omega f q\_h \\, dx \\, dy, & \\forall q\_h \\in P\_h\n\\end{cases}$
Uh = RT0_FES; ord_Uh = 1;
Ph = P0_FES; ord_Ph = 0;

% Uh = BDM1_FES; ord_Uh = 1;
% Ph = P0_FES; ord_Ph = 0;

trls = [Uh, Ph]; tsts = [Uh, Ph];
iu = 1; ip = 2; iv = 1; iq = 2;
d0_v = d0_u; div_v = div_u; d0_q = d0_p;
Auv = DLF(msh, 2, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D2T", ord_Uh * 2));
Bpv = DLF(msh, 2, Fcn.cst(1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1));
Buq = DLF(msh, 2, Fcn.cst(-1), div_u, d0_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh - 1 + ord_Ph));
Gv = SLF(msh, 1, g_D .* UNV, d0_v, "iTst", iv, "EntIdx", DirEdge,"GInt", GInt("D2L", 1 + ord_Uh));
Fq = SLF(msh, 2, f, d0_q, "iTst", iq, "GInt", GInt("D2T", 1 + ord_Ph));
%%
%[text] Solution
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Buq], [Gv, Fq]);
[uh, ph] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (ord_Uh + 1) * 2));
divNorm_Uh = Norm(msh, 2, div_u, "form",  @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow, "GInt", GInt("D2T", ord_Uh * 2));
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (ord_Ph + 1) * 2));
fprintf('|u-uh|_L2: %e, |u-uh|_div: %e, |p-ph|_L2: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, divNorm_Uh), eNorm(msh, p, ph, L2Norm_Ph))

Pu = Uh.proj(u);
Pp = Ph.proj(p);
fprintf('|Pu-uh|_L2: %e, |Pu-uh|_div: %e, |Pp-ph|_L2: %e\n', eNorm(msh, Fcn.cst(0), subt(Pu, uh), L2Norm_Uh), eNorm(msh, Fcn.cst(0), subt(Pu, uh), divNorm_Uh), eNorm(msh, Fcn.cst(0), subt(Pp, ph), L2Norm_Ph))
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
