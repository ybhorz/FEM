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
%[text] Scalar $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $q|\_c$ at element vertex \
%[text] Scalar $P\_1$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h \\}$ \
SP1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], d0_p, "share", false));
SP1_FES = FES(msh, SP1_FE);
%%
%[text] Vector $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1^2$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Vector $P\_1$ element space
%[text] - $U\_h = \\{ v \\in L^2 (\\Omega; R^2) : v|\_{K} \\in P\_1 (K)^2, \\forall K \\in T\_h \\}$ \
d0_ux = [0, nan; 0, nan]; d0_uy = [nan, 0; nan, 0];
VP1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), [NdDoF("D2", MshEnt("D2T").msh, 0, [], d0_ux, "share", false), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], d0_uy, "share", false)]);
VP1_FES = FES(msh, VP1_FE);
%%
%[text] Trace $P\_1$ element
%[text] - Element: reference edge \[0, 1\]
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $\\hat{\\mu} |\_c$ at edge vertex \
%[text] Trace $P\_1$ element space
%[text] - $\\Lambda\_{h, g} = \\{ \\mu \\in L^2 (E\_h) : \\mu |\_e \\in P\_1 (e), \\forall e \\in E\_h; \\ \\mu |\_e = g\_D, \\forall e \\in \\Gamma^D\_h \\}$ \
d0_l = 0;
TP1_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], d0_l, "share", false));
TP1_BC = BC(g_D, "edge", DirEdge);
TP1_FES = FES(msh, TP1_FE, TP1_BC);
%%
%[text] HDG scheme: find $u\_h \\in U\_h$, $p\_h \\in P\_h$ and $\\Lambda\_{h, g}$ such that
%[text] $\\begin{cases}\n  \\int\_\\Omega u\_h \\cdot v\_h \\, dx \\, dy + \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla \\cdot v\_h \\, dx \\, dy - \\sum\_{e \\in E\_h} \\int\_e \\lambda\_h \[v\_h \\cdot n\] \\, ds = 0, & \\forall v\_h \\in U\_h \\\\\n  \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla q\_h \\, dx \\, dy - \\sum\_{e \\in E\_h} \\int\_e (u\_h \\cdot n + \\tau (p\_h - \\lambda\_h)) \[q\_h\] \\, ds  = \\int\_\\Omega f q\_h \\, dx \\, dy, & \\forall q\_h \\in P\_h \\\\\n  \\sum\_{e \\in E\_h} \\int\_e \[u\_h \\cdot n + \\tau (p\_h - \\lambda\_h)\] \\mu\_h \\,ds = \\sum\_{e \\in \\Gamma^N\_h} \\int\_e g\_N \\mu\_h \\, ds, & \\forall \\mu\_h \\in \\Lambda\_{h, 0}\n\\end{cases}$
tau = 1;
Uh = VP1_FES; ord_Uh = 1;
Ph = SP1_FES; ord_Ph = 1;
Lh = TP1_FES; ord_Lh = 1;

trls = [Uh, Ph, Lh]; tsts = [Uh, Ph, Lh];
iu = 1; ip = 2; il = 3; iv = 1; iq = 2; im = 3;
d0_v = d0_u; div_v = div_u; d0_q = d0_p; grad_q = grad_p; d0_m = d0_l;

Auv = DLF(msh, 2, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D2T", ord_Uh * 2));
Bpv = DLF(msh, 2, Fcn.cst(1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1));
Blv = DLF.interface(msh, 1, -UNV, d0_l, d0_v, "iTrl", il, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D2L",ord_Lh + ord_Uh));

Buq = [DLF(msh, 2, Fcn.cst(1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh + ord_Ph - 1)), ...
    DLF.interface(msh, 1, -UNV, d0_u, d0_q, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L",ord_Uh + ord_Ph))];
Apq = DLF.interface(msh, 1, Fcn.cst(-tau), d0_p, d0_q, "iTrl", ip, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", ord_Ph *2));
Blq = DLF.interface(msh, 1, Fcn.cst(tau), d0_l, d0_q, "iTrl", il, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", ord_Lh + ord_Ph));
Fq = SLF(msh, 2, f, d0_q, "iTst", iq, "GInt", GInt("D2T", 1 + ord_Ph));

Bum = DLF.interface(msh, 1, UNV, d0_u, d0_m, "iTrl", iu, "iTst", im, "trlOpr", "jump", "GInt", GInt("D2L", ord_Uh + ord_Lh));
Bpm = DLF.interface(msh, 1, Fcn.cst(tau), d0_p, d0_m, "iTrl", ip, "iTst", im, "trlOpr", "jump", "GInt", GInt("D2L", ord_Ph + ord_Lh));
Alm = DLF(msh, 1, -Fcn.cst(tau), d0_l, d0_m, "iTrl", il, "iTst", im, "EntIdx", NeuEdge, "GInt", GInt("D2L", ord_Lh * 2));
Gm = SLF(msh, 1, g_N, d0_m, "iTst", im,"EntIdx", NeuEdge, "GInt", GInt("D2L", 1 + ord_Lh));
%%
%[text] Solution
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Blv, Buq, Apq, Blq, Bum, Bpm, Alm], [Fq, Gm]);
[uh, ph, lh] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (ord_Uh + 1) * 2));
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (ord_Ph + 1) * 2));
L2Norm_Lh = Norm(msh, 1, d0_l, "GInt", GInt("D2L", (ord_Lh + 1) * 2));
fprintf('|u-uh|_L2: %e, |p-ph|_L2: %e, |l-l_h|_L2: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, p, ph, L2Norm_Ph), eNorm(msh, p, lh, L2Norm_Lh));
%%
%[text] Plot
plotFEF(uh);
plotFcn(msh, u);
%%
plotFEF(ph);
plotFcn(msh, p);
%%
plotFEF(lh);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
