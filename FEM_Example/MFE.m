%[text] Poisson equation:
%[text] $\\begin{cases}  - \\Delta p = f & \\text{in } \\Omega, \\\\  p = g\_D & \\text{on}\n\\Gamma\_D, \\\\  \\frac{\\partial p}{\\partial n} = g\_N & \\text{on } \\Gamma\_N.  \\end{cases}$
%[text] Mixed formulation:
%[text] $\\begin{cases}  u + \\nabla p = 0 & \\text{in } \\Omega, \\\\   \\nabla \\cdot u\n= f & \\text{in } \\Omega, \\\\  p = g\_D & \\text{on } \\Gamma\_D, \\\\ u \\cdot n = -\ng\_N & \\text{on } \\Gamma\_N.  \\end{cases}$
p = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y");
d0_p = [0; 0]; dx_p = [1; 0]; dy_p = [0; 1]; grad_p = cat(3, dx_p, dy_p);
u = - dif(p, grad_p);
d0_u = [0, 0; 0, 0]; dxdy_u = [1, 0; 0, 1];
f = sum(dif(u, dxdy_u));
g_D = p;
UNV = MshEnt("D2L").UNV;
g_N = - dot(u, UNV);
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$.
%[text] Dirichlet boundary: lower and right boundary.
%[text] Neumann boundary: upper and left boundary.
%[text] Mesh: structured triangulation.
msh = mshD2TS([0, 1, 0, 1], 16);
DirBd = [1, 2];
NeuBd = [3, 4];
DirEdge = find(ismember(msh.edge.type, DirBd));
NeuEdge = find(ismember(msh.edge.type, NeuBd));
%%
%[text] RT0 element.
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'",...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", UNV, "orien", true));
RT0_BC = BC(- g_N, "edge", NeuEdge);
RT0_FES = FES(msh, RT0_FE, RT0_BC);
%%
%[text] P0 element.
P0_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], d0_p));
P0_FES = FES(msh, P0_FE);
%%
%[text] MFE scheme: find $u\_h \\in U\_h$ and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}   \\int\_\\Omega u\_h \\cdot v\_h \\, dx - \\int\_\\Omega p\_h \\nabla\n\\cdot v\_h \\, dx = - \\int\_{\\Gamma\_D} g\_D v\_h \\cdot n \\, ds, & \\forall v\_h \\in\nV\_h. \\\\  \\int\_\\Omega \\nabla \\cdot u\_h q\_h \\, dx = \\int\_\\Omega f q\_h \\, dx, &\n\\forall q\_h \\in Q\_h, \\end{cases}$
Uh = RT0_FES; ord_Uh = 1;
Ph = P0_FES; ord_Ph = 0;
trls = [Uh, Ph]; tsts = [Uh, Ph];
iu = 1; ip = 2; iv = 1; iq = 2;
d0_v = d0_u; dxdy_v = dxdy_u; d0_q = d0_p;
Auv = DLF(msh, 2, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D2T", ord_Uh * 2));
Bpv = DLF(msh, 2, Fcn.cst(-1), d0_p, dxdy_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1));
Buq = DLF(msh, 2, Fcn.cst(1), dxdy_u, d0_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh - 1 + ord_Ph));
Gv = SLF(msh, 1, -g_D .* UNV, d0_v, "iTst", iv, "EntIdx", DirEdge,"GInt", GInt("D2L", 1 + ord_Uh));
Fq = SLF(msh, 2, f, d0_q, "iTst", iq, "GInt", GInt("D2T", 1 + ord_Ph));
%%
%[text] Solution.
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Buq], [Gv, Fq]);
[uh, ph] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error.
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (ord_Uh + 1) * 2));
divNorm_Uh = Norm(msh, 2, dxdy_u, "form",  @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow, "GInt", GInt("D2T", ord_Uh * 2));
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (ord_Ph + 1) * 2));
fprintf('L2 norm of uh error: %e, div norm of uh error: %e, L2 norm of ph error: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, divNorm_Uh), eNorm(msh, p, ph, L2Norm_Ph))
%%
%[text] Plot.
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
