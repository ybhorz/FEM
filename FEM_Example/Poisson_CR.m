%[text] Poisson equation
%[text] $\\begin{cases}\n  - \\Delta u = f & \\text{in } \\Omega \\\\\n  u = g\_D & \\text{on } \\Gamma\_D \\\\\n  \\frac{\\partial u}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
u = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y");
d0_u = [0; 0]; dx_u = [1; 0]; dy_u = [0; 1]; grad_u = cat(3, dx_u, dy_u); dxx_u = [2; 0]; dyy_u = [0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u);
g_D = u;
g_N = dot(MshEnt("D2L").UNV, dif(u, grad_u));
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured triangulation
%[text] Dirichlet boundary: lower and right boundary
%[text] Neumann boundary: upper and left boundary
msh = mshD2TS([0, 1, 0, 1], 32);
DirBd = [1, 2, -1, -2, -3];
NeuBd = [3, 4];
DirEdge = find(ismember(msh.edge.type, DirBd));
NeuEdge = find(ismember(msh.edge.type, NeuBd));
%%
%[text] CR element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at edge midpoint \
%[text] CR element space
%[text] - $U\_{h, g} = \\{ v \\in L^2 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h; \\ \[v\]\_{c(e)} = 0, \\forall e \\in E^o\_h; \\ v|\_{\\Gamma\_D} = g \\}$ \
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, d0_u));
CR_BC = BC(g_D, "edge", DirEdge);
CR_FES = FES(msh, CR_FE, CR_BC);
%%
%[text] CR scheme: find $u\_h \\in U\_{h, g}$ such that
%[text] $\\int\_{\\Omega} \\nabla u\_h \\cdot \\nabla v\_h \\, dx \\, dy = \\int\_{\\Omega} f v\_h \\, dx \\, dy + \\int\_{\\Gamma\_N} g\_N v\_h \\, ds , \\quad \\forall v\_h \\in U\_{h, 0}$
Uh = CR_FES; ord = 1;

d0_v = d0_u; grad_v = grad_u;
Auv = DLF(msh, 2, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D2T", (ord - 1) * 2));
Fv = [SLF(msh, 2, f, d0_v, "GInt", GInt("D2T", 1 + ord)), ...
    SLF(msh, 1, g_N, d0_v, "EntIdx", NeuEdge, "GInt", GInt("D2L", 1 + ord))];
%%
%[text] Solution
[Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv);
uh = FEF(Uh, Stiff \ Load);
%%
%[text] Error
%[text] $\\| v \\|^2\_{L^2} = \\int\_{\\Omega} |v|^2 \\, dx \\, dy \\\\\n\\| v \\|^2\_{H^1, h} = \\sum\_{K \\in T\_h} \\int\_K |\\nabla v|^2 \\, dx \\, dy$
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
H1Norm = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", (ord - 1) * 2))];
fprintf('|u-uh|_L2: %e, |u-uh|_H1: %e\n', eNorm(msh, u, uh, L2Norm), eNorm(msh, u, uh, H1Norm));
%%
%[text] Plot
plotFEF(uh);
plotFcn(msh, u);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
