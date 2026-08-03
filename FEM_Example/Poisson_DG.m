%[text] Poisson equation
%[text] $\\begin{cases}\n  - \\Delta u = f & \\text{in } \\Omega \\\\\n  u = g\_D & \\text{on} \\Gamma\_D \\\\\n  \\frac{\\partial u}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
u = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y");
d0_u = [0; 0]; dx_u = [1; 0]; dy_u = [0; 1]; grad_u = cat(3, dx_u, dy_u); dxx_u = [2; 0]; dyy_u = [0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u);
g_D = u;
UNV = MshEnt("D2L").UNV;
g_N = dot(UNV, dif(u, grad_u));
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured triangulation
%[text] Dirichlet boundary: lower and right boundary
%[text] Neumann boundary: upper and left boundary
msh = mshD2TS([0, 1, 0, 1], 16);
DirBd = [1, 2, -1, -2, -3];
NeuBd = [3, 4];
DirNode = find(ismember(msh.node.type, DirBd));
DirEdge = find(ismember(msh.edge.type, DirBd));
NeuEdge = find(ismember(msh.edge.type, NeuBd));
IntEdge = find(ismember(msh.edge.type, 0));
%%
%[text] $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Discontinuous $P\_1$ element space
%[text] - $U\_{h, g} = \\{ v \\in L^2 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h; \\ v|\_{\\Gamma\_D} = g \\}$ \
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], d0_u, "share", false));
DG1_BC = BC(g_D, "node", DirNode);
DG1_FES = FES(msh, DG1_FE, DG1_BC);
%%
%[text] DG scheme (SIPG): find $u\_h \\in U\_{h, g}$ such that
%[text] $\\sum\_{K \\in T\_h} \\int\_{K} \\nabla u\_h \\cdot \\nabla v\_h \\, dx \\, dy - \\sum\_{e \\in E^o\_h} \\int\_{e} \\{ \\frac{\\partial u\_h}{\\partial n} \\} \[ v\_h \] \\, ds - \\sum\_{e \\in E^o\_h} \\int\_{e} \\{ \\frac{\\partial v\_h}{\\partial n} \\} \[ u\_h \] \\, ds + \\sum\_{e \\in E^o\_h} \\int\_{e} \\frac{\\gamma}{h\_e} \[ u\_h \] \[ v\_h \] \\, ds =  \\sum\_{K \\in T\_h} \\int\_{K} f v\_h \\, dx \\, dy + \\sum\_{e \\in \\Gamma^N\_h} \\int\_e g\_N v\_h \\, ds, \\quad \\forall v\_h \\in U\_{h, 0}$
Uh = DG1_FES; ord = 1; gamma = 10;
grad_v = grad_u; d0_v = d0_u;
Auv = [DLF(msh, 2, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D2T", (ord - 1) * 2)), ...
    DLF.interface(msh, 1, - UNV, grad_u, d0_v, "EntIdx", IntEdge, "trlOpr", "aver", "tstOpr", "jump", "GInt", GInt("D2L", ord - 1 + ord)), ...
    DLF.interface(msh, 1, - UNV, d0_u, grad_v, "EntIdx", IntEdge, "trlOpr", "jump", "tstOpr", "aver", "GInt", GInt("D2L", ord + ord - 1)), ...
    DLF.interface(msh, 1, MshEnt("D2L").len \ gamma, d0_u, d0_v, "EntIdx", IntEdge, "trlOpr", "jump", "tstOpr", "jump", "GInt", GInt("D2L", ord * 2))];
Fv = [SLF(msh, 2, f, d0_v, "GInt", GInt("D2T", 1 + ord)), ...
    SLF(msh, 1, g_N, d0_v, "EntIdx", NeuEdge, "GInt", GInt("D2L", 1 + ord))];
%%
%[text] Solution.
[Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv);
uh = FEF(Uh, Stiff \ Load);
%%
%[text] Error.
%[text] $|| v ||\_{H^1, h}^2 =  \\sum\_{K \\in T\_h} \\int\_{K}  | \\nabla v |^2 \\, dx \\, dy + \\sum\_{e \\in E^o\_h} \\frac{1}{h\_e} \\int\_{e} \[ e \]^2 \\, ds.$
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
H1Norm = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", (ord - 1) * 2)), ...
    Norm(msh, 1, d0_u, "coef", MshEnt("D2L").len \ 1, "fcnOpr", "jump", "GInt", GInt("D2L", ord * 2))];
fprintf('|u-uh|_L2: %e, |u-uh|_H1: %e\n', eNorm(msh, u, uh, L2Norm), eNorm(msh, u, uh, H1Norm));
%%
%[text] Plot.
plotFEF(uh);
plotFcn(msh, u);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
