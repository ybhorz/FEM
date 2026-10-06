%[text] Poisson equation (3D)
%[text] $\\begin{cases}\n  - \\Delta u = f & \\text{in } \\Omega \\\\\n  u = g\_D & \\text{on } \\Gamma\_D \\\\\n  \\frac{\\partial u}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
u = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z)+x*y*z");
d0_u = [0; 0; 0]; dx_u = [1; 0; 0]; dy_u = [0; 1; 0]; dz_u = [0; 0; 1]; grad_u = cat(3, dx_u, dy_u, dz_u);
dxx_u = [2; 0; 0]; dyy_u = [0; 2; 0]; dzz_u = [0; 0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u) - dif(u, dzz_u);
g_D = u;
UNV = MshEnt("D3F").UNV;
g_N = dot(UNV, dif(u, grad_u));
%%
%[text] Domain: cube $\[0, 1\] \\times \[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured tetrahedral mesh (each cube split into 6 tetrahedra)
%[text] Dirichlet boundary: faces $x = 0$, $y = 0$, $z = 0$
%[text] Neumann boundary: faces $x = 1$, $y = 1$, $z = 1$
msh = mshD3TS([0, 1, 0, 1, 0, 1], 4);
DirBd = [1, 3, 5];
NeuBd = [2, 4, 6];
DirNode = msh.bdEnt(0, DirBd);
NeuFace = msh.bdEnt(2, NeuBd);
IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, [DirBd, NeuBd]));
%%
%[text] $P\_1$ element
%[text] - Element: tetrahedron
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Discontinuous $P\_1$ element space
%[text] - $U\_{h, g} = \\{ v \\in L^2 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h; \\ v|\_{\\Gamma\_D} = g \\}$ \
DG1_FE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_u, "share", false), "map", "affine");
DG1_BC = BC(g_D, "node", DirNode);
DG1_FES = FES(msh, DG1_FE, DG1_BC);
%%
%[text] DG scheme (SIPG): find $u\_h \\in U\_{h, g}$ such that
%[text] $\\sum\_{K \\in T\_h} \\int\_{K} \\nabla u\_h \\cdot \\nabla v\_h \\, dx - \\sum\_{F \\in F^o\_h} \\int\_{F} \\{ \\frac{\\partial u\_h}{\\partial n} \\} \[ v\_h \] \\, ds - \\sum\_{F \\in F^o\_h} \\int\_{F} \\{ \\frac{\\partial v\_h}{\\partial n} \\} \[ u\_h \] \\, ds + \\sum\_{F \\in F^o\_h} \\int\_{F} \\frac{\\gamma}{h\_F} \[ u\_h \] \[ v\_h \] \\, ds =  \\sum\_{K \\in T\_h} \\int\_{K} f v\_h \\, dx + \\sum\_{F \\in \\Gamma^N\_h} \\int\_F g\_N v\_h \\, ds, \\quad \\forall v\_h \\in U\_{h, 0}$
%[text] Face size: $h\_F = |F|^{1/2}$.
Uh = DG1_FES; ord = 1; gamma = 20;
hF = MshEnt("D3F").area .^ (1/2);
grad_v = grad_u; d0_v = d0_u;
Auv = [DLF(msh, 3, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D3T", (ord - 1) * 2)), ...
    DLF.interface(msh, 2, - UNV, grad_u, d0_v, "EntIdx", IntFace, "trlOpr", "aver", "tstOpr", "jump", "GInt", GInt("D3F", ord - 1 + ord)), ...
    DLF.interface(msh, 2, - UNV, d0_u, grad_v, "EntIdx", IntFace, "trlOpr", "jump", "tstOpr", "aver", "GInt", GInt("D3F", ord + ord - 1)), ...
    DLF.interface(msh, 2, hF \ gamma, d0_u, d0_v, "EntIdx", IntFace, "trlOpr", "jump", "tstOpr", "jump", "GInt", GInt("D3F", ord * 2))];
Fv = [SLF(msh, 3, f, d0_v, "GInt", GInt("D3T", 1 + ord)), ...
    SLF(msh, 2, g_N, d0_v, "EntIdx", NeuFace, "GInt", GInt("D3F", 1 + ord))];
%%
%[text] Solution.
[Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv);
uh = FEF(Uh, Stiff \ Load);
%%
%[text] Error.
%[text] $|| v ||\_{H^1, h}^2 =  \\sum\_{K \\in T\_h} \\int\_{K}  | \\nabla v |^2 \\, dx + \\sum\_{F \\in F^o\_h} \\frac{1}{h\_F} \\int\_{F} \[ v \]^2 \\, ds.$
L2Norm = Norm(msh, 3, d0_u, "GInt", GInt("D3T", ord * 2));
H1Norm = [Norm(msh, 3, grad_u, "GInt", GInt("D3T", ord * 2)), ...
    Norm(msh, 2, d0_u, "EntIdx", IntFace, "coef", hF \ 1, "fcnOpr", "jump", "GInt", GInt("D3F", ord * 2))];
fprintf('|u-uh|_L2: %e, |u-uh|_H1: %e\n', eNorm(msh, u, uh, L2Norm), eNorm(msh, u, uh, H1Norm));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
