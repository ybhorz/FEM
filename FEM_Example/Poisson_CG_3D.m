%[text] Poisson equation (3D)
%[text] $\\begin{cases}\n  - \\Delta u = f & \\text{in } \\Omega \\\\\n  u = g\_D & \\text{on } \\Gamma\_D \\\\\n  \\frac{\\partial u}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
u = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z)+x*y*z");
d0_u = [0; 0; 0]; dx_u = [1; 0; 0]; dy_u = [0; 1; 0]; dz_u = [0; 0; 1]; grad_u = cat(3, dx_u, dy_u, dz_u);
dxx_u = [2; 0; 0]; dyy_u = [0; 2; 0]; dzz_u = [0; 0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u) - dif(u, dzz_u);
g_D = u;
g_N = dot(MshEnt("D3F").UNV, dif(u, grad_u));
%%
%[text] Domain: cube $\[0, 1\] \\times \[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured tetrahedral mesh (each cube split into 6 tetrahedra)
%[text] Dirichlet boundary: faces $x = 0$, $y = 0$, $z = 0$
%[text] Neumann boundary: faces $x = 1$, $y = 1$, $z = 1$
msh = mshD3TS([0, 1, 0, 1, 0, 1], 4);
DirBd = [1, 3, 5];
NeuBd = [2, 4, 6];
DirNode = msh.bdEnt(0, DirBd);
DirEdge = msh.bdEnt(1, DirBd);
NeuFace = msh.bdEnt(2, NeuBd);
%%
%[text] $P\_1$ element
%[text] - Element: tetrahedron
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Continuous $P\_1$ element space
%[text] - $U\_{h, g} = \\{ v \\in H^1 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h; \\ v|\_{\\Gamma\_D} = g \\}$ \
%[text] Base functions are generated on the reference element and mapped affinely (`"map", "affine"`).
P1_FE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_u), "map", "affine");
P1_BC = BC(g_D, "node", DirNode);
P1_FES = FES(msh, P1_FE, P1_BC);
%%
%[text] $P\_2$ element
%[text] - Element: tetrahedron
%[text] - Function space: $P\_2$
%[text] - Nodal DoF: $v|\_c$ at element vertex and edge midpoint \
%[text] Continuous $P\_2$ element space
%[text] - $U\_{h, g} = \\{ v \\in H^1 (\\Omega) : v|\_{K} \\in P\_2 (K), \\forall K \\in T\_h; \\ v|\_{\\Gamma\_D} = g \\}$ \
P2_FE = FE("D3T", "[1,x,y,z,x^2,y^2,z^2,x*y,y*z,z*x]", [NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_u), ...
    NdDoF("D3", MshEnt("D3T").msh, 1, 0.5, d0_u)], "map", "affine");
P2_BC = BC(g_D, "node", DirNode, "edge", DirEdge);
P2_FES = FES(msh, P2_FE, P2_BC);
%%
%[text] CG scheme: find $u\_h \\in U\_{h, g}$ such that
%[text] $\\int\_{\\Omega} \\nabla u\_h \\cdot \\nabla v\_h \\, dx = \\int\_{\\Omega} f v\_h \\, dx + \\int\_{\\Gamma\_N} g\_N v\_h \\, ds , \\quad \\forall v\_h \\in U\_{h, 0}$
Uh = P1_FES; ord = 1;
% Uh = P2_FES; ord = 2;

d0_v = d0_u; grad_v = grad_u;
Auv = DLF(msh, 3, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D3T", (ord - 1) * 2));
Fv = [SLF(msh, 3, f, d0_v, "GInt", GInt("D3T", 1 + ord)), ...
    SLF(msh, 2, g_N, d0_v, "EntIdx", NeuFace, "GInt", GInt("D3F", 1 + ord))];
%%
%[text] Solution
[Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv);
uh = FEF(Uh, Stiff \ Load);
%%
%[text] Error
%[text] $\\| v \\|\_{L^2} = \\Big( \\int\_{\\Omega} |v|^2 \\, dx \\Big)^{1/2} \\\\\n\\| v \\|\_{H^1} = \\Big( \\int\_{\\Omega} |\\nabla v|^2 \\, dx \\Big)^{1/2}$
L2Norm = Norm(msh, 3, d0_u, "GInt", GInt("D3T", ord * 2));
H1Norm = [Norm(msh, 3, grad_u, "GInt", GInt("D3T", ord * 2))];
fprintf('|u-uh|_L2: %e, |u-uh|_H1: %e\n', eNorm(msh, u, uh, L2Norm), eNorm(msh, u, uh, H1Norm));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
