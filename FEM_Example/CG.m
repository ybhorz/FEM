%[text] Poisson equation:
%[text] $\\begin{cases}  - \\Delta u = f & \\text{in } \\Omega, \\\\  u = g\_D & \\text{on \n} \\Gamma\_D, \\\\  \\frac{\\partial u}{\\partial n} = g\_N & \\text{on } \\Gamma\_N.\\end{cases}$
u = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y");
d0_u = [0; 0]; dx_u = [1; 0]; dy_u = [0; 1]; grad_u = cat(3, dx_u, dy_u); dxx_u = [2; 0]; dyy_u = [0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u);
g_D = u;
g_N = dot(MshEnt("D2L").UNV, dif(u, grad_u));
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$.
%[text] Dirichlet boundary: lower and right boundary.
%[text] Neumann boundary: upper and left boundary.
%[text] Mesh: structured triangulation.
msh = mshD2TS([0, 1, 0, 1], 16);
DirBd = [1, 2, -1, -2, -3];
NeuBd = [3, 4];
DirNode = find(ismember(msh.node.type, DirBd));
DirEdge = find(ismember(msh.edge.type, DirBd));
NeuEdge = find(ismember(msh.edge.type, NeuBd));
%%
%[text] P1 Lagrange element.
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
P1_BC = BC(g_D, "node", DirNode);
P1_FES = FES(msh, P1_FE, P1_BC);
%%
%[text] P2 Lagrange element.
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 0.5, [0; 0])]);
P2_BC = BC(g_D, "node", DirNode, "edge", DirEdge);
P2_FES = FES(msh, P2_FE, P2_BC);
%%
%[text] FEM scheme: find $u\_h \\in U\_h$ such that
%[text] $\\int\_{\\Omega} \\nabla u\_h \\cdot \\nabla v\_h \\, dx = \\int\_{\\Omega} f v\_h \\, \ndx + \\int\_{\\Gamma\_N} g\_N v\_h \\, ds , \\quad \\forall v\_h \\in U\_h.$
Uh = P1_FES; ord = 1;
% Uh = P2_FES; ord = 2;

d0_v = d0_u; grad_v = grad_u;
Auv = DLF(msh, 2, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D2T", (ord - 1) * 2));
Fv = [SLF(msh, 2, f, d0_v, "GInt", GInt("D2T", 1 + ord)), ...
    SLF(msh, 1, g_N, d0_v, "EntIdx", NeuEdge, "GInt", GInt("D2L", 1 + ord))];
%%
%[text] Solution.
[Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv);
uh = FEF(Uh, Stiff \ Load);
%%
%[text] Error.
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
H1Norm = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", (ord - 1) * 2))];
fprintf('L2 norm of error: %e, H1 norm of error: %e\n', eNorm(msh, u, uh, L2Norm), eNorm(msh, u, uh, H1Norm));
%%
%[text] Plot.
plotFEF(uh);
plotFcn(msh, u);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
