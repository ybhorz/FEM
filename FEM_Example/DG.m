%[text] Poisson equation:
%[text] $\\begin{cases}  - \\Delta u = f & \\text{in } \\Omega, \\\\  u = g\_D & \\text{on\n} \\Gamma\_D, \\\\  \\frac{\\partial u}{\\partial n} = g\_N & \\text{on } \\Gamma\_N.\\end{cases}$
u = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y");
d0_u = [0; 0]; dx_u = [1; 0]; dy_u = [0; 1]; grad_u = cat(3, dx_u, dy_u); dxx_u = [2; 0]; dyy_u = [0; 2];
f =- dif(u, dxx_u) - dif(u, dyy_u);
g_D = u;
UNV = MshEnt("D2L").UNV;
g_N = dot(UNV, dif(u, grad_u));
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
IntEdge = find(ismember(msh.edge.type, 0));
%%
%[text] DG P1 element.
DG1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "share", false));
DG1_BC = BC(g_D, "node", DirNode);
DG1_FES = FES(msh, DG1_FE, DG1_BC);
%%
%[text] DG scheme (SIPG): find $u\_h \\in U\_h$ such that
%[text] $ \\sum\_{K \\in T\_h} \\int\_{K} \\nabla u\_h \\cdot \\nabla v\_h \\, dx - \\sum\_{e \\in\n\\Gamma^o\_h} \\int\_{e} \\{ \\frac{\\partial u\_h}{\\partial n} \\} \[ v\_h \] \\, ds - \\sum\_{e\n\\in \\Gamma^o\_h} \\int\_{e} \\{ \\frac{\\partial v\_h}{\\partial n} \\} \[ u\_h \] \\, ds\n+ \\sum\_{e \\in \\Gamma^o\_h} \\int\_{e} \\frac{\\gamma}{h\_e} \[ u\_h \] \[ v\_h \] \\, ds\n=  \\sum\_{K \\in T\_h} \\int\_{K} f v\_h \\, dx + \\sum\_{e \\in \\Gamma^N\_h} \\int\_e g\_N\nv\_h \\, ds \\quad \\forall v\_h \\in U\_h.$
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
%[text] $|| e ||\_{H^1, h}^2 =  \\sum\_{K \\in T\_h} \\int\_{K}  | \\nabla e |^2 \\, dx +\n\\sum\_{e \\in \\Gamma^o\_h} \\int\_{e} \\frac{1}{h\_e} \[ e \]^2 \\, ds.$
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
H1Norm = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", (ord - 1) * 2)), ...
    Norm(msh, 1, d0_u, "coef", MshEnt("D2L").len \ 1, "fcnOpr", "jump", "GInt", GInt("D2L", ord * 2))];
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
