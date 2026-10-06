%[text] Poisson equation in mixed form, HDG method (3D)
%[text] $\\begin{cases}\n  \\mathbf{u} = \\nabla p & \\text{in } \\Omega \\\\\n  - \\nabla \\cdot \\mathbf{u} = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on } \\Gamma\_D \\\\\n  \\mathbf{u} \\cdot \\mathbf{n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
p = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z)+x*y*z");
d0_p = [0; 0; 0]; dx_p = [1; 0; 0]; dy_p = [0; 1; 0]; dz_p = [0; 0; 1]; grad_p = cat(3, dx_p, dy_p, dz_p);
u = dif(p, grad_p);
d0_u = zeros(3); div_u = eye(3);
f = -sum(dif(u, div_u));
g_D = p;
UNV = MshEnt("D3F").UNV;
g_N = dot(u, UNV);
%%
%[text] Domain: cube $\[0, 1\] \\times \[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured tetrahedral mesh (each cube split into 6 tetrahedra)
%[text] Dirichlet boundary: faces $x = 0$, $y = 0$, $z = 0$
%[text] Neumann boundary: faces $x = 1$, $y = 1$, $z = 1$
msh = mshD3TS([0, 1, 0, 1, 0, 1], 4);
DirBd = [1, 3, 5];
NeuBd = [2, 4, 6];
DirFace = msh.bdEnt(2, DirBd);
NeuFace = msh.bdEnt(2, NeuBd);
%%
%[text] Discontinuous $P\_1$ element (scalar)
SP1_FE = FE("D3T", "[1,x,y,z]", NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_p, "share", false), "map", "affine");
SP1_FES = FES(msh, SP1_FE);
%%
%[text] Discontinuous $P\_1$ element (vector)
d0_ux = [0, nan, nan; 0, nan, nan; 0, nan, nan];
d0_uy = [nan, 0, nan; nan, 0, nan; nan, 0, nan];
d0_uz = [nan, nan, 0; nan, nan, 0; nan, nan, 0];
VP1_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), [NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_ux, "share", false), ...
    NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_uy, "share", false), ...
    NdDoF("D3", MshEnt("D3T").msh, 0, [], d0_uz, "share", false)], "map", "affine");
VP1_FES = FES(msh, VP1_FE);
%%
%[text] Trace $P\_1$ element
%[text] - Element: reference triangle (face), parameterized by global face vertices
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at face vertex \
d0_l = [0; 0];
TP1_FE = FE("D3FR", "[1,s,t]", NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], d0_l, "share", false));
TP1_BC = BC(g_D, "face", DirFace);
TP1_FES = FES(msh, TP1_FE, TP1_BC);
%%
%[text] HDG scheme: find $\\mathbf{u}\_h \\in U\_h$, $p\_h \\in P\_h$ and $\\lambda\_h \\in \\Lambda\_{h, g}$ such that
%[text] $\\begin{cases}\n  \\int\_\\Omega \\mathbf{u}\_h \\cdot \\mathbf{v}\_h \\, dx + \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla \\cdot \\mathbf{v}\_h \\, dx - \\sum\_{K \\in T\_h} \\int\_{\\partial K} \\lambda\_h \\mathbf{v}\_h \\cdot \\mathbf{n}\_K \\, ds = 0, & \\forall \\mathbf{v}\_h \\in U\_h \\\\\n  \\sum\_{K \\in T\_h} \\int\_K \\mathbf{u}\_h \\cdot \\nabla q\_h \\, dx - \\sum\_{K \\in T\_h} \\int\_{\\partial K} \\hat{\\mathbf{u}}\_h \\cdot \\mathbf{n}\_K \\, q\_h \\, ds = \\int\_\\Omega f q\_h \\, dx, & \\forall q\_h \\in P\_h \\\\\n  \\sum\_{K \\in T\_h} \\int\_{\\partial K} \\hat{\\mathbf{u}}\_h \\cdot \\mathbf{n}\_K \\, \\mu\_h \\, ds = \\int\_{\\Gamma\_N} g\_N \\mu\_h \\, ds, & \\forall \\mu\_h \\in \\Lambda\_{h, 0}\n\\end{cases}$
%[text] Numerical flux on each element: $\\hat{\\mathbf{u}}\_h \\cdot \\mathbf{n}\_K = \\mathbf{u}\_h|\_K \\cdot \\mathbf{n}\_K - \\tau (p\_h|\_K - \\lambda\_h)$ ($\\tau > 0$; minus sign because $\\mathbf{u} = \\nabla p$).
%[text] On a face with positive and negative elements, $\\mathbf{n}\_K = \\mathbf{n}$ and $-\\mathbf{n}$ respectively: sign of `iTrl / iTst` selects the side.
tau = 1;
Uh = VP1_FES; ord_Uh = 1;
Ph = SP1_FES; ord_Ph = 1;
Lh = TP1_FES; ord_Lh = 1;
IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, [DirBd, NeuBd]));

trls = [Uh, Ph, Lh]; tsts = [Uh, Ph, Lh];
iu = 1; ip = 2; il = 3; iv = 1; iq = 2; im = 3;
d0_v = d0_u; div_v = div_u; d0_q = d0_p; grad_q = grad_p; d0_m = d0_l;

Auv = DLF(msh, 3, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D3T", ord_Uh * 2));
Bpv = DLF(msh, 3, Fcn.cst(1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D3T", ord_Ph + ord_Uh - 1));
Blv = DLF.interface(msh, 2, -UNV, d0_l, d0_v, "iTrl", il, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D3F", ord_Lh + ord_Uh));

Buq = [DLF(msh, 3, Fcn.cst(1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3T", ord_Uh + ord_Ph - 1)), ...
    DLF(msh, 2, -UNV, d0_u, d0_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3F", ord_Uh + ord_Ph)), ...
    DLF(msh, 2, UNV, d0_u, d0_q, "iTrl", -iu, "iTst", -iq, "GInt", GInt("D3F", ord_Uh + ord_Ph))];
Apq = [DLF(msh, 2, Fcn.cst(tau), d0_p, d0_q, "iTrl", ip, "iTst", iq, "GInt", GInt("D3F", ord_Ph * 2)), ...
    DLF(msh, 2, Fcn.cst(tau), d0_p, d0_q, "iTrl", -ip, "iTst", -iq, "GInt", GInt("D3F", ord_Ph * 2))];
Blq = [DLF(msh, 2, Fcn.cst(-tau), d0_l, d0_q, "iTrl", il, "iTst", iq, "GInt", GInt("D3F", ord_Lh + ord_Ph)), ...
    DLF(msh, 2, Fcn.cst(-tau), d0_l, d0_q, "iTrl", il, "iTst", -iq, "GInt", GInt("D3F", ord_Lh + ord_Ph))];
Fq = SLF(msh, 3, f, d0_q, "iTst", iq, "GInt", GInt("D3T", 1 + ord_Ph));

Bum = DLF.interface(msh, 2, UNV, d0_u, d0_m, "iTrl", iu, "iTst", im, "trlOpr", "jump", "GInt", GInt("D3F", ord_Uh + ord_Lh));
Bpm = [DLF(msh, 2, Fcn.cst(-tau), d0_p, d0_m, "iTrl", ip, "iTst", im, "GInt", GInt("D3F", ord_Ph + ord_Lh)), ...
    DLF(msh, 2, Fcn.cst(-tau), d0_p, d0_m, "iTrl", -ip, "iTst", im, "GInt", GInt("D3F", ord_Ph + ord_Lh))];
% Number of elements sharing the face: 1 on boundary, 2 in interior.
Alm = [DLF(msh, 2, Fcn.cst(tau), d0_l, d0_m, "iTrl", il, "iTst", im, "EntIdx", NeuFace, "GInt", GInt("D3F", ord_Lh * 2)), ...
    DLF(msh, 2, Fcn.cst(2 * tau), d0_l, d0_m, "iTrl", il, "iTst", im, "EntIdx", IntFace, "GInt", GInt("D3F", ord_Lh * 2))];
Gm = SLF(msh, 2, g_N, d0_m, "iTst", im, "EntIdx", NeuFace, "GInt", GInt("D3F", 1 + ord_Lh));
%%
%[text] Solution: element-local DoFs of $U\_h$ and $P\_h$ are eliminated by static condensation, only the trace system is solved.
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Blv, Buq, Apq, Blq, Bum, Bpm, Alm], [Fq, Gm]);
[uh, ph, lh] = FEF.multi(trls, condSolve(Stiff, Load, trls, [iu, ip]));
%%
%[text] Error
L2Norm_Uh = Norm(msh, 3, d0_u, "GInt", GInt("D3T", (ord_Uh + 1) * 2));
L2Norm_Ph = Norm(msh, 3, d0_p, "GInt", GInt("D3T", (ord_Ph + 1) * 2));
L2Norm_Lh = Norm(msh, 2, d0_l, "GInt", GInt("D3F", 4));
fprintf('|u-uh|_L2: %e, |p-ph|_L2: %e, |l-l_h|_L2: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, p, ph, L2Norm_Ph), eNorm(msh, p, lh, L2Norm_Lh));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
