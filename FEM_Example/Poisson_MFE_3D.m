%[text] Poisson equation in mixed form (3D)
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
%[text] $RT\_0$ element
%[text] - Element: tetrahedron
%[text] - Function space: $\[P\_0\]^3 + \\mathbf{x} P\_0$
%[text] - Moment DoF: $\\int\_F \\mathbf{v} \\cdot \\mathbf{n} \\, ds$ on each face \
%[text] Base functions are generated on the reference element and mapped by the contravariant Piola transformation (`"map", "piolaDiv"`).
RT0_FE = FE("D3T", "[1,0,0; 0,1,0; 0,0,1; x,y,z].'", ...
    MoDoF("D3", MshEnt("D3T").msh, 2, Fcn.cst(1), d0_u, "coef", UNV, "orien", true, "GInt", GInt("D3F", 1)), "map", "piolaDiv");
RT0_BC = BC(g_N, "face", NeuFace);
RT0_FES = FES(msh, RT0_FE, RT0_BC);
%%
%[text] $BDM\_1$ element
%[text] - Element: tetrahedron
%[text] - Function space: $\[P\_1\]^3$
%[text] - Moment DoF: $\int\_F \mathbf{v} \cdot \mathbf{n} \, \lambda\_i \, ds$ on each face, $\lambda\_i$: barycentric coordinates of face \
BDM1_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), ...
    MoDoF("D3", MshEnt("D3T").msh, 2, [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")], d0_u, "coef", UNV, ...
    "orien", true, "GInt", GInt("D3F", 2)), "map", "piolaDiv");
BDM1_BC = BC(g_N, "face", NeuFace);
BDM1_FES = FES(msh, BDM1_FE, BDM1_BC);
%%
%[text] $P\_0$ element
P0_FE = FE("D3T", "1", NdDoF("D3", MshEnt("D3T").msh, 3, [1/4; 1/4; 1/4], d0_p), "map", "affine");
P0_FES = FES(msh, P0_FE);
%%
%[text] MFE scheme: find $(\\mathbf{u}\_h, p\_h) \\in U\_{h, g} \\times P\_h$ such that
%[text] $\\begin{cases}\n  \\int\_{\\Omega} \\mathbf{u}\_h \\cdot \\mathbf{v}\_h \\, dx + \\int\_{\\Omega} p\_h \\nabla \\cdot \\mathbf{v}\_h \\, dx = \\int\_{\\Gamma\_D} g\_D \\mathbf{v}\_h \\cdot \\mathbf{n} \\, ds, & \\forall \\mathbf{v}\_h \\in U\_{h, 0} \\\\\n  - \\int\_{\\Omega} \\nabla \\cdot \\mathbf{u}\_h q\_h \\, dx = \\int\_{\\Omega} f q\_h \\, dx, & \\forall q\_h \\in P\_h\n\\end{cases}$
Uh = RT0_FES; ord_Uh = 1;
Ph = P0_FES; ord_Ph = 0;

% Uh = BDM1_FES; ord_Uh = 1;
% Ph = P0_FES; ord_Ph = 0;

trls = [Uh, Ph]; tsts = [Uh, Ph];
iu = 1; ip = 2; iv = 1; iq = 2;
d0_v = d0_u; div_v = div_u; d0_q = d0_p;
Auv = DLF(msh, 3, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D3T", ord_Uh * 2));
Bpv = DLF(msh, 3, Fcn.cst(1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D3T", ord_Ph + ord_Uh - 1));
Buq = DLF(msh, 3, Fcn.cst(-1), div_u, d0_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3T", ord_Uh - 1 + ord_Ph));
Gv = SLF(msh, 2, g_D .* UNV, d0_v, "iTst", iv, "EntIdx", DirFace, "GInt", GInt("D3F", 1 + ord_Uh));
Fq = SLF(msh, 3, f, d0_q, "iTst", iq, "GInt", GInt("D3T", 1 + ord_Ph));
%%
%[text] Solution
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Buq], [Gv, Fq]);
[uh, ph] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error
L2Norm_Uh = Norm(msh, 3, d0_u, "GInt", GInt("D3T", (ord_Uh + 1) * 2));
divNorm_Uh = Norm(msh, 3, div_u, "form",  @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow, "GInt", GInt("D3T", ord_Uh * 2));
L2Norm_Ph = Norm(msh, 3, d0_p, "GInt", GInt("D3T", (ord_Ph + 1) * 2));
fprintf('|u-uh|_L2: %e, |u-uh|_div: %e, |p-ph|_L2: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, divNorm_Uh), eNorm(msh, p, ph, L2Norm_Ph));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
