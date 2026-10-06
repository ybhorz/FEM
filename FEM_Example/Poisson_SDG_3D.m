%[text] Poisson equation (3D)
%[text] $\\begin{cases}\n  - \\Delta p = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on }\\Gamma\_D \\\\\n  \\frac{\\partial p}{\\partial n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  \\mathbf{u} - \\nabla p = 0 & \\text{in } \\Omega \\\\ \n  -\\nabla \\cdot \\mathbf{u} = f & \\text{in } \\Omega \\\\\n  p = g\_D & \\text{on } \\Gamma\_D \\\\\n  \\mathbf{u} \\cdot \\mathbf{n} = g\_N & \\text{on } \\Gamma\_N\n\\end{cases}$
p = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z)+x*y*z");
d0_p = [0; 0; 0]; dx_p = [1; 0; 0]; dy_p = [0; 1; 0]; dz_p = [0; 0; 1]; grad_p = cat(3, dx_p, dy_p, dz_p);
u = dif(p, grad_p);
d0_u = zeros(3); div_u = eye(3);
f = -sum(dif(u, div_u));
g_D = p;
UNV = MshEnt("D3F").UNV;
g_N = dot(u, UNV);
hF = MshEnt("D3F").area .^ (1/2);
%%
%[text] Domain: cube $\[0, 1\] \\times \[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured tetrahedral mesh + Alfeld splitting (each tetrahedron is split into 4 by its barycenter)
%[text] Primal faces: faces of the original mesh (type 0 or boundary type); dual faces: new faces (type "1i")
%[text] Dirichlet boundary: faces $x = 0$, $y = 0$, $z = 0$
%[text] Neumann boundary: faces $x = 1$, $y = 1$, $z = 1$
msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], 4));
DirBd = [1, 3, 5];
NeuBd = [2, 4, 6];
DirFace = msh.bdEnt(2, DirBd);
NeuFace = msh.bdEnt(2, NeuBd);
PrFace = find(ismember(msh.face.type, 0:6));
PrOFace = find(ismember(msh.face.type, 0));
DlFace = find(ismember(msh.face.type, 1i));
%%
%[text] Sub-element: tetrahedron $\[v\_1, v\_2, v\_3, c\]$, face 4 is a primal face, faces 1, 2, 3 are dual faces.
%[text] Scalar $SDG\_1$ element
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $q|\_c$ at vertices of face 4 (shared), and at vertex 4 \
%[text] Scalar $SDG\_1$ element space
%[text] - $P\_{h, g} = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_1 (K), \\forall K \\in T\_h; \\ \[q\]\_F = 0, \\forall F \\in F^{pr, o}\_h; \\ q|\_{\\Gamma\_D} = g \\}$ \
D3TElem = MshEnt("D3T").msh;
SDG_1S_FE = FE("D3T", "[1,x,y,z]", [NdDoF("D3", D3TElem, 2, [0, 1, 0; 0, 0, 1], d0_p, "EntIdx", 4), ...
    NdDoF("D3", D3TElem, 0, [], d0_p, "EntIdx", 4, "share", false)], "map", "affine");
SDG_1S_BC = BC(g_D, "node", nan, "face", DirFace);
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);
%[text] Vector $SDG\_1$ element
%[text] - Function space: $P^3\_1$
%[text] - Moment DoF: $\\int\_F \\mathbf{v} \\cdot \\mathbf{n} \\, \\lambda\_i \\, ds$ on faces 1, 2, 3 (shared) and on face 4, $\\lambda\_i$: barycentric coordinates of face \
%[text] Vector $SDG\_1$ element space
%[text] - $U\_h = \\{ \\mathbf{v} \\in L^2 (\\Omega; R^3) : \\mathbf{v}|\_K \\in P\_1 (K)^3, \\forall K \\in T\_h; \\ \[\\mathbf{v} \\cdot \\mathbf{n}\]\_F = 0, \\forall F \\in F^{dl}\_h \\}$ \
FcBar = [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")];
SDG_1V_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), ...
    [MoDoF("D3", D3TElem, 2, FcBar, d0_u, "coef", UNV, "EntIdx", [1, 2, 3], "orien", true, "GInt", GInt("D3F", 2)), ...
    MoDoF("D3", D3TElem, 2, FcBar, d0_u, "coef", UNV, "EntIdx", 4, "orien", true, "share", false, "GInt", GInt("D3F", 2))], ...
    "map", "piolaDiv");
SDG_1V_FES = FES(msh, SDG_1V_FE);
%%
%[text] Scalar $SDG\_0$ element
%[text] - Function space: $P\_0$
%[text] - Nodal DoF: $q|\_c$ at barycenter of face 4 (shared) \
SDG_0S_FE = FE("D3T", "1", NdDoF("D3", D3TElem, 2, [1/3; 1/3], d0_p, "EntIdx", 4), "map", "affine");
SDG_0S_BC = BC(g_D, "face", DirFace);
SDG_0S_FES = FES(msh, SDG_0S_FE, SDG_0S_BC);
%[text] Vector $SDG\_0$ element
%[text] - Function space: $P^3\_0$
%[text] - Moment DoF: $\\int\_F \\mathbf{v} \\cdot \\mathbf{n} \\, ds$ on faces 1, 2, 3 (shared) \
SDG_0V_FE = FE("D3T", FE.repFS("1", [3, 1]), ...
    MoDoF("D3", D3TElem, 2, Fcn.cst(1), d0_u, "coef", UNV, "EntIdx", [1, 2, 3], "orien", true, "GInt", GInt("D3F", 1)), ...
    "map", "piolaDiv");
SDG_0V_FES = FES(msh, SDG_0V_FE);
%%
%[text] SDG scheme: find $\\mathbf{u}\_h \\in U\_h$ and $p\_h \\in P\_{h, g}$ such that
%[text] $\\begin{cases}\n  \\sum\_{K \\in T\_h} \\int\_K \\mathbf{u}\_h \\cdot \\mathbf{v}\_h \\, dx + \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla\\cdot \\mathbf{v}\_h \\, dx - \\sum\_{F \\in F\_h^{pr}} \\int\_F p\_h \[\\mathbf{v}\_h \\cdot \\mathbf{n}\] \\, ds  = 0, & \\forall \\mathbf{v}\_h \\in U\_h \\\\\n  \\sum\_{K \\in T\_h} \\int\_K \\mathbf{u}\_h \\cdot \\nabla q\_h \\, dx - \\sum\_{F \\in F\_h^{dl}} \\int\_F \\mathbf{u}\_h \\cdot \\mathbf{n} \[q\_h\] \\, ds = \\sum\_{K \\in T\_h} \\int\_K f q\_h \\, dx + \\sum\_{F \\in \\Gamma\_h^N} \\int\_F g\_N q\_h \\, ds, &\\forall q\_h \\in P\_{h, 0}\n\\end{cases}$
Uh = SDG_1V_FES; ord_Uh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;

% Uh = SDG_0V_FES; ord_Uh = 1;
% Ph = SDG_0S_FES; ord_Ph = 1;

trls = [Uh, Ph]; tsts = [Uh, Ph];
iu = 1; ip = 2; iv = 1; iq = 2;
d0_v = d0_u; div_v = div_u; d0_q = d0_p; grad_q = grad_p;
Auv = DLF(msh, 3, Fcn.cst(1), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D3T", ord_Uh * 2));
Bpv = [DLF(msh, 3, Fcn.cst(1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D3T", ord_Ph + ord_Uh)), ...
    DLF.interface(msh, 2, -UNV, d0_p, d0_v, "EntIdx", PrFace, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D3F", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 3, Fcn.cst(1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3T", ord_Uh + ord_Ph)), ...
    DLF.interface(msh, 2, -UNV, d0_u, d0_q, "EntIdx", DlFace, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D3F", ord_Uh + ord_Ph))];
Fq = [SLF(msh, 3, f, d0_q, "iTst", iq, "GInt", GInt("D3T", 2 + ord_Ph)), ...
    SLF(msh, 2, g_N, d0_q, "iTst", iq, "EntIdx", NeuFace, "GInt", GInt("D3F", 2 + ord_Ph))];
%%
%[text] Solution
%[text] DoFs of $\mathbf{u}\_h$ on dual faces are shared only inside a macro element (original tetrahedron), hence the mass matrix of $\mathbf{u}\_h$ is block diagonal and $\mathbf{u}\_h$ is eliminated macro element by macro element (static condensation).
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Bpv, Buq], Fq);
[uh, ph] = FEF.multi(trls, condSolve(Stiff, Load, trls, iu));
%%
%[text] Error
%[text] $|| \\mathbf{v} ||^2\_{div, h} = \\sum\_{K \\in T\_h} || \\nabla \\cdot \\mathbf{v} ||^2\_K + \\sum\_{F \\in F\_h^{pr, o}} h^{-1}\_F || \[\\mathbf{v} \\cdot \\mathbf{n}\]||^2\_F \\\\\n|| q ||^2\_{H^1, h} = \\sum\_{K \\in T\_h} ||\\nabla q ||^2\_K + \\sum\_{F \\in F^{dl}\_h} h^{-1}\_F || \[q\] ||^2\_F$
%[text] For $SDG\_0$ ($\\nabla \\cdot \\mathbf{v}\_h = 0$ and $\\nabla q\_h = 0$ in each element), the broken norms do not converge.
L2Norm_Uh = Norm(msh, 3, d0_u, "GInt", GInt("D3T", (1 + ord_Uh) * 2));
divNorm_Uh = [Norm(msh, 3, div_u, "form", @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow, "GInt", GInt("D3T", ord_Uh * 2)), ...
    Norm(msh, 2, d0_u, "EntIdx", PrOFace, "coef", [hF \ 1, UNV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(coef(2), fcn)) .^ pow, "fcnOpr", "jump", "GInt", GInt("D3F", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 3, d0_p, "GInt", GInt("D3T", (1 + ord_Ph) * 2));
H1Norm_Ph = [Norm(msh, 3, grad_p, "GInt", GInt("D3T", ord_Ph * 2)), ...
    Norm(msh, 2, d0_p, "EntIdx", DlFace, "coef", hF \ 1, "fcnOpr", "jump", "GInt", GInt("D3F", ord_Ph * 2))];
fprintf('|u-uh|_L2: %e, |u-uh|_div: %e, |p-ph|_L2: %e, |p-ph|_H1: %e\n', eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, divNorm_Uh), eNorm(msh, p, ph, L2Norm_Ph), eNorm(msh, p, ph, H1Norm_Ph));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
