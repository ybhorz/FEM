%[text] Brinkman equation (3D)
%[text] $\\begin{cases}\n  - \\nu \\Delta \\mathbf{u} + \\alpha \\mathbf{u} + \\nabla p = \\mathbf{f} &\\text{in } \\Omega \\\\\n  \\nabla \\cdot \\mathbf{u} = 0 &\\text{in } \\Omega \\\\\n  \\mathbf{u} = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  \\sigma - \\nabla \\mathbf{u} = 0 &\\text{in } \\Omega \\\\\n  -\\nu \\nabla \\cdot \\sigma + \\alpha \\mathbf{u} + \\nabla p = \\mathbf{f} &\\text{in } \\Omega \\\\\n  \\nabla \\cdot \\mathbf{u} = 0 &\\text{in } \\Omega \\\\\n  \\mathbf{u} = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
%[text] Exact solution: $\\mathbf{u} = \\nabla \\times (\\phi, \\phi, \\phi)$ with $\\phi = (\\sin \\pi x \\sin \\pi y \\sin \\pi z)^2$ (divergence free, zero on boundary)
v = 1; a = 10;
d0_p = [0; 0; 0]; grad_p = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
phi = Fcn("D3", "sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2").dif(grad_p).fun;
u = Fcn("D3", [phi(2) - phi(3); phi(3) - phi(1); phi(1) - phi(2)]);
p = Fcn("D3", "10*(x-1/2)*(y-1/2)*(z-1/2)");
d0_u = zeros(3); div_u = eye(3);
grad_u = cat(3, [1, 1, 1; 0, 0, 0; 0, 0, 0], [0, 0, 0; 1, 1, 1; 0, 0, 0], [0, 0, 0; 0, 0, 0; 1, 1, 1]);
s = dif(u, grad_u);
d0_s = zeros(3, 9); div_s = repelem(eye(3), 1, 3);
f = - sum(dif(s, div_s), 2) * v + u * a + dif(p, grad_p);
d0_g = zeros(2, 2); d0_l = zeros(2, 3);

UNV = MshEnt("D3F").UNV;
FcParm = MshEnt("D3F").parm;
% Tangent vectors of face (columns, from global face vertices): tangential field T * [g1; g2] on face.
Tmat = Fcn("D3F", [FcParm(:, 2) - FcParm(:, 1), FcParm(:, 3) - FcParm(:, 1)]);
hF = MshEnt("D3F").area .^ (1/2);
%%
%[text] Domain: cube $\[0, 1\] \\times \[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured tetrahedral mesh + Alfeld splitting
%[text] Sub-element: tetrahedron $\[v\_1, v\_2, v\_3, c\]$, face 4 is a primal face, faces 1, 2, 3 are dual faces.
msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], 4));
PrFace = find(ismember(msh.face.type, 0:6));
PrOFace = find(ismember(msh.face.type, 0));
DlFace = find(ismember(msh.face.type, 1i));
BdFace = msh.bdEnt(2, 1:6);
D3TElem = MshEnt("D3T").msh;
FcBar = [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")];
%%
%[text] Scalar $SDG\_1$ element
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $q|\_c$ at vertices of face 4 (shared) and at vertex 4 \
%[text] Scalar $SDG\_1$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_1 (K), \\forall K \\in T\_h; \\ \[q\]\_F = 0, \\forall F \\in F^{pr, o}\_h \\}$ (pressure fixed at one DoF) \
SDG_1S_FE = FE("D3T", "[1,x,y,z]", [NdDoF("D3", D3TElem, 2, [0, 1, 0; 0, 0, 1], d0_p, "EntIdx", 4), ...
    NdDoF("D3", D3TElem, 0, [], d0_p, "EntIdx", 4, "share", false)], "map", "affine");
SDG_1S_BC = BC(p, "DoF", 1);
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);
%[text] Vector $SDG\_1$ element
%[text] - Function space: $P^3\_1$
%[text] - Moment DoF: $\\int\_F \\mathbf{v} \\cdot \\mathbf{n} \\, \\lambda\_k \\, ds$ on faces 1, 2, 3 (shared) and on face 4, $\\lambda\_k$: barycentric coordinates of face (BDM1, mapped by contravariant Piola transformation) \
%[text] Vector $SDG\_1$ element space
%[text] - $U\_h = \\{ \\mathbf{v} \\in L^2 (\\Omega, R^3) : \\mathbf{v}|\_K \\in P\_1 (K)^3, \\forall K \\in T\_h; \\ \[\\mathbf{v} \\cdot \\mathbf{n}\]\_F = 0, \\forall F \\in F^{dl}\_h \\}$ \
SDG_1V_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), ...
    [MoDoF("D3", D3TElem, 2, FcBar, d0_u, "coef", UNV, "EntIdx", [1, 2, 3], "orien", true, "GInt", GInt("D3F", 2)), ...
    MoDoF("D3", D3TElem, 2, FcBar, d0_u, "coef", UNV, "EntIdx", 4, "orien", true, "share", false, "GInt", GInt("D3F", 2))], ...
    "map", "piolaDiv");
SDG_1V_FES = FES(msh, SDG_1V_FE);
%[text] Matrix $DG\_1$ element (the SDG matrix element has no simple basis in 3D, its continuity is imposed weakly below)
%[text] - Function space: $P^{3 \\times 3}\_1$
%[text] - Nodal DoF: each component at vertices (not shared) \
%[text] Matrix $DG\_1$ element space
%[text] - $\\Sigma\_h = \\{ \\tau \\in L^2 (\\Omega, R^{3 \\times 3}) : \\tau|\_K \\in P\_1 (K)^{3 \\times 3}, \\forall K \\in T\_h \\}$ \
DG_1M_DoFs = NdDoF.empty;
for iComp = 1:9
    d0_si = nan(3, 9); d0_si(:, iComp) = 0;
    DG_1M_DoFs(end + 1) = NdDoF("D3", D3TElem, 0, [], d0_si, "share", false);
end
DG_1M_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 3]), DG_1M_DoFs, "map", "affine");
DG_1M_FES = FES(msh, DG_1M_FE);
%[text] Multipliers (trace spaces on all faces; faces where a multiplier is not used are set to zero)
%[text] - Tangential velocity $T \\gamma\_h$ on dual faces: $\\gamma\_h$ with 2 components, linear, each component at face vertices
%[text] - Velocity $\\lambda\_h$ on interior primal faces: 3 components, linear, each component at face vertices ($\\lambda\_h = 0$ on $\\partial \\Omega$ imposes $\\mathbf{u} = 0$) \
TP_1_DoFs = @(nComp) arrayfun(@(iComp) NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], zeros(2, nComp), ...
    "coef", Fcn.cst(double(1:nComp == iComp)'), "share", false), 1:nComp);
TP_1V2_FES = FES(msh, FE("D3FR", FE.repFS("[1,s,t]", [2, 1]), TP_1_DoFs(2)), BC(Fcn("D3", "0"), "face", PrFace));
TP_1V3_FES = FES(msh, FE("D3FR", FE.repFS("[1,s,t]", [3, 1]), TP_1_DoFs(3)), BC(Fcn("D3", "0"), "face", [DlFace, BdFace]));
%%
%[text] Hybridized SDG scheme (Zhao, Chung and Lam, CMAME 364, 2020): $\\sigma\_h$ is discontinuous; $\[\\sigma\_h \\mathbf{n}\] = 0$ on interior primal faces and $\[\\sigma\_h \\mathbf{n}\] \\times \\mathbf{n} = 0$ on dual faces are imposed weakly by $\\lambda\_h$ and $\\gamma\_h$. Find $\\sigma\_h \\in \\Sigma\_h$, $\\mathbf{u}\_h \\in U\_h$, $\\gamma\_h$, $\\lambda\_h$, and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}\n  \\sum\_{K} \\int\_K \\sigma\_h : \\tau\_h \\, dx + \\sum\_{K} \\int\_K \\mathbf{u}\_h \\cdot \\nabla \\cdot \\tau\_h \\, dx - \\sum\_{F \\in F\_h^{dl}} \\int\_F ((\\mathbf{u}\_h \\cdot \\mathbf{n}) \\mathbf{n} + T \\gamma\_h) \\cdot \[\\tau\_h \\mathbf{n}\] \\, ds - \\sum\_{F \\in F\_h^{pr, o}} \\int\_F \\lambda\_h \\cdot \[\\tau\_h \\mathbf{n}\] \\, ds = 0 &\\forall \\tau\_h \\in \\Sigma\_h \\\\\n  \\alpha \\sum\_{K} \\int\_K \\mathbf{u}\_h \\cdot \\mathbf{v}\_h \\, dx - \\nu \\sum\_{K} \\int\_K \\nabla \\cdot \\sigma\_h \\cdot \\mathbf{v}\_h \\, dx + \\nu \\sum\_{F \\in F\_h^{dl}} \\int\_F \[\\sigma\_h \\mathbf{n}\] \\cdot \\mathbf{n} \\, \\mathbf{v}\_h \\cdot \\mathbf{n} \\, ds - \\sum\_{K} \\int\_K p\_h \\nabla \\cdot \\mathbf{v}\_h \\, dx + \\sum\_{F \\in F\_h^{pr}} \\int\_F p\_h \[\\mathbf{v}\_h \\cdot \\mathbf{n}\] \\, ds = \\sum\_{K} \\int\_K \\mathbf{f} \\cdot \\mathbf{v}\_h \\, dx &\\forall \\mathbf{v}\_h \\in U\_h \\\\\n  \\nu \\sum\_{F \\in F\_h^{dl}} \\int\_F \[\\sigma\_h \\mathbf{n}\] \\cdot T \\eta\_h \\, ds = 0, \\quad \\nu \\sum\_{F \\in F\_h^{pr, o}} \\int\_F \[\\sigma\_h \\mathbf{n}\] \\cdot \\mu\_h \\, ds = 0 &\\forall \\eta\_h, \\mu\_h \\\\\n  - \\sum\_{K} \\int\_K \\mathbf{u}\_h \\cdot \\nabla q\_h \\, dx + \\sum\_{F \\in F\_h^{dl}} \\int\_F \\mathbf{u}\_h \\cdot \\mathbf{n} \[ q\_h \] \\, ds = 0 &\\forall q\_h \\in P\_h\n\\end{cases}$
Sh = DG_1M_FES; ord_Sh = 1;
Uh = SDG_1V_FES; ord_Uh = 1;
Gh = TP_1V2_FES; ord_Gh = 1;
Lh = TP_1V3_FES; ord_Lh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;

trls = [Sh, Uh, Gh, Lh, Ph]; tsts = [Sh, Uh, Gh, Lh, Ph];
is = 1; iu = 2; ig = 3; il = 4; ip = 5; it = 1; iv = 2; ih = 3; im = 4; iq = 5;

d0_t = d0_s; div_t = div_s; d0_v = d0_u; div_v = div_u; grad_v = grad_u; d0_h = d0_g; d0_m = d0_l; d0_q = d0_p; grad_q = grad_p;
Auv = DLF(msh, 3, Fcn.cst(a), d0_u, d0_v, "iTrl", iu, "iTst", iv, "GInt", GInt("D3T", ord_Uh * 2));
Ast = DLF(msh, 3, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D3T", ord_Sh * 2));
But = [DLF(msh, 3, Fcn.cst(1), d0_u, div_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D3T", ord_Uh + ord_Sh - 1)), ...
    DLF.interface(msh, 2, [-UNV, UNV, UNV], d0_u, d0_t, "EntIdx", DlFace, "iTrl", iu, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, coef(1)) .* dot(tst * coef(2), coef(3)), "GInt", GInt("D3F", ord_Uh + ord_Sh))];
Bgt = DLF.interface(msh, 2, [-UNV, Tmat], d0_g, d0_t, "EntIdx", DlFace, "iTrl", ig, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(coef(2) * trl, tst * coef(1)), "GInt", GInt("D3F", ord_Gh + ord_Sh));
Blt = DLF.interface(msh, 2, -UNV, d0_l, d0_t, "EntIdx", PrOFace, "iTrl", il, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D3F", ord_Lh + ord_Sh));
Bsv = [DLF(msh, 3, Fcn.cst(-v), div_s, d0_v, "iTrl", is, "iTst", iv, "form", @(coef, trl, tst) coef .* dot(sum(trl, 2), tst), "GInt", GInt("D3T", ord_Sh + ord_Uh - 1)), ...
    DLF.interface(msh, 2, [UNV * v, UNV, UNV], d0_s, d0_v, "EntIdx", DlFace, "iTrl", is, "iTst", iv, "trlOpr", "jump", "form", @(coef, trl, tst) dot(tst, coef(1)) .* dot(trl * coef(2), coef(3)), "GInt", GInt("D3F", ord_Sh + ord_Uh))];
Bsh = DLF.interface(msh, 2, [UNV * v, Tmat], d0_s, d0_h, "EntIdx", DlFace, "iTrl", is, "iTst", ih, "trlOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef(1), coef(2) * tst), "GInt", GInt("D3F", ord_Sh + ord_Gh));
Bsm = DLF.interface(msh, 2, UNV * v, d0_s, d0_m, "EntIdx", PrOFace, "iTrl", is, "iTst", im, "trlOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D3F", ord_Sh + ord_Lh));
Bpv = [DLF(msh, 3, Fcn.cst(-1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D3T", ord_Ph + ord_Uh - 1)), ...
    DLF.interface(msh, 2, UNV, d0_p, d0_v, "EntIdx", PrFace, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D3F", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 3, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3T", ord_Uh + ord_Ph - 1)), ...
    DLF.interface(msh, 2, UNV, d0_u, d0_q, "EntIdx", DlFace, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D3F", ord_Uh + ord_Ph))];
Fv = SLF(msh, 3, f, d0_v, "iTst", iv, "GInt", GInt("D3T", 2 + ord_Uh));
%%
%[text] Solution
%[text] Local elimination: $\\sigma\_h$ element by element, then $(\\mathbf{u}\_h, \\gamma\_h)$ macro element by macro element (original tetrahedra); the global system is for $(\\lambda\_h, p\_h)$.
[Stiff, Load] = assemble(msh, trls, tsts, [Auv, Ast, But, Bgt, Blt, Bsv, Bsh, Bsm, Bpv, Buq], Fv);
[sh, uh, ~, ~, ph] = FEF.multi(trls, condSolve(Stiff, Load, trls, {is, [iu, ig]}));
%%
%[text] Error
%[text] $\\| \\mathbf{v} \\|\_{H^1, h}^2 = \\sum\_{K} \\| \\nabla \\mathbf{v} \\|\_K^2 + \\sum\_{F \\in F\_h} h^{-1}\_F \\| \[\\mathbf{v}\] \\|\_F^2$ (on dual faces only the tangential part jumps)
L2Norm_Sh = Norm(msh, 3, d0_s, "GInt", GInt("D3T", (1 + ord_Sh) * 2));
L2Norm_Uh = Norm(msh, 3, d0_u, "GInt", GInt("D3T", (1 + ord_Uh) * 2));
H1Norm_Uh = [Norm(msh, 3, grad_u, "GInt", GInt("D3T", ord_Uh * 2)), ...
    Norm(msh, 2, d0_u, "coef", hF \ 1, "form", @(coef, fcn, pow) coef .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D3F", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 3, d0_p, "GInt", GInt("D3T", (1 + ord_Ph) * 2));
% Pressure is determined up to a constant (fixed at one DoF): error is measured with mean removed,
% |p - ph - c|_L2^2 = |p - ph|_L2^2 - |Omega| c^2, c = mean(p - ph), |Omega| = 1.
MeanNorm_Ph = Norm(msh, 3, d0_p, "pow", 1, "form", @(coef, fcn, pow) fcn, "GInt", GInt("D3T", (1 + ord_Ph) * 2));
eMean_Ph = eNorm(msh, p, ph, MeanNorm_Ph);
fprintf('|s-sh|_L2: %e, |u-uh|_L2: %e, |u-uh|_H1: %e, |p-ph|_L2: %e\n', eNorm(msh, s, sh, L2Norm_Sh), eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, H1Norm_Uh), sqrt(eNorm(msh, p, ph, L2Norm_Ph) ^ 2 - eMean_Ph ^ 2));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
