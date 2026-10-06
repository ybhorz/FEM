%[text] Stokes equation (3D)
%[text] $\\begin{cases}\n  - \\nu \\Delta \\mathbf{u} + \\nabla p = \\mathbf{f} &\\text{in } \\Omega \\\\\n  \\nabla \\cdot \\mathbf{u} = 0 &\\text{in } \\Omega \\\\\n  \\mathbf{u} = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
%[text] Mixed formulation
%[text] $\\begin{cases}\n  \\sigma - \\nabla \\mathbf{u} = 0 &\\text{in } \\Omega \\\\\n  -\\nu \\nabla \\cdot \\sigma + \\nabla p = \\mathbf{f} &\\text{in } \\Omega \\\\\n  \\nabla \\cdot \\mathbf{u} = 0 &\\text{in } \\Omega \\\\\n  \\mathbf{u} = 0 &\\text{on } \\partial \\Omega\n\\end{cases}$
%[text] Exact solution: $\\mathbf{u} = \\nabla \\times (\\phi, \\phi, \\phi)$ with $\\phi = (\\sin \\pi x \\sin \\pi y \\sin \\pi z)^2$ (divergence free, zero on boundary)
v = 1;
d0_p = [0; 0; 0]; grad_p = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
phi = Fcn("D3", "sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2").dif(grad_p).fun;
u = Fcn("D3", [phi(2) - phi(3); phi(3) - phi(1); phi(1) - phi(2)]);
p = Fcn("D3", "10*(x-1/2)*(y-1/2)*(z-1/2)");
d0_u = zeros(3); div_u = eye(3);
grad_u = cat(3, [1, 1, 1; 0, 0, 0; 0, 0, 0], [0, 0, 0; 1, 1, 1; 0, 0, 0], [0, 0, 0; 0, 0, 0; 1, 1, 1]);
s = dif(u, grad_u);
d0_s = zeros(3, 9); div_s = repelem(eye(3), 1, 3);
f = - sum(dif(s, div_s), 2) * v + dif(p, grad_p);

UNV = MshEnt("D3F").UNV;
hF = MshEnt("D3F").area .^ (1/2);
%%
%[text] Domain: cube $\[0, 1\] \\times \[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured tetrahedral mesh + Alfeld splitting
%[text] Sub-element: tetrahedron $\[v\_1, v\_2, v\_3, c\]$, face 4 is a primal face, faces 1, 2, 3 are dual faces, edges 3, 5, 6 ($\[v\_i, c\]$) are dual edges.
msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], 4));
PrFace = find(ismember(msh.face.type, 0:6));
DlFace = find(ismember(msh.face.type, 1i));
BdFace = msh.bdEnt(2, 1:6);
DlEdge = find(ismember(msh.edge.type, 1i));
D3TElem = MshEnt("D3T").msh;
FcBar = [Fcn("D3R2", "1 - s - t"), Fcn("D3R2", "s"), Fcn("D3R2", "t")];
%%
%[text] Scalar $SDG\_1$ element
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $q|\_c$ at vertex 4 and at midpoints of dual edges \
%[text] Scalar $SDG\_1$ element space
%[text] - $P\_h = \\{ q \\in L^2 (\\Omega) : q|\_K \\in P\_1 (K), \\forall K \\in T\_h; \\ \[q\]\_F = 0, \\forall F \\in F^{dl}\_h \\}$ (pressure fixed on one dual edge) \
SDG_1S_FE = FE("D3T", "[1,x,y,z]", [NdDoF("D3", D3TElem, 0, [], d0_p, "EntIdx", 4), ...
    NdDoF("D3", D3TElem, 1, 1/2, d0_p, "EntIdx", [3, 5, 6])], "map", "affine");
SDG_1S_BC = BC(p, "node", nan, "edge", DlEdge(1));
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);
%[text] Vector $SDG\_1$ element
%[text] - Function space: $P^3\_1$
%[text] - Nodal DoF: each component at vertices of face 4 (shared), and at vertex 4 \
%[text] Vector $SDG\_1$ element space
%[text] - $U\_{h, 0} = \\{ \\mathbf{v} \\in L^2 (\\Omega, R^3) : \\mathbf{v}|\_K \\in P\_1 (K)^3, \\forall K \\in T\_h; \\ \[\\mathbf{v}\]\_F = 0, \\forall F \\in F^{pr, o}\_h; \\ \\mathbf{v}|\_{\\partial \\Omega} = 0 \\}$ \
SDG_1V_DoFs = NdDoF.empty;
for iComp = 1:3
    d0_ui = nan(3); d0_ui(:, iComp) = 0;
    SDG_1V_DoFs(end + 1) = NdDoF("D3", D3TElem, 2, [0, 1, 0; 0, 0, 1], d0_ui, "EntIdx", 4);
end
for iComp = 1:3
    d0_ui = nan(3); d0_ui(:, iComp) = 0;
    SDG_1V_DoFs(end + 1) = NdDoF("D3", D3TElem, 0, [], d0_ui, "EntIdx", 4, "share", false);
end
SDG_1V_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), SDG_1V_DoFs, "map", "affine");
SDG_1V_BC = BC(Fcn("D3", "0"), "node", nan, "face", BdFace);
SDG_1V_FES = FES(msh, SDG_1V_FE, SDG_1V_BC);
%[text] Matrix $SDG\_1$ element
%[text] - Function space: $P^{3 \\times 3}\_1$
%[text] - Moment DoF: $\\int\_F (\\tau \\mathbf{n})\_i \\, \\lambda\_k \\, ds$ on faces 1, 2, 3 (shared) and on face 4, $i = 1, 2, 3$, $\\lambda\_k$: barycentric coordinates of face
%[text] - Rows are mapped by the contravariant Piola transformation (`"map", "piolaDiv"`), which preserves these moments. \
%[text] Matrix $SDG\_1$ element space
%[text] - $\\Sigma\_h = \\{ \\tau \\in L^2 (\\Omega, R^{3 \\times 3}) : \\tau|\_K \\in P\_1 (K)^{3 \\times 3}, \\forall K \\in T\_h; \\ \[\\tau \\mathbf{n}\]\_F = 0, \\forall F \\in F^{dl}\_h \\}$ \
SDG_1M_DoFs = MoDoF.empty;
for EntIdx = {[1, 2, 3], 4}
    for iRow = 1:3
        d0_si = nan(3, 9); d0_si(:, iRow:3:9) = 0;
        SDG_1M_DoFs(end + 1) = MoDoF("D3", D3TElem, 2, FcBar, d0_si, "coef", UNV, "EntIdx", EntIdx{1}, "orien", true, ...
            "share", isequal(EntIdx{1}, [1, 2, 3]), "form", @(coef, fcn, tst) sum(fcn * coef(1)) .* tst, "GInt", GInt("D3F", 2));
    end
end
SDG_1M_FE = FE("D3T", FE.repFS("[1,x,y,z]", [3, 3]), SDG_1M_DoFs, "map", "piolaDiv");
SDG_1M_FES = FES(msh, SDG_1M_FE);
%%
%[text] SDG scheme: find $\\sigma\_h \\in \\Sigma\_h$, $\\mathbf{u}\_h \\in U\_{h, 0}$, and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}\n  \\sum\_{K} \\int\_K \\sigma\_h : \\tau\_h \\, dx + \\sum\_{K} \\int\_K \\mathbf{u}\_h \\cdot \\nabla \\cdot \\tau\_h \\, dx - \\sum\_{F \\in F\_h^{pr}} \\int\_F \\mathbf{u}\_h \\cdot \[ \\tau\_h \\mathbf{n}\] \\, ds = 0 &\\forall \\tau\_h \\in \\Sigma\_h \\\\\n  \\nu \\sum\_{K} \\int\_K \\sigma\_h : \\nabla \\mathbf{v}\_h \\, dx - \\nu \\sum\_{F \\in F\_h^{dl}} \\int\_F \\sigma\_h \\mathbf{n} \\cdot \[ \\mathbf{v}\_h \] \\, ds - \\sum\_{K} \\int\_K p\_h \\nabla \\cdot \\mathbf{v}\_h \\, dx + \\sum\_{F \\in F\_h^{dl}} \\int\_F p\_h \[ \\mathbf{v}\_h \\cdot \\mathbf{n} \] \\, ds = \\sum\_{K} \\int\_K \\mathbf{f} \\cdot \\mathbf{v}\_h \\, dx &\\forall \\mathbf{v}\_h \\in U\_{h, 0} \\\\\n  - \\sum\_{K} \\int\_K \\mathbf{u}\_h \\cdot \\nabla q\_h \\, dx + \\sum\_{F \\in F\_h^{pr}} \\int\_F \\mathbf{u}\_h \\cdot \\mathbf{n} \[ q\_h \] \\, ds = 0 &\\forall q\_h \\in P\_h\n\\end{cases}$
Sh = SDG_1M_FES; ord_Sh = 1;
Uh = SDG_1V_FES; ord_Uh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;

trls = [Sh, Uh, Ph]; tsts = [Sh, Uh, Ph];
is = 1; iu = 2; ip = 3; it = 1; iv = 2; iq = 3;

d0_t = d0_s; div_t = div_s; d0_v = d0_u; div_v = div_u; grad_v = grad_u; d0_q = d0_p; grad_q = grad_p;
Ast = DLF(msh, 3, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D3T", ord_Sh * 2));
But = [DLF(msh, 3, Fcn.cst(1), d0_u, div_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D3T", ord_Uh + ord_Sh)), ...
    DLF.interface(msh, 2, -UNV, d0_u, d0_t, "EntIdx", PrFace, "iTrl", iu, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D3F", ord_Uh + ord_Sh))];
Bsv = [DLF(msh, 3, Fcn.cst(v), d0_s, grad_v, "iTrl", is, "iTst", iv, "GInt", GInt("D3T", ord_Sh + ord_Uh)), ...
    DLF.interface(msh, 2, -UNV * v, d0_s, d0_v, "EntIdx", DlFace, "iTrl", is, "iTst", iv, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D3F", ord_Sh + ord_Uh))];
Bpv = [DLF(msh, 3, Fcn.cst(-1), d0_p, div_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D3T", ord_Ph + ord_Uh)), ...
    DLF.interface(msh, 2, UNV, d0_p, d0_v, "EntIdx", DlFace, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D3F", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 3, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D3T", ord_Uh + ord_Ph)), ...
    DLF.interface(msh, 2, UNV, d0_u, d0_q, "EntIdx", PrFace, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D3F", ord_Uh + ord_Ph))];
Fv = SLF(msh, 3, f, d0_v, "iTst", iv, "GInt", GInt("D3T", 2 + ord_Uh));
%%
%[text] Solution
%[text] DoFs of $\\sigma\_h$ are shared only inside macro elements (original tetrahedra), so $\\sigma\_h$ is eliminated macro element by macro element (static condensation).
[Stiff, Load] = assemble(msh, trls, tsts, [Ast, But, Bsv, Bpv, Buq], Fv);
[sh, uh, ph] = FEF.multi(trls, condSolve(Stiff, Load, trls, is));
%%
%[text] Error
%[text] $\\| \\mathbf{v} \\|\_{H^1, h}^2 = \\sum\_{K} \\| \\nabla \\mathbf{v} \\|\_K^2 + \\sum\_{F \\in F\_h^{dl}} h^{-1}\_F \\| \[\\mathbf{v}\] \\|\_F^2$
L2Norm_Sh = Norm(msh, 3, d0_s, "GInt", GInt("D3T", (1 + ord_Sh) * 2));
L2Norm_Uh = Norm(msh, 3, d0_u, "GInt", GInt("D3T", (1 + ord_Uh) * 2));
H1Norm_Uh = [Norm(msh, 3, grad_u, "GInt", GInt("D3T", ord_Uh * 2)), ...
    Norm(msh, 2, d0_u, "EntIdx", DlFace, "coef", hF \ 1, "form", @(coef, fcn, pow) coef .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D3F", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 3, d0_p, "GInt", GInt("D3T", (1 + ord_Ph) * 2));
fprintf('|s-sh|_L2: %e, |u-uh|_L2: %e, |u-uh|_H1: %e, |p-ph|_L2: %e\n', eNorm(msh, s, sh, L2Norm_Sh), eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, H1Norm_Uh), eNorm(msh, p, ph, L2Norm_Ph));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
