%[text] Nonlinear Poisson equation
%[text] $\\begin{cases}\n  - \\Delta u - \\alpha u^2 = f & \\text{in } \\Omega \\\\\n  u = g & \\text{on } \\Gamma \\\\\n\\end{cases}$
u = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y"); alph = 3;
d0_u = [0; 0]; dx_u = [1; 0]; dy_u = [0; 1]; grad_u = cat(3, dx_u, dy_u); dxx_u = [2; 0]; dyy_u = [0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u) - (u.^2) * alph;
g = u;
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$
%[text] Mesh: structured triangulation
msh = mshD2TS([0, 1, 0, 1], 16);
DirBd = [1, 2, 3, 4, -1, -2, -3, -4];
DirNode = find(ismember(msh.node.type, DirBd));
%%
%[text] $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
%[text] Continuous $P\_1$ element space
%[text] - $U\_{h, g} = \\{ v \\in H^1 (\\Omega) : v|\_{K} \\in P\_1 (K), \\forall K \\in T\_h; \\ v|\_{\\Gamma\_D} = g \\}$ \
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], d0_u));
P1_BC = BC(g, "node", DirNode);
P1_FES = FES(msh, P1_FE, P1_BC);
%%
%[text] Iterative FEM scheme 1: find $u\_h \\in U\_{h, g}$ such that
%[text] $\\int\_{\\Omega} \\nabla u^{n+1}\_h \\cdot \\nabla v\_h \\, dx \\, dy - \\alpha \\int\_{\\Omega} u^n\_h u^{n+1}\_h v\_h \\, dx \\, dy = \\int\_{\\Omega} f v\_h \\, dx \\, dy , \\quad \\forall v\_h \\in U\_{h, 0}$
%[text] Iterative FEM scheme 2: find $u\_h \\in U\_{h, g}$ such that
%[text] $\\int\_{\\Omega} \\nabla u^{n+1}\_h \\cdot \\nabla v\_h \\, dx \\, dy = \\int\_{\\Omega} f v\_h \\, dx \\, dy + \\alpha \\int\_{\\Omega} (u^n\_h)^2 v\_h \\, dx \\, dy, \\quad \\forall v\_h \\in U\_{h, 0}$
Uh = P1_FES; ord = 1;

d0_v = d0_u; grad_v = grad_u;
Auv = DLF(msh, 2, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D2T", (ord - 1) * 2));
Awuv = LDLF(msh, 2, Fcn.cst(-alph), d0_u, d0_u, d0_v, "GInt", GInt("D2T", ord * 3));
Fv = SLF(msh, 2, f, d0_v, "GInt", GInt("D2T", 1 + ord));
Fwv = LSLF(msh, 2, Fcn.cst(alph), d0_u, d0_v, "form", @(load, pre, tst) load .* (pre).^2 .* tst, "GInt", GInt("D2T", ord * 3));
%%
%[text] Scheme 1
% L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
% uh0 = FEF(Uh, zeros(Uh.nGlDoF, 1));
% 
% while true
%     [Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv, "preSol", uh0, "Awuv", Awuv);
%     uh = FEF(Uh, Stiff \ Load);
%     if eNorm(msh, Fcn.cst(0), subt(uh, uh0), L2Norm) < 1e-4
%         break
%     else
%         uh0 = uh;
%     end
% end
%[text] Scheme 2
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
uh0 = FEF(Uh, zeros(Uh.nGlDoF, 1));

while true
    [Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv, "preSol", uh0, "Fwv", Fwv);
    uh = FEF(Uh, Stiff \ Load);
    if eNorm(msh, Fcn.cst(0), subt(uh, uh0), L2Norm) < 1e-4
        break
    else
        uh0 = uh;
    end
end
%%
%[text] Error
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
H1Norm = Norm(msh, 2, grad_u, "GInt", GInt("D2T", (ord - 1) * 2));
fprintf('|u-uh|_L2: %e, |u-uh|_H1: %e\n', eNorm(msh, u, uh, L2Norm), eNorm(msh, u, uh, H1Norm));
%%
%[text] Plot
plotFEF(uh);
plotFcn(msh, u);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
