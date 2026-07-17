%[text] Nonlinear PDE:
%[text] $- \\Delta u - \\alpha u^2 = f$
u = Fcn("D2", "sin(pi*x)*sin(pi*y)+x*y"); alph = 3;
d0_u = [0; 0]; dx_u = [1; 0]; dy_u = [0; 1]; grad_u = cat(3, dx_u, dy_u); dxx_u = [2; 0]; dyy_u = [0; 2];
f = - dif(u, dxx_u) - dif(u, dyy_u) - (u.^2) * alph;
g_D = u;
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$.
%[text] Dirichlet boundary.
%[text] Mesh: structured triangulation.
msh = mshD2TS([0, 1, 0, 1], 16);
DirBd = [1, 2, 3, 4, -1, -2, -3, -4];
DirNode = find(ismember(msh.node.type, DirBd));
%%
%[text] P1 Lagrange element.
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
P1_BC = BC(g_D, "node", DirNode);
P1_FES = FES(msh, P1_FE, P1_BC);
%%
%[text] Iterative FEM scheme 1: find $u\_h \\in U\_h$ such that
%[text] $\\int\_{\\Omega} \\nabla u^{n+1}\_h \\cdot \\nabla v\_h \\, dx - \\alpha \\int\_{\\Omega} u^n\_h u^{n+1}\_h v\_h = \\int\_{\\Omega} f v\_h \\, \ndx , \\quad \\forall v\_h \\in U\_h.$
%[text] Iterative FEM scheme 2: find $u\_h \\in U\_h$ such that
%[text] $\\int\_{\\Omega} \\nabla u^{n+1}\_h \\cdot \\nabla v\_h \\, dx = \\int\_{\\Omega} f v\_h + \\alpha \\int\_{\\Omega} (u^n\_h)^2 v\_h \\, \ndx , \\quad \\forall v\_h \\in U\_h.$
Uh = P1_FES; ord = 1;

d0_v = d0_u; grad_v = grad_u;
Auv = DLF(msh, 2, Fcn.cst(1), grad_u, grad_v, "GInt", GInt("D2T", (ord - 1) * 2));
Awuv = LDLF(msh, 2, Fcn.cst(-alph), d0_u, d0_u, d0_v, "GInt", GInt("D2T", ord * 2));
Fv = SLF(msh, 2, f, d0_v, "GInt", GInt("D2T", 1 + ord));
Fwv = LSLF(msh, 2, Fcn.cst(alph), d0_u, d0_v, "form", @(load, pre, tst) load .* (pre).^2 .* tst, "GInt", GInt("D2T", ord *3));
%%
%[text] Iteration solution.
%[text] Scheme 1.
L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
uh0 = FEF(Uh, zeros(Uh.nGlDoF, 1));

while true
    [Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv, "preSol", uh0, "Awuv", Awuv);
    uh = FEF(Uh, Stiff \ Load);
    plotFEF(uh);
    if eNorm(msh, Fcn.cst(0), subt(uh, uh0), L2Norm) < 1e-3
        break
    else
        uh0 = uh;
    end
end
%[text] Scheme 2.
% L2Norm = Norm(msh, 2, d0_u, "GInt", GInt("D2T", ord * 2));
% uh0 = FEF(Uh, zeros(Uh.nGlDoF, 1));
% 
% while true
%     [Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv, "preSol", uh0, "Fwv", Fwv);
%     uh = FEF(Uh, Stiff \ Load);
%     % plotFEF(uh);
%     if eNorm(msh, Fcn.cst(0), subt(uh, uh0), L2Norm) < 1e-3
%         break
%     else
%         uh0 = uh;
%     end
% end
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
%   data: {"layout":"onright","rightPanelPercent":40}
%---
