%[text] Stokes equation:
%[text] $\\begin{cases}\n- \\Delta u + \\nabla p = f &\\text{in } \\Omega, \\\\\n\\nabla \\cdot u = 0 &\\text{in } \\Omega, \\\\\nu = 0 &\\text{on } \\partial \\Omega.\n\\end{cases}$
%[text] Mixed formulation:
%[text] $\\begin{cases}\n\\sigma + \\nabla u = 0 &\\text{in } \\Omega, \\\\\n\\nabla \\cdot \\sigma + \\nabla p = f &\\text{in } \\Omega, \\\\\n\\nabla \\cdot u = 0 &\\text{in } \\Omega, \\\\\nu = 0 &\\text{on } \\partial \\Omega,\n\\end{cases}$
u = Fcn("D2", "[sin(pi*x)^2*sin(2*pi*y); - sin(2*pi*x)*sin(pi*y)^2]");
p = Fcn("D2", "10*(x-1/2)*(y-1/2)");
d0_u = [0, 0; 0, 0]'; dxdy_u = [1, 0; 0, 1]; dxdx_u = [1, 0; 1, 0]'; dydy_u = [0, 1; 0, 1]'; grad_u = cat(3, dxdx_u, dydy_u);
s = - dif(u, grad_u);
d0_s = [0, 0; 0, 0; 0, 0; 0, 0]'; dxdxdydy_s = [1, 0; 1, 0; 0, 1; 0 ,1]';
d0_p = [0, 0]'; dx_p = [1, 0]'; dy_p = [0, 1]'; grad_p = cat(3, dx_p, dy_p);
f = sum(dif(s, dxdxdydy_s), 2) + dif(p, grad_p);

UNV = MshEnt("D2L").UNV; UTV = MshEnt("D2L").UTV;
len = MshEnt("D2L").len;
%%
%[text] Domain: rectangular $\[0, 1\] \\times \[0, 1\]$.
%[text] Dirichlet boundary.
%[text] Mesh: structured triangulation.
msh = mshSplit(mshD2TS([0, 1, 0, 1], 32));
% load("mshD2T_10P.mat"); msh = mshSplit(msh);
PrEdge = find(ismember(msh.edge.type, [0, 1, 2, 3, 4]));
DlEdge = find(ismember(msh.edge.type, 1i));
%%
%[text] SDG P1 element.
D2TElem = MshEnt("D2T").msh;
SDG_1S_FE = FE("D2T", "[1,x,y]", [NdDoF("D2", D2TElem, 1, [0, 1], d0_p, "EntIdx", 1), ...
    NdDoF("D2", D2TElem, 0, [], d0_p, "EntIdx", 3, "share", false)]);
SDG_1S_BC = BC(p, "node", nan, "edge", 1);
SDG_1S_FES = FES(msh, SDG_1S_FE, SDG_1S_BC);

SDG_1V_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", D2TElem, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true), ...
    NdDoF("D2", D2TElem, 1, [0, 1], d0_u, "coef", UNV, "EntIdx", 1, "orien", true, "share", false)]);
SDG_1V_FES = FES(msh, SDG_1V_FE);

load("SDG_1M_FE.mat");
SDG_1M_FES = FES(msh, SDG_1M_FE);
%%
%[text] SDG P0 element.
SDG_0S_FE = FE("D2T", "1", NdDoF("D2", D2TElem, 1, 1/2, d0_p, "EntIdx", 1));
SDG_0S_BC = BC(p, "edge", 1);
SDG_0S_FES = FES(msh, SDG_0S_FE, SDG_0S_BC);

SDG_0V_FE = FE("D2T", FE.repFS("1", [2, 1]), ...
    NdDoF("D2", D2TElem, 1, 1/2, d0_u, "coef", UNV, "EntIdx", [2, 3], "orien", true));
SDG_0V_FES = FES(msh, SDG_0V_FE);

SDG_0M_FE = FE("D2T", FE.repFS("1", [2, 2]), ...
    [NdDoF("D2", D2TElem, 1, 1/2, d0_s, "coef", [UNV, UNV], "EntIdx", 1, "form", @(coef, fcn) dot(fcn * coef(1), coef(2))), ...
    NdDoF("D2", D2TElem, 1, 1/2, d0_s, "coef", [UNV, UTV], "form", @(coef, fcn) dot(fcn * coef(1), coef(2)))]);
SDG_0M_FES = FES(msh, SDG_0M_FE);
%%
%[text] SDG scheme: find $\\sigma\_h \\in \\Sigma\_h$, $u\_h \\in U\_h$, and $p\_h \\in P\_h$ such that
%[text] $\\begin{cases}\n\\sum\_{K \\in T\_h} \\int\_K \\sigma\_h \\tau\_h \\, dx - \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla \\cdot \\tau\_h \\, dx + \\sum\_{e \\in \\Gamma\_h^{dl}} \\int\_e u\_h \\cdot n \[ \\tau\_h n \\cdot n \] \\, ds = 0 &\\forall \\tau\_h \\in \\Sigma\_h, \\\\\n- \\sum\_{K \\in T\_h} \\int\_K \\sigma\_h \\nabla v\_h \\, dx + \\sum\_{e \\in \\Gamma\_h^{dl}} \\int\_e \\sigma\_h n \\cdot t \[ v\_h \\cdot t\] \\, ds + \\sum\_{e \\in \\Gamma\_h^{pr}} \\int\_e \\sigma\_h n \[ v\_h \] \\, ds - \\sum\_{K \\in T\_h} \\int\_K p\_h \\nabla \\cdot v\_h \\, dx + \\sum\_{e \\in \\Gamma\_h^{pr}} \\int\_e p\_h \[ v\_h \\cdot n \] \\, ds = \\sum\_{K \\in T\_h} \\int\_K f v\_h \\, dx &\\forall v\_h \\in U\_h, \\\\\n- \\sum\_{K \\in T\_h} \\int\_K u\_h \\nabla q\_h \\, dx + \\sum\_{e \\in \\Gamma\_h^{dl}} \\int\_e u\_h \\cdot n \[ q\_h \] \\, ds = 0 &\\forall q\_h \\in Q\_h.\n\\end{cases}$
Sh = SDG_1M_FES; ord_Sh = 1;
Uh = SDG_1V_FES; ord_Uh = 1;
Ph = SDG_1S_FES; ord_Ph = 1;
% Sh = SDG_0M_FES; ord_Sh = 1;
% Uh = SDG_0V_FES; ord_Uh = 1;
% Ph = SDG_0S_FES; ord_Ph = 1;

trls = [Sh, Uh, Ph]; tsts = [Sh, Uh, Ph];
is = 1; iu = 2; ip = 3; it = 1; iv = 2; iq = 3;

d0_t = d0_s; dxdxdydy_t = dxdxdydy_s; d0_v = d0_u; dxdy_v = dxdy_u; grad_v = grad_u; d0_q = d0_p; grad_q = grad_p;
Ast = DLF(msh, 2, Fcn.cst(1), d0_s, d0_t, "iTrl", is, "iTst", it, "GInt", GInt("D2T", ord_Sh * 2));
But = [DLF(msh, 2, Fcn.cst(-1), d0_u, dxdxdydy_t, "iTrl", iu, "iTst", it, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D2T", ord_Uh + ord_Sh - 1)), ...
    DLF.interface(msh, 1, [UNV, UNV, UNV], d0_u, d0_t, "EntIdx", DlEdge, "iTrl", iu, "iTst", it, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl, coef(1)) .* dot(tst * coef(2), coef(3)), "GInt", GInt("D2L", ord_Uh + ord_Sh))];
Bsv = [DLF(msh, 2, Fcn.cst(-1), d0_s, grad_v, "iTrl", is, "iTst", iv, "GInt", GInt("D2T", ord_Sh + ord_Uh - 1)), ...
    DLF.interface(msh, 1, [UNV, UTV, UTV], d0_s, d0_v, "EntIdx", DlEdge, "iTrl", is, "iTst", iv, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef(1), coef(2)) .* dot(tst, coef(3)), "GInt", GInt("D2L", ord_Sh + ord_Uh)), ...
    DLF.interface(msh, 1, UNV, d0_s, d0_v, "EntIdx", PrEdge, "iTrl", is, "iTst", iv, "tstOpr", "jump", "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D2L", ord_Sh + ord_Uh))];
Bpv = [DLF(msh, 2, Fcn.cst(-1), d0_p, dxdy_v, "iTrl", ip, "iTst", iv, "GInt", GInt("D2T", ord_Ph + ord_Uh - 1)), ...
    DLF.interface(msh, 1, UNV, d0_p, d0_v, "EntIdx", PrEdge, "iTrl", ip, "iTst", iv, "tstOpr", "jump", "GInt", GInt("D2L", ord_Ph + ord_Uh))];
Buq = [DLF(msh, 2, Fcn.cst(-1), d0_u, grad_q, "iTrl", iu, "iTst", iq, "GInt", GInt("D2T", ord_Uh + ord_Ph - 1)), ...
    DLF.interface(msh, 1, UNV, d0_u, d0_q, "EntIdx", DlEdge, "iTrl", iu, "iTst", iq, "tstOpr", "jump", "GInt", GInt("D2L", ord_Uh + ord_Ph))];
Fv = SLF(msh, 2, f, d0_v, "iTst", iv, "GInt", GInt("D2T", 2 + ord_Uh));
%%
%[text] Solution.
[Stiff, Load] = assemble(msh, trls, tsts, [Ast, But, Bsv, Bpv, Buq], Fv);
[sh, uh, ph] = FEF.multi(trls, Stiff \ Load);
%%
%[text] Error.
%[text] $\\| \\sigma\_h \\|\_{div, h}^2 = \\sum\_{k \\in T\_h} \\| \\nabla \\cdot \\sigma\_h \\|\_K^2 + \\sum\_{e \\in \\Gamma\_h^{dl}} h^{-1}\_e \\| \[\\sigma\_h n \\cdot n \] \\|\_e^2, \\\\\n\\| u\_h \\|\_{1, h}^2 = \\sum\_{K \\in T\_h} \\| \\nabla u\_h \\|\_K^2 + \\sum\_{e \\in \\Gamma\_h^{pr}} h^{-1}\_e \\| \[ u\_h \] \\|\_e^2 + \\sum\_{e \\in \\Gamma\_h^{dl}} h^{-1}\_e \\| \[ u\_h \\cdot t \] \\|\_e^2.$
L2Norm_Sh = Norm(msh, 2, d0_s, "GInt", GInt("D2T", (1 + ord_Sh) * 2));
divNorm_Sh = [Norm(msh, 2, dxdxdydy_s, "form", @(coef, fcn, pow) sum(abs(sum(fcn, 2)) .^ pow), "GInt", GInt("D2T", ord_Sh * 2)), ...
    Norm(msh, 1, d0_s, "EntIdx", DlEdge, "coef", [len \ 1, UNV, UNV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(fcn * coef(2), coef(3))) .^ pow, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Sh * 2))];
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (1 + ord_Uh) * 2));
H1Norm_Uh = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", PrEdge, "coef", len \ 1, "form", @(coef, fcn, pow) coef(1) .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", DlEdge, "coef", [len \ 1, UTV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(fcn, coef(2))) .^ pow, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (1 + ord_Ph) * 2));
fprintf('|s-s_h|_L2: %e, |s-s_h|_div: %e, |u-u_h|_L2: %e, |u-u_h|_H1: %e, |p-p_h|_L2: %e\n', eNorm(msh, s, sh, L2Norm_Sh), eNorm(msh, s, sh, divNorm_Sh), eNorm(msh, u, uh, L2Norm_Uh), eNorm(msh, u, uh, H1Norm_Uh), eNorm(msh, p, ph, L2Norm_Ph));
%%
Ps = SDG_0M_FES.proj(s);
Pu = SDG_0V_FES.proj(u);
Pp = SDG_0S_FES.proj(p);
L2Norm_Sh = Norm(msh, 2, d0_s, "GInt", GInt("D2T", (1 + ord_Sh) * 2));
divNorm_Sh = [Norm(msh, 2, dxdxdydy_s, "form", @(coef, fcn, pow) sum(abs(sum(fcn, 2)) .^ pow), "GInt", GInt("D2T", ord_Sh * 2)), ...
    Norm(msh, 1, d0_s, "EntIdx", DlEdge, "coef", [len \ 1, UNV, UNV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(fcn * coef(2), coef(3))) .^ pow, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Sh * 2))];
L2Norm_Uh = Norm(msh, 2, d0_u, "GInt", GInt("D2T", (1 + ord_Uh) * 2));
H1Norm_Uh = [Norm(msh, 2, grad_u, "GInt", GInt("D2T", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", PrEdge, "coef", len \ 1, "form", @(coef, fcn, pow) coef(1) .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2)), ...
    Norm(msh, 1, d0_u, "EntIdx", DlEdge, "coef", [len \ 1, UTV], "form", @(coef, fcn, pow) coef(1) .* abs(dot(fcn, coef(2))) .^ pow, "fcnOpr", "jump", "GInt", GInt("D2L", ord_Uh * 2))];
L2Norm_Ph = Norm(msh, 2, d0_p, "GInt", GInt("D2T", (1 + ord_Ph) * 2));
fprintf('|s-s_h|_L2: %e, |s-s_h|_div: %e, |u-u_h|_L2: %e, |u-u_h|_H1: %e, |p-p_h|_L2: %e\n', eNorm(msh, Fcn.cst(0), subt(Ps, sh), L2Norm_Sh), eNorm(msh, Fcn.cst(0), subt(Ps, sh), divNorm_Sh), eNorm(msh, Fcn.cst(0), subt(Pu, uh), L2Norm_Uh), eNorm(msh, Fcn.cst(0), subt(Pu, uh), H1Norm_Uh), eNorm(msh, Fcn.cst(0), subt(Pp, ph), L2Norm_Ph));
%%
%[text] Plot.
plotFEF(sh);
plotFcn(msh, s);
%%
plotFEF(uh);
plotFcn(msh, u);
%%
plotFEF(ph);
plotFcn(msh, p);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
