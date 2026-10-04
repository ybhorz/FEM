%[text] # FE - Finite Element
%%
%[text] ## $P\_0$ element
%[text] - Domain: triangle
%[text] - Function space: $P\_0$
%[text] - Nodal DoF: $v|\_c$ at centroid \
P0_FE = FE("D2T", "1", NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0]));
disp(P0_FE.base.fun)
%%
%[text] ## $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at element vertex \
P1_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]));
bases = P1_FE.base;

base = bases(1);
disp(base.fun);
disp([base.eval("[x1;y1]"), base.eval("[x2;y2]"), base.eval("[x3;y3]")]);
%%
%[text] ## $P\_2$ element
%[text] - Element: triangle
%[text] - Function space: $P\_2$
%[text] - Nodal DoF: $v|\_c$ at element vertex and edge midpoint \
P2_FE = FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0])]);
bases = P2_FE.base;

base = bases(1);
disp(base.fun);
disp([base.eval("[x1;y1]"), base.eval("[x2;y2]"), base.eval("[x3;y3]"), ...
    base.eval("[(x1+x2)/2;(y1+y2)/2]"), base.eval("[(x2+x3)/2;(y2+y3)/2]"), base.eval("[(x3+x1)/2;(y3+y1)/2]")]);
%%
%[text] ## CR element
%[text] - Element: triangle
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $v|\_c$ at edge midpoint \
CR_FE = FE("D2T", "[1,x,y]", NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0]));
bases = CR_FE.base;

base = bases(1);
disp(base.fun);
disp([base.eval("[(x1+x2)/2;(y1+y2)/2]"), base.eval("[(x2+x3)/2;(y2+y3)/2]"), base.eval("[(x3+x1)/2;(y3+y1)/2]")]);
%%
%[text] ## $RT\_0$ element
%[text] - Element: triangle
%[text] - Function space: $P^2\_0 + x P\_0$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ on edge \
RT0_FE = FE("D2T", "[1,0; 0,1; x,y].'", ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true));
bases = RT0_FE.base;

base = bases(1);
disp(base.fun);
disp([int(dot(base, MshEnt("D2T").UNV(1)), "D2L", 1).fun, ...
    int(dot(base, MshEnt("D2T").UNV(2)), "D2L", 2).fun, ...
    int(dot(base, MshEnt("D2T").UNV(3)), "D2L", 3).fun]);
%%
%[text] ## $BDM\_1$ element
%[text] - Domain: triangle
%[text] - Function space: $P^2\_1$
%[text] - Moment DoF: $\\int\_e v \\cdot n \\, ds$ and $\\int\_e (v \\cdot n) (s - 1/2) \\, ds$ on edge \
BDM1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [MoDoF("D2", MshEnt("D2T").msh, 1, Fcn.cst(1), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true), ...
    MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "orien", true)]);
bases = BDM1_FE.base;

base = bases(1);
disp(base.fun);

disp([int(dot(base, MshEnt("D2T").UNV(1)), "D2L", 1).fun, ...
    int(dot(base, MshEnt("D2T").UNV(2)), "D2L", 2).fun, ...
    int(dot(base, MshEnt("D2T").UNV(3)), "D2L", 3).fun, ...
    int(dot(base, MshEnt("D2T").UNV(1)) .* Fcn("D2T", "(x-x1+y-y1)/(x2-x1+y2-y1)-1/2"), "D2L", 1).fun, ...
    int(dot(base, MshEnt("D2T").UNV(2)) .* Fcn("D2T", "(x-x2+y-y2)/(x3-x2+y3-y2)-1/2"), "D2L", 2).fun, ...
    int(dot(base, MshEnt("D2T").UNV(3)) .* Fcn("D2T", "(x-x3+y-y3)/(x1-x3+y1-y3)-1/2"), "D2L", 3).fun]);
%%
%[text] ## Vector $P\_1$ element
%[text] - Element: triangle
%[text] - Function space: $P\_1^2$
%[text] - Nodal DoF: $v\_x \\, |\_c$ and $v\_y \\, |\_c$ at element vertex \
VP1_FE = FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
    [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]), ...
    NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0])]);
bases = VP1_FE.base;

base = bases(1);
disp(base.fun);
disp([base.eval("[x1;y1]"), base.eval("[x2;y2]"), base.eval("[x3;y3]")]);
%%
%[text] ## Trace P1 element
%[text] - Element: reference edge \[0, 1\]
%[text] - Function space: $P\_1$
%[text] - Nodal DoF: $\\hat{\\mu} |\_c$ at edge vertex \
TP1_FE = FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0));
bases = TP1_FE.base;

base = bases(1);
disp(base.fun)
disp([base.eval(0), base.eval(1)]);
%%
%[text] ## FE.Method
%[text] ### FE.repFS
disp(FE.repFS("[1,x,y]", [2, 1]))
disp(FE.repFS("[1,x,y]", [2, 2]))

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
