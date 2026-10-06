%[text] # DoF - Degree of Freedom
%%
%[text] ## Nodal DoF
%%
%[text] ### Nodal DoF for scalar function at element vertex
%[text] - $v|\_c$ at vertex \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]);
disp(eval(ndDoF, Fcn("D2", "x+y")));
%[text] - $v|\_c$ at second vertex \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0], "EntIdx", 2);
disp(eval(ndDoF, Fcn("D2", "x+y")));
%[text] - $\\partial\_x v |\_c$ and $\\partial\_y v |\_c$ at vertex \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [1; 0]);
disp(eval(ndDoF, Fcn("D2", "x^2 + x*y + y^2")));

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 1]);
disp(eval(ndDoF, Fcn("D2", "x^2 + x*y + y^2")));
%%
%[text] ### Nodal DoF for vector function at vertex
%[text] - $v\_x \\, |\_c$ and $v\_y \\, |\_c$ at vertex \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]);
disp(eval(ndDoF, Fcn("D2", "[x+y; x-y]")));

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [nan, 0; nan, 0]);
disp(eval(ndDoF, Fcn("D2", "[x+y; x-y]")));
%%
%[text] ### Nodal DoF for scalar function on edge
%[text] - $v|\_c$ at midpoint of edge \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0; 0]);
disp(eval(ndDoF, Fcn("D2", "x")));
%[text] - $v|\_c$ at vertex of edge \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0; 0]);
disp(eval(ndDoF, Fcn("D2", "x")));
%%
%[text] ### Nodal DoF for vector function on edge
%[text] - $v \\cdot n |\_c$ at midpoint of edge \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, 1/2, [0, 0; 0, 0], "coef", MshEnt("D2L").UNV);
disp(eval(ndDoF, Fcn("D2", "[x;y]")));
%[text] - $v \\cdot n |\_c$ at vertex of edge \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, [0, 1], [0, 0; 0, 0], "coef", MshEnt("D2L").UNV);
disp(eval(ndDoF, Fcn("D2", "[x;y]")));
%%
%[text] ### Nodal DoF for scalar function in element 
%[text] - $v|\_c$ at centroid of element \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [1/3; 1/3], [0; 0]);
disp(eval(ndDoF, Fcn("D2", "x")));
%[text] - $v|\_c$ at vertex of element \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [0, 0; 1, 0; 0, 1]', [0; 0]);
disp(eval(ndDoF, Fcn("D2", "x")));
%[text] - $v|\_c$ at edge-midpoint of element \
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [1/2, 0; 1/2, 1/2; 0, 1/2]', [0; 0]);
disp(eval(ndDoF, Fcn("D2", "x")));
%%
%[text] ### Nodal DoF for trace function
%[text] - $\\mu |\_c$ at vertex of edge \
ndDoF = NdDoF("D2R1", MshEnt("D2T").msh, 1, [0, 1], 0);
eval(ndDoF, Fcn("D2", "x"))
%[text] - $\\hat{\\mu} |\_c$ at vertex of reference edge \
ndDoF = NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0);
disp(eval(ndDoF, Fcn("D2R1", "s")));
%%
%[text] ## Moment DoF
%%
%[text] ### Moment DoF for scalar function on edge
%[text] - $\\int\_e v \\, ds$ and $\\int\_e v (s-1/2) \\, ds$ on edge \
moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), [0; 0]);
disp(eval(moDoF, Fcn("D2", "1")));

moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0; 0]);
disp(eval(moDoF, Fcn("D2", "1")));
%%
%[text] ### Moment DoF for vector function on edge
%[text] - $\\int\_e v \\cdot n \\, ds$ and $\\int\_e v \\cdot n (s - 1/2) \\, ds$ on edge \
moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV);
disp(eval(moDoF, Fcn("D2", "[1;1]")));

moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV);
disp(eval(moDoF, Fcn("D2", "[1;1]")));
%%
%[text] ### Moment DoF for scalar function on element
%[text] - $\\int\_k v \\, dx \\, dy$ on element \
moDoF = MoDoF("D2", MshEnt("D2T").msh, 2, Fcn("D2R", "1"), [0; 0]);
disp(eval(moDoF, Fcn("D2", "1")));
%%
%[text] ### Moment DoF for trace function
%[text] - $\\int\_e \\mu \\, ds$ and $\\int\_e \\mu (s-1/2) \\, ds$ on edge \
moDoF = MoDoF("D2R1", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), 0);
disp(eval(moDoF, Fcn("D2", "1")));

moDoF = MoDoF("D2R1", MshEnt("D2T").msh, 1, Fcn("D2R1", "s-1/2"), 0);
disp(eval(moDoF, Fcn("D2", "1")));
%[text] - $\\int\_{\\hat{e}} \\hat{\\mu} \\, d \\hat{s}$ and $\\int\_\\hat{e} \\hat{\\mu} (\\hat{s}-1/2) \\, d \\hat{s}$ \
moDoF = MoDoF("D2R1", MshEnt("D2LR").msh, 1, Fcn("D2R1", "1"), 0);
disp(eval(moDoF, Fcn("D2R1", "1")));

moDoF = MoDoF("D2R1", MshEnt("D2LR").msh, 1, Fcn("D2R1", "s-1/2"), 0);
disp(eval(moDoF, Fcn("D2R1", "1")));
%%
%[text] ## DoF.Method
%[text] Example: nodal DoF
%[text] - $v|\_c$ at vertex
%[text] - $v|\_c$ at 1/3-points of edge
%[text] - $v|\_c$ at 1/3-points of element \
dof = [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 1, [1/3, 2/3], [0; 0]), ...
    NdDoF("D2", MshEnt("D2T").msh, 2, [1/2, 1/4; 1/4, 1/2; 1/4, 1/4]', [0; 0])];
%[text] ### DoF.`cumDoF`
disp(dof.cumDoF); % expect: 12
disp(dof.cumDoF(1)); % expect: 3
disp(dof.cumDoF(2)); % expect: 9
disp(dof.cumDoF(3)); % expect: 12
%[text] ### DoF.`sub2ind`
disp(sub2ind(dof, 1)); % expect: 1:3 (group 1)
disp(sub2ind(dof, 2)); % expect: 4:9 (group 2)
disp(sub2ind(dof, 3)); % expect: 10:12 (group 3)

disp(sub2ind(dof, 1, 1:3, 1)); % expect: 1:3 (group 1, entity 1:3, sample 1)
disp(sub2ind(dof, 2, 1:3, 1)); % expect: 4:6 (group 2, entity 1:3, sample 1)
disp(sub2ind(dof, 2, 1:3, 2)); % expect: 7:9 (group 2, entity 1:3, sample 2)
disp(sub2ind(dof, 3, 1, 1:3)); % expect: 10:12 (group 3, entity 1, sample 1:3)
%[text] ### DoF.`ind2sub`
[iDoF, iEnt, iSamp] = ind2sub(dof, 2);
disp([iDoF, iEnt, iSamp]); % expect: 1, 2, 1
[iDoF, iEnt, iSamp] = ind2sub(dof, 5);
disp([iDoF, iEnt, iSamp]); % expect: 2, 2, 1
[iDoF, iEnt, iSamp] = ind2sub(dof, 7);
disp([iDoF, iEnt, iSamp]); % expect: 2, 1, 2
[iDoF, iEnt, iSamp] = ind2sub(dof, 11);
disp([iDoF, iEnt, iSamp]); % expect: 3, 1, 2

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
