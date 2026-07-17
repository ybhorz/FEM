%[text] ## NdDoF
%[text] domn = D2, EntDim = 0
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]);
eval(ndDoF, Fcn("D2", "1"))
eval(ndDoF, Fcn("D2", "x"))
eval(ndDoF, Fcn("D2", "y"))

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [1; 0]);
eval(ndDoF, Fcn("D2", "x^2"))
%%
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, 0; 0, 0]);
eval(ndDoF, Fcn("D2", "[x;y]"))
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [0, nan; 0, nan]);
eval(ndDoF, Fcn("D2", "[x;y]"))
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 0, [], [1, 0; 0, 1]);
eval(ndDoF, Fcn("D2", "[x^2;y^2]"))
%%
%[text] domn = D2, EntDim = 1
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, [1/3, 2/3], [0; 0]);
eval(ndDoF, Fcn("D2", "1"))
eval(ndDoF, Fcn("D2", "x"))
eval(ndDoF, Fcn("D2", "y"))

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, [1/3, 2/3], [0, 0; 0, 0], "coef", MshEnt("D2L").UNV);
eval(ndDoF, Fcn("D2", "[x;0]"))

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 1, [1/3, 2/3], [0; 0]);
eval(ndDoF, dot(Fcn("D2", "[x;0]"), MshEnt("D2L").UNV))
%%
%[text] domn = D2, EntDim = 2
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2; 1/3, 1/3]', [0; 0]);
eval(ndDoF, Fcn("D2", "1"))
eval(ndDoF, Fcn("D2", "x"))
eval(ndDoF, Fcn("D2", "y"))

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [1/3, 2/3; 2/3, 1/3]', [0, 0; 0, 0], "coef", MshEnt("D2T").UNV(2));
eval(ndDoF, Fcn("D2", "[x;0]"))

ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [1/3, 2/3; 2/3, 1/3]', [0; 0]);
eval(ndDoF, dot(Fcn("D2", "[x;0]"), MshEnt("D2T").UNV(2)))
%%
ndDoF = NdDoF("D2", MshEnt("D2T").msh, 2, [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2; 1/3, 1/3]', [0; 0]);
ndDoF.coord
ndDoF.loc
%%
%[text] domn = D2R1, msh.type = D1
ndDoF = NdDoF("D2R1", MshEnt("D2LR").msh, 1, [1/3, 2/3], 0);
eval(ndDoF, Fcn("D2R1", "1"))
eval(ndDoF, Fcn("D2R1", "s"))
%%
%[text] domn = D2R1, msh.type = D2T
ndDoF = NdDoF("D2R1", MshEnt("D2T").msh, 1, [1/3, 2/3], 0);
eval(ndDoF, Fcn("D2", "1"))
eval(ndDoF, Fcn("D2", "x"))

ndDoF = NdDoF("D2R1", MshEnt("D2T").msh, 1, [1/3, 2/3], 0);
eval(ndDoF, dot(MshEnt("D2L").UNV, Fcn("D2", "[x;0]")))
%%
%[text] ## MoDoF
%[text] domn = D2, EntDim = 1
moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), [0; 0]);
eval(moDoF, Fcn("D2", "1"))

moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV);
eval(moDoF, Fcn("D2", "[1;1]"))

moDoF = MoDoF("D2", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), [0; 0]);
eval(moDoF, dot(Fcn("D2", "[1;1]"), MshEnt("D2L").UNV))
%%
%[text] domn = D2, EntDim = 2
moDoF = MoDoF("D2", MshEnt("D2T").msh, 2, Fcn("D2", "1"), [0; 0]);
eval(moDoF, Fcn("D2", "1"))
%%
%[text] domn = D2R1, msh.type = D2R1
moDoF = MoDoF("D2R1", MshEnt("D2LR").msh, 1, Fcn("D2R1", "1"), 0);
eval(moDoF, Fcn("D2R1", "1"))
eval(moDoF, Fcn("D2R1", "s"))
%%
%[text] domn = D2R1, msh.type = D2T
moDoF = MoDoF("D2R1", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), 0);
eval(moDoF, Fcn("D2", "1"))

moDoF = MoDoF("D2R1", MshEnt("D2T").msh, 1, Fcn("D2R1", "1"), 0);
eval(moDoF, dot(Fcn("D2", "[1;1]"), MshEnt("D2L").UNV))
%%
%[text] ## DoF
%check
% dof = [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
%         NdDoF("D2R1", MshEnt("D2T").msh, 1, [1/3, 2/3], 0)];
% dof.getDomn

% node = Node([0, 0; 1, 0; 1, 1; 0, 1]', [5, 6, 7, 8]);
% elem = Elem([2, 4, 1; 4, 2, 3]', [-5, 3, 1; 5, 4, 2]', [5, 7]);
% edge = Edge([1, 2; 3, 4; 4, 1; 2, 3; 4, 2]', [1, 0; 2, 0; 1, 0; 2, 0; -1, 2]', [3, 1, 4, 2, 0]);
% msh = Msh("D2T", node, elem, edge);
% dof = [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
%         NdDoF("D2", msh, 1, [0, 0.5, 1], [0; 0])];
% dof.getMsh
%%
dof = [NdDoF("D2", MshEnt("D2T").msh, 0, [], [0; 0]), ...
        NdDoF("D2", MshEnt("D2T").msh, 1, [1/3, 2/3], [0; 0]), ...
        NdDoF("D2", MshEnt("D2T").msh, 2, [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2; 1/3, 1/3]', [0; 0])];

dof.cumDoF
dof.cumDoF(1)
dof.cumDoF(2)
dof.cumDoF(3)

sub2ind(dof, 1)
sub2ind(dof, 2)

sub2ind(dof, 1, 2, 1)
sub2ind(dof, 2, 2, 1)
sub2ind(dof, 2, 1, 2)

[iDoF, iEnt, iSamp] = ind2sub(dof, 2)
[iDoF, iEnt, iSamp] = ind2sub(dof, 5)
[iDoF, iEnt, iSamp] = ind2sub(dof, 7)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
