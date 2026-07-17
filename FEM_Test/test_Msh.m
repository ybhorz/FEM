%[text] ## Node
Node(rand([2, 8]))
Node(rand([2, 8]), zeros(1, 8))
%%
%[text] ## Elem
Elem(zeros([3, 8]))
Elem(zeros([3, 8]), zeros([3, 8]))
Elem(zeros([3, 8]), zeros([3, 8]), zeros([1, 8]))
%%
%[text] ## Edge
Edge(zeros([2, 8]))
Edge(zeros([2, 8]), zeros([2, 8]))
Edge(zeros([2, 8]), zeros([2, 8]), zeros([1, 8]))
%%
%[text] ## Msh
Msh("D2T", Node(zeros(2, 8)), Elem(zeros(3, 8)), Edge(zeros(2, 8)))
Msh("D2T", Node(zeros(2, 8)), Elem(zeros(3, 8), zeros(3, 8)), Edge(zeros(2, 8), zeros(2, 8)))

Msh("D1", Node(zeros(1, 3)), Elem(zeros(2, 2)))
%%
%[text] msh.auto
msh1 = mshD2TS([0, 1, 3, 4], [2, 3]);
msh2 = Msh.auto("D2T", msh1.node, msh1.elem.node);

msh2.node.coord
msh2.node.type

msh2.elem.node
msh2.elem.edge
msh2.elem.type

msh2.edge.node
msh2.edge.elem
msh2.edge.type

msh2.figure
%%
% PDE toolbox.
NdCrd = p;
NdType = zeros(1, size(p, 2));
NdType(unique(e(1:2, e(5, :) == 3))) = 1;
NdType(unique(e(1:2, e(5, :) == 2))) = 2;
NdType(unique(e(1:2, e(5, :) == 1))) = 3;
NdType(unique(e(1:2, e(5, :) == 4))) = 4;
NdType(4) = -1; NdType(3) = -2; NdType(2) = -3; NdType(1) = -4;
ElNode = t(1:3,:);
msh = Msh.auto("D2T", Node(NdCrd, NdType), ElNode);
msh.figure;
%%
%[text] ## mshD2TS
msh = mshD2TS([0, 1, 3, 4], [2, 3]);

msh.node.coord
msh.node.type

msh.elem.node
msh.elem.edge
msh.elem.type

msh.edge.node
msh.edge.elem
msh.edge.type

msh.figure
%%
%[text] ## mshSplit
msh1 = mshD2TS([0, 1, 3, 4], [2, 3]);
msh2 = mshSplit(msh1);
msh1.figure;
msh2.figure
%%
%[text] ## MshEnt
MshEnt
%%
MshEnt("D2")
MshEnt("D2").var
%%
MshEnt("D2R")
MshEnt("D2R").var
%%
MshEnt("D2R1")
MshEnt("D2R1").var
%%
MshEnt("D2T")
MshEnt("D2T").msh
MshEnt("D2T").node.coord
MshEnt("D2T").elem.node
MshEnt("D2T").edge.node
MshEnt("D2T").var
MshEnt("D2T").parm
%%
MshEnt("D2TR")
MshEnt("D2TR").msh
MshEnt("D2TR").node.coord
MshEnt("D2TR").elem.node
MshEnt("D2TR").edge.node
MshEnt("D2TR").var
MshEnt("D2TR").parm
%%
MshEnt("D2L")
MshEnt("D2L").node.coord
MshEnt("D2L").edge.node
MshEnt("D2L").var
MshEnt("D2L").parm
%%
MshEnt("D2LR")
MshEnt("D2LR").msh
MshEnt("D2LR").node.coord
MshEnt("D2LR").elem.node
MshEnt("D2LR").var
MshEnt("D2LR").parm
%%
%[text] UNV
MshEnt("D2L").UNV.domn
MshEnt("D2L").UNV.fun

MshEnt("D2T").UNV.domn
MshEnt("D2T").UNV.fun
%%
%[text] UTV
MshEnt("D2L").UTV.domn
MshEnt("D2L").UTV.fun

MshEnt("D2T").UTV.domn
MshEnt("D2T").UTV.fun
%%
%[text] len
MshEnt("D2L").len.domn
MshEnt("D2L").len.fun
MshEnt("D2T").len.domn
MshEnt("D2T").len.fun

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
