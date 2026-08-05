%[text] # Msh - mesh data structure
%%
%[text] ## D2T mesh
%[text] ### Node
node = Node([0, 1; 1, 1; 0, 2; 1, 2]', [-1, -2, -4, -3]);
disp(node.dim);
disp(node.nNode);
%[text] ### Elem
elem = Elem([2, 3, 1; 3, 2, 4]', [5, -3, 1; -5, 4, -2]', [5, 7]);
disp(elem.nElem);
disp(elem.nNode);
disp(elem.nEdge);
%[text] ### Edge
edge = Edge([1, 2; 3, 4; 1, 3; 2, 4; 2, 3]', [1, 0; -2, 0; -1, 0; 2, 0; 1, -2]', [1, 3, 4, 2, 0]);
disp(edge.nEdge);
disp(edge.nNode);
disp(edge.nElem);
%[text] ### Msh
msh = Msh("D2T", node, elem, edge);
msh.figure;

disp(msh.dim);

disp(msh.nNode);
disp(msh.nElem);
disp(msh.nEdge);

disp(msh.nEnt(0));
disp(msh.nEnt(1));
disp(msh.nEnt(2));
%%
%[text] ## D1 mesh
msh1D = Msh("D1", Node([0, 1, 2]), Elem([0, 1; 1, 2]'));
disp(msh1D.dim);

disp(msh.nNode);
disp(msh.nElem);

disp(msh1D.nEnt(0));
disp(msh1D.nEnt(1));
%%
%[text] # mshD2TS - structured triangulation of a rectangle
msh = mshD2TS([0, 1, 3, 4], [2, 3]);
msh.figure;

disp(msh.dim);

disp(msh.nNode);
disp(msh.nElem);
disp(msh.nEdge);

disp(msh.node.coord);
disp(msh.node.type);

disp(msh.elem.node);
disp(msh.elem.edge);
disp(msh.elem.type);

disp(msh.edge.node);
disp(msh.edge.elem);
disp(msh.edge.type);
%%
%[text] # Msh.auto - automatically complete `mesh data`
msh0 = mshD2TS([0, 1, 3, 4], [2, 3]);
msh = Msh.auto("D2T", msh0.node, msh0.elem.node);
msh.figure;

disp(msh.dim);

disp(msh.nNode);
disp(msh.nElem);
disp(msh.nEdge);

disp(msh.node.coord);
disp(msh.node.type);

disp(msh.elem.node);
disp(msh.elem.edge);
disp(msh.elem.type);

disp(msh.edge.node);
disp(msh.edge.elem);
disp(msh.edge.type);
%%
%[text] # mshSplit - Alfeld refinement
msh0 = mshD2TS([0, 1, 3, 4], [2, 3]);
msh = mshSplit(msh0);
msh.figure;

disp(msh.dim);

disp(msh.nNode);
disp(msh.nElem);
disp(msh.nEdge);

disp(msh.node.coord);
disp(msh.node.type);

disp(msh.elem.node);
disp(msh.elem.edge);
disp(msh.elem.type);

disp(msh.edge.node);
disp(msh.edge.elem);
disp(msh.edge.type);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
