%[text] # MshEnt - symbolic mesh entity
%[text] ## D2, D2R, D2R1
disp(MshEnt("D2").var);
disp(MshEnt("D2R").var);
disp(MshEnt("D2R1").var);
%%
%[text] ## D2T
disp(MshEnt("D2T").var);
disp(MshEnt("D2T").parm);

disp(MshEnt("D2T").msh);
disp(MshEnt("D2T").node.coord);
disp(MshEnt("D2T").elem.node);
disp(MshEnt("D2T").edge.node);
%[text] ### UNV
disp(MshEnt("D2T").UNV(2).domn);
disp(MshEnt("D2T").UNV(2).fun);
%[text] ### **UTV**
disp(MshEnt("D2T").UTV(2).domn);
disp(MshEnt("D2T").UTV(2).fun);
%[text] ### len
disp(MshEnt("D2T").len(2).domn);
disp(MshEnt("D2T").len(2).fun);
%%
%[text] ## D2TR
disp(MshEnt("D2TR").var);
disp(MshEnt("D2TR").parm);

disp(MshEnt("D2T").msh);
disp(MshEnt("D2T").node.coord);
disp(MshEnt("D2T").elem.node);
disp(MshEnt("D2T").edge.node);
%%
%[text] ## D2L
disp(MshEnt("D2L").var);
disp(MshEnt("D2L").parm);

disp(MshEnt("D2L").node.coord);
disp(MshEnt("D2L").edge.node);
%[text] ### UNV
disp(MshEnt("D2L").UNV.domn);
disp(MshEnt("D2L").UNV.fun);
%[text] ### UTV
disp(MshEnt("D2L").UTV.domn);
disp(MshEnt("D2L").UTV.fun);
%[text] ### len
disp(MshEnt("D2L").len.domn);
disp(MshEnt("D2L").len.fun);
%%
%[text] ## D2LR
disp(MshEnt("D2LR").var);
disp(MshEnt("D2LR").parm);

disp(MshEnt("D2LR").node.coord);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
