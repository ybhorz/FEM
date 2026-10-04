%[text] # Tfm - Affine Transform
%%
%[text] ## D2T $\\leftrightarrow$ D2TR
%[text] variable
disp(Tfm("D2T").orgVar);
disp(Tfm("D2T").refVar);
%[text] transform
disp(Tfm("D2T").toRef);
disp(Tfm("D2T").toOrg);
%[text] Jacobian determinant
disp(Tfm("D2T").JDet.fun);
%%
%[text] ### Function: D2T $\\rightarrow$ D2TR
tfmFcn = Tfm("D2T").refTfm;
disp(tfmFcn.fun);

disp(tfmFcn.eval("[x1;y1]"));
disp(tfmFcn.eval("[x2;y2]"));
disp(tfmFcn.eval("[x3;y3]"));
%%
%[text] ### Fucntion: D2TR $\\rightarrow$ D2T
tfmFcn = Tfm("D2T").orgTfm;
disp(tfmFcn.fun);

disp(tfmFcn.eval("[0;0]"));
disp(tfmFcn.eval("[1;0]"));
disp(tfmFcn.eval("[0;1]"));
%%
%[text] ## D2LR $\\rightarrow$ D2L
%[text] variable
disp(Tfm("D2L").orgVar);
disp(Tfm("D2L").refVar);
%[text] transform
disp(Tfm("D2L").toOrg);
%[text] Jacobian Norm
disp(Tfm("D2L").JNorm.fun);
%%
%[text] ### Function: D2LR $\\rightarrow$ D2L
tfmFcn = Tfm("D2L").orgTfm;
disp(tfmFcn.fun);

disp(tfmFcn.eval("0"));
disp(tfmFcn.eval("1/2"));
disp(tfmFcn.eval("1"));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
