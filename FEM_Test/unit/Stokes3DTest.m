classdef Stokes3DTest < matlab.unittest.TestCase
    % Stokes3DTest: convergence of 3D Stokes solver (SDG on Alfeld split mesh) and smoke test of the 3D example.
    % Exact solution u = curl(phi, phi, phi), phi = (sin(pi*x)*sin(pi*y)*sin(pi*z))^2, p = 10*(x-1/2)*(y-1/2)*(z-1/2)
    % on unit cube, u = 0 on boundary, nu = 1, mesh mshSplit(mshD3TS) with nSub subdivisions.

    % Scheme | Element | nSub  | Expected rate (|s-sh|_L2, |u-uh|_L2, |u-uh|_H1, |p-ph|_L2) | Lower bound of rate
    % -------|---------|-------|----------------------------------------------------------|--------------------
    % SDG    | SDG_1   | 4, 8  | 2, 2, 1, 2                                               | 1.7, 1.7, 0.8, 1.5

    methods (Test, TestTags = {'Slow'})
        function SDG(tc)
            nSubs = [4, 8];
            err = zeros(2, 4);
            for iSub = 1:2
                tic;
                err(iSub, :) = solveStokes(nSubs(iSub));
                fprintf("Stokes SDG_1, nSub = %d: |s-sh|_L2 = %e, |u-uh|_L2 = %e, |u-uh|_H1 = %e, |p-ph|_L2 = %e (%.1f s)\n", ...
                    nSubs(iSub), err(iSub, :), toc);
            end
            rate = log2(err(1, :) ./ err(2, :));
            fprintf("Stokes SDG_1: rate |s-sh|_L2 = %.2f, |u-uh|_L2 = %.2f, |u-uh|_H1 = %.2f, |p-ph|_L2 = %.2f\n", rate);
            tc.verifyGreaterThan(rate, [1.7, 1.7, 0.8, 1.5]);
            % Example (nSub = 4) prints the same errors.
            exDir = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "FEM_Example");
            tc.verifyEqual(runExample(exDir, "Stokes_SDG_3D"), err(1, :), "RelTol", 1e-5);
        end
    end
end
%% Local functions.
function err = solveStokes(nSub)
    % solveStokes: solve 3D Stokes equation by SDG_1 (same scheme as Stokes_SDG_3D) and return errors
    % [|s-sh|_L2, |u-uh|_L2, |u-uh|_H1, |p-ph|_L2] (broken H1 norm, see the example).
    d0_p = [0; 0; 0]; grad_p = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
    phi = Fcn("D3", "sin(pi*x)^2*sin(pi*y)^2*sin(pi*z)^2").dif(grad_p).fun;
    u = Fcn("D3", [phi(2) - phi(3); phi(3) - phi(1); phi(1) - phi(2)]);
    p = Fcn("D3", "10*(x-1/2)*(y-1/2)*(z-1/2)");
    d0_u = zeros(3); div_u = eye(3);
    grad_u = cat(3, [1, 1, 1; 0, 0, 0; 0, 0, 0], [0, 0, 0; 1, 1, 1; 0, 0, 0], [0, 0, 0; 0, 0, 0; 1, 1, 1]);
    s = dif(u, grad_u);
    d0_s = zeros(3, 9); div_s = repelem(eye(3), 1, 3);
    f = -sum(dif(s, div_s), 2) + dif(p, grad_p);
    UNV = MshEnt("D3F").UNV;
    hF = MshEnt("D3F").area .^ (1/2);
    msh = mshSplit(mshD3TS([0, 1, 0, 1, 0, 1], nSub));
    PrFace = find(ismember(msh.face.type, 0:6));
    DlFace = find(ismember(msh.face.type, 1i));
    DlEdge = find(ismember(msh.edge.type, 1i));
    Sh = FES(msh, stdFE("StSDG1M"));
    Uh = FES(msh, stdFE("StSDG1V"), BC(Fcn("D3", "0"), "node", nan, "face", msh.bdEnt(2, 1:6)));
    Ph = FES(msh, stdFE("StSDG1S"), BC(p, "node", nan, "edge", DlEdge(1)));
    trls = [Sh, Uh, Ph];
    Auv = [DLF(msh, 3, Fcn.cst(1), d0_s, d0_s, "iTrl", 1, "iTst", 1, "GInt", GInt("D3T", 2)), ...
        DLF(msh, 3, Fcn.cst(1), d0_u, div_s, "iTrl", 2, "iTst", 1, "form", @(coef, trl, tst) coef .* dot(trl, sum(tst, 2)), "GInt", GInt("D3T", 2)), ...
        DLF.interface(msh, 2, -UNV, d0_u, d0_s, "EntIdx", PrFace, "iTrl", 2, "iTst", 1, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl, tst * coef), "GInt", GInt("D3F", 2)), ...
        DLF(msh, 3, Fcn.cst(1), d0_s, grad_u, "iTrl", 1, "iTst", 2, "GInt", GInt("D3T", 2)), ...
        DLF.interface(msh, 2, -UNV, d0_s, d0_u, "EntIdx", DlFace, "iTrl", 1, "iTst", 2, "tstOpr", "jump", ...
        "form", @(coef, trl, tst) dot(trl * coef, tst), "GInt", GInt("D3F", 2)), ...
        DLF(msh, 3, Fcn.cst(-1), d0_p, div_u, "iTrl", 3, "iTst", 2, "GInt", GInt("D3T", 2)), ...
        DLF.interface(msh, 2, UNV, d0_p, d0_u, "EntIdx", DlFace, "iTrl", 3, "iTst", 2, "tstOpr", "jump", "GInt", GInt("D3F", 2)), ...
        DLF(msh, 3, Fcn.cst(-1), d0_u, grad_p, "iTrl", 2, "iTst", 3, "GInt", GInt("D3T", 2)), ...
        DLF.interface(msh, 2, UNV, d0_u, d0_p, "EntIdx", PrFace, "iTrl", 2, "iTst", 3, "tstOpr", "jump", "GInt", GInt("D3F", 2))];
    Fv = SLF(msh, 3, f, d0_u, "iTst", 2, "GInt", GInt("D3T", 3));
    [Stiff, Load] = assemble(msh, trls, trls, Auv, Fv);
    [sh, uh, ph] = FEF.multi(trls, condSolve(Stiff, Load, trls, 1));
    err = [eNorm(msh, s, sh, Norm(msh, 3, d0_s, "GInt", GInt("D3T", 4))), ...
        eNorm(msh, u, uh, Norm(msh, 3, d0_u, "GInt", GInt("D3T", 4))), ...
        eNorm(msh, u, uh, [Norm(msh, 3, grad_u, "GInt", GInt("D3T", 2)), ...
        Norm(msh, 2, d0_u, "EntIdx", DlFace, "coef", hF \ 1, "form", @(coef, fcn, pow) coef .* sum(abs(fcn) .^ pow), ...
        "fcnOpr", "jump", "GInt", GInt("D3F", 2))]), ...
        eNorm(msh, p, ph, Norm(msh, 3, d0_p, "GInt", GInt("D3T", 4)))];
end
function err = runExample(exDir, exName)
    % runExample: run example script and extract printed error norms, e.g. "|u-uh|_L2: 1.234567e-03".
    oldDir = cd(exDir);
    cleanup = onCleanup(@() cd(oldDir));
    exOut = evalc("runScript(fullfile(exDir, exName + "".m""))");
    close all;
    tok = regexp(exOut, "\|[^|]+\|_\w+:\s*([-+.\deEinfINFnaN]+)", "tokens");
    err = cellfun(@(t) str2double(t{1}), tok);
end
function runScript(exFile)
    % runScript: run script in a separate workspace.
    run(exFile);
end
