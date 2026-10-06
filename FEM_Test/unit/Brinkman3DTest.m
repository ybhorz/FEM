classdef Brinkman3DTest < matlab.unittest.TestCase
    % Brinkman3DTest: convergence of hybridized 3D Brinkman SDG solver (Zhao, Chung and Lam, CMAME 2020) on Alfeld split
    % meshes, and smoke test of the 3D example.
    % Exact solution u = curl(phi, phi, phi), phi = (sin(pi*x)*sin(pi*y)*sin(pi*z))^2, p = 10*(x-1/2)*(y-1/2)*(z-1/2)
    % on unit cube, u = 0 on boundary, alpha = 10, s = grad(u), mesh mshSplit(mshD3TS) with nSub subdivisions.
    % Pressure error is measured with mean removed (pressure is fixed at one DoF).

    % nu   | nSub | Expected rate (|s-sh|_L2, |u-uh|_L2, |u-uh|_H1, |p-ph|_L2) | Lower bound of rate
    % -----|------|-------------------------------------------------------------|--------------------
    % 1    | 4, 8 | 2, 2, 1, 2                                                  | 1.7, 1.7, 0.8, 1.7
    % 1e-4 | 4, 8 | -, 2, -, 2 (robust in the Darcy limit for u and p)          | -, 1.5, -, 1.5

    methods (Test, TestTags = {'Slow'})
        function SDG(tc)
            nSubs = [4, 8];
            vLst = [1, 1e-4];
            minRates = {[1.7, 1.7, 0.8, 1.7], [nan, 1.5, nan, 1.5]};
            for iV = 1:2
                err = zeros(2, 4);
                for iSub = 1:2
                    tic;
                    err(iSub, :) = solveBrinkman(nSubs(iSub), vLst(iV));
                    fprintf("Brinkman SDG_1, nu = %g, nSub = %d: |s-sh|_L2 = %e, |u-uh|_L2 = %e, |u-uh|_H1 = %e, |p-ph|_L2 = %e (%.1f s)\n", ...
                        vLst(iV), nSubs(iSub), err(iSub, :), toc);
                end
                rate = log2(err(1, :) ./ err(2, :));
                fprintf("Brinkman SDG_1, nu = %g: rate |s-sh|_L2 = %.2f, |u-uh|_L2 = %.2f, |u-uh|_H1 = %.2f, |p-ph|_L2 = %.2f\n", vLst(iV), rate);
                isChk = ~isnan(minRates{iV});
                tc.verifyGreaterThan(rate(isChk), minRates{iV}(isChk), sprintf("nu = %g", vLst(iV)));
                if iV == 1
                    % Example (nu = 1, nSub = 4) prints the same errors.
                    exDir = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "FEM_Example");
                    tc.verifyEqual(runExample(exDir, "Brinkman_SDG_3D"), err(1, :), "RelTol", 1e-5);
                end
            end
        end
    end
end
%% Local functions.
function err = solveBrinkman(nSub, v)
    % solveBrinkman: solve 3D Brinkman equation by hybridized SDG_1 (same scheme as Brinkman_SDG_3D) and
    % return errors [|s-sh|_L2, |u-uh|_L2, |u-uh|_H1, |p-ph|_L2].
    sys = brinkmanSDG3D(nSub, v);
    msh = sys.msh;
    [sh, uh, ~, ~, ph] = FEF.multi(sys.trls, condSolve(sys.Stiff, sys.Load, sys.trls, {1, [2, 3]}));
    hF = MshEnt("D3F").area .^ (1/2);
    grad_u = cat(3, [1, 1, 1; 0, 0, 0; 0, 0, 0], [0, 0, 0; 1, 1, 1; 0, 0, 0], [0, 0, 0; 0, 0, 0; 1, 1, 1]);
    H1Norm = [Norm(msh, 3, grad_u, "GInt", GInt("D3T", 2)), ...
        Norm(msh, 2, zeros(3), "coef", hF \ 1, "form", @(coef, fcn, pow) coef .* sum(abs(fcn) .^ pow), "fcnOpr", "jump", "GInt", GInt("D3F", 2))];
    err = [eNorm(msh, sys.s, sh, Norm(msh, 3, zeros(3, 9), "GInt", GInt("D3T", 4))), ...
        eNorm(msh, sys.u, uh, Norm(msh, 3, zeros(3), "GInt", GInt("D3T", 4))), ...
        eNorm(msh, sys.u, uh, H1Norm), ...
        pErr(msh, sys.p, ph)];
end
function err = pErr(msh, p, ph)
    % pErr: L2 error of pressure with mean removed (pressure fixed at one DoF; |Omega| = 1).
    eMean = eNorm(msh, p, ph, Norm(msh, 3, [0; 0; 0], "pow", 1, "form", @(coef, fcn, pow) fcn, "GInt", GInt("D3T", 4)));
    err = sqrt(eNorm(msh, p, ph, Norm(msh, 3, [0; 0; 0], "GInt", GInt("D3T", 4))) ^ 2 - eMean ^ 2);
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
