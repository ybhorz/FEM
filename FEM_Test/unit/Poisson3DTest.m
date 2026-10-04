classdef Poisson3DTest < matlab.unittest.TestCase
    % Poisson3DTest: convergence of 3D Poisson solvers (CG, CR, DG) and smoke test of 3D examples.
    % Exact solution u = sin(pi*x)*sin(pi*y)*sin(pi*z) + x*y*z on unit cube, mesh mshD3TS with nSub subdivisions.
    % Dirichlet boundary: faces x = 0, y = 0, z = 0; Neumann boundary: faces x = 1, y = 1, z = 1.

    % Scheme | Element | nSub  | Expected rate (L2, H1) | Lower bound of rate (L2, H1)
    % -------|---------|-------|------------------------|-----------------------------
    % CG     | P1      | 8, 16 | 2, 1                   | 1.7, 0.8
    % CG     | P2      | 4, 8  | 3, 2                   | 2.7, 1.8
    % CR     | P1      | 8, 16 | 2, 1                   | 1.7, 0.8
    % DG     | P1      | 8, 16 | 2, 1                   | 1.7, 0.8
    % MFE    | RT0-P0  | 4, 8  | 1, 1, 1 (|u-uh|_L2, |u-uh|_div, |p-ph|_L2) | 0.8, 0.8, 0.8
    % MFE    | BDM1-P0 | 4, 8  | 2, 1, 1                                    | 1.7, 0.8, 0.8
    % HDG    | P1-P1-P1| 4, 8  | 2, 2, 1.5 (|u-uh|_L2, |p-ph|_L2, |p-lh|_L2) | 1.7, 1.7, 1.3
    % |p-lh|_L2 is the unscaled L2 norm on all faces (total area ~ 1/h), hence order 2 - 1/2.
    % Measured on nSub = 4, 8 (R2025a): CG P1 1.70, 0.87; CG P2 2.97, 1.87; CR 1.91, 0.98; DG P1 1.74, 0.92.

    methods (Test, TestTags = {'Slow'})
        function CG1(tc)
            verifyRate(tc, "CG", 1, [8, 16], [1.7, 0.8]);
        end
        function CG2(tc)
            verifyRate(tc, "CG", 2, [4, 8], [2.7, 1.8]);
        end
        function CR(tc)
            verifyRate(tc, "CR", 1, [8, 16], [1.7, 0.8]);
        end
        function DG1(tc)
            verifyRate(tc, "DG", 1, [8, 16], [1.7, 0.8]);
        end
        function MFE(tc)
            elems = ["RT0", "BDM1"];
            minRates = [0.8, 0.8, 0.8; 1.7, 0.8, 0.8];
            nSubs = [4, 8];
            for iElem = 1:2
                err = zeros(2, 3);
                for iSub = 1:2
                    tic;
                    err(iSub, :) = solveMFE(elems(iElem), nSubs(iSub));
                    fprintf("MFE %s-P0, nSub = %d: |u-uh|_L2 = %e, |u-uh|_div = %e, |p-ph|_L2 = %e (%.1f s)\n", ...
                        elems(iElem), nSubs(iSub), err(iSub, :), toc);
                end
                rate = log2(err(1, :) ./ err(2, :));
                fprintf("MFE %s-P0: rate |u-uh|_L2 = %.2f, |u-uh|_div = %.2f, |p-ph|_L2 = %.2f\n", elems(iElem), rate);
                tc.verifyGreaterThan(rate, minRates(iElem, :), elems(iElem));
                if iElem == 1
                    % Example (RT0, nSub = 4) prints the same errors.
                    exDir = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "FEM_Example");
                    tc.verifyEqual(runExample(exDir, "Poisson_MFE_3D"), err(1, :), "RelTol", 1e-5);
                end
            end
        end
        function HDG(tc)
            nSubs = [4, 8];
            err = zeros(2, 3);
            for iSub = 1:2
                tic;
                err(iSub, :) = solveHDG(nSubs(iSub));
                fprintf("HDG, nSub = %d: |u-uh|_L2 = %e, |p-ph|_L2 = %e, |p-lh|_L2 = %e (%.1f s)\n", nSubs(iSub), err(iSub, :), toc);
            end
            rate = log2(err(1, :) ./ err(2, :));
            fprintf("HDG: rate |u-uh|_L2 = %.2f, |p-ph|_L2 = %.2f, |p-lh|_L2 = %.2f\n", rate);
            tc.verifyGreaterThan(rate, [1.7, 1.7, 1.3]);
            % Example (nSub = 4) prints the same errors.
            exDir = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "FEM_Example");
            tc.verifyEqual(runExample(exDir, "Poisson_HDG_3D"), err(1, :), "RelTol", 1e-5);
        end
        function example(tc)
            % 3D examples (nSub = 4) run and print the same errors as `solvePoisson` with nSub = 4.
            exDir = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "FEM_Example");
            exName = ["Poisson_CG_3D", "Poisson_CR_3D", "Poisson_DG_3D"];
            scheme = ["CG", "CR", "DG"];
            for iEx = 1:length(exName)
                err = runExample(exDir, exName(iEx));
                fprintf("%s: %s\n", exName(iEx), sprintf("%e ", err));
                refErr = zeros(1, 2);
                [refErr(1), refErr(2)] = solvePoisson(scheme(iEx), 1, 4);
                tc.verifyEqual(err, refErr, "RelTol", 1e-5, exName(iEx));
            end
        end
    end
end
%% Local functions.
function verifyRate(tc, scheme, ord, nSubs, minRate)
    % verifyRate: verify convergence rates of L2 and H1 errors between two meshes.
    err = zeros(length(nSubs), 2);
    for iSub = 1:length(nSubs)
        tic;
        [err(iSub, 1), err(iSub, 2)] = solvePoisson(scheme, ord, nSubs(iSub));
        fprintf("%s P%d, nSub = %d: |u-uh|_L2 = %e, |u-uh|_H1 = %e (%.1f s)\n", scheme, ord, nSubs(iSub), err(iSub, 1), err(iSub, 2), toc);
    end
    rate = log2(err(1, :) ./ err(2, :)) / log2(nSubs(2) / nSubs(1));
    fprintf("%s P%d: rate L2 = %.2f, rate H1 = %.2f\n", scheme, ord, rate(1), rate(2));
    tc.verifyGreaterThan(rate(1), minRate(1), sprintf("%s P%d: L2 rate.", scheme, ord));
    tc.verifyGreaterThan(rate(2), minRate(2), sprintf("%s P%d: H1 rate.", scheme, ord));
end
function [errL2, errH1] = solvePoisson(scheme, ord, nSub)
    % solvePoisson: solve 3D Poisson equation and return error norms (same schemes as 3D examples).
    u = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z) + x*y*z");
    d0 = [0; 0; 0]; grad = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
    f = - dif(u, [2; 0; 0]) - dif(u, [0; 2; 0]) - dif(u, [0; 0; 2]);
    UNV = MshEnt("D3F").UNV;
    g_N = dot(UNV, dif(u, grad));
    msh = mshD3TS([0, 1, 0, 1, 0, 1], nSub);
    DirBd = [1, 3, 5];
    NeuBd = [2, 4, 6];
    switch scheme
        case "CG"
            switch ord
                case 1
                    fE = stdFE("P1");
                    bC = BC(u, "node", msh.bdEnt(0, DirBd));
                case 2
                    fE = stdFE("P2");
                    bC = BC(u, "node", msh.bdEnt(0, DirBd), "edge", msh.bdEnt(1, DirBd));
            end
        case "CR"
            fE = stdFE("CR");
            bC = BC(u, "face", msh.bdEnt(2, DirBd));
        case "DG"
            fE = stdFE("DG1");
            bC = BC(u, "node", msh.bdEnt(0, DirBd));
    end
    Uh = FES(msh, fE, bC);
    Auv = DLF(msh, 3, Fcn.cst(1), grad, grad, "GInt", GInt("D3T", (ord - 1) * 2));
    H1Norm = Norm(msh, 3, grad, "GInt", GInt("D3T", ord * 2));
    if isequal(scheme, "DG")
        IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, [DirBd, NeuBd]));
        hF = MshEnt("D3F").area .^ (1/2);
        gamma = 20;
        Auv = [Auv, ...
            DLF.interface(msh, 2, - UNV, grad, d0, "EntIdx", IntFace, "trlOpr", "aver", "tstOpr", "jump", "GInt", GInt("D3F", 2 * ord - 1)), ...
            DLF.interface(msh, 2, - UNV, d0, grad, "EntIdx", IntFace, "trlOpr", "jump", "tstOpr", "aver", "GInt", GInt("D3F", 2 * ord - 1)), ...
            DLF.interface(msh, 2, hF \ gamma, d0, d0, "EntIdx", IntFace, "trlOpr", "jump", "tstOpr", "jump", "GInt", GInt("D3F", 2 * ord))];
        H1Norm = [H1Norm, Norm(msh, 2, d0, "EntIdx", IntFace, "coef", hF \ 1, "fcnOpr", "jump", "GInt", GInt("D3F", 2 * ord))];
    end
    Fv = [SLF(msh, 3, f, d0, "GInt", GInt("D3T", 1 + ord)), ...
        SLF(msh, 2, g_N, d0, "EntIdx", msh.bdEnt(2, NeuBd), "GInt", GInt("D3F", 1 + ord))];
    [Stiff, Load] = assemble(msh, Uh, Uh, Auv, Fv);
    uh = FEF(Uh, Stiff \ Load);
    errL2 = eNorm(msh, u, uh, Norm(msh, 3, d0, "GInt", GInt("D3T", ord * 2)));
    errH1 = eNorm(msh, u, uh, H1Norm);
end
function err = solveMFE(elem, nSub)
    % solveMFE: solve 3D Poisson equation in mixed form by RT0-P0 or BDM1-P0 (same scheme as Poisson_MFE_3D) and return
    % errors [|u-uh|_L2, |u-uh|_div, |p-ph|_L2].
    p = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z) + x*y*z");
    d0_p = [0; 0; 0];
    u = dif(p, cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]));
    d0_u = zeros(3); div_u = eye(3);
    f = -sum(dif(u, div_u));
    UNV = MshEnt("D3F").UNV;
    msh = mshD3TS([0, 1, 0, 1, 0, 1], nSub);
    DirFace = msh.bdEnt(2, [1, 3, 5]);
    NeuFace = msh.bdEnt(2, [2, 4, 6]);
    Uh = FES(msh, stdFE(elem), BC(dot(u, UNV), "face", NeuFace));
    Ph = FES(msh, stdFE("P0"));
    Auv = [DLF(msh, 3, Fcn.cst(1), d0_u, d0_u, "iTrl", 1, "iTst", 1, "GInt", GInt("D3T", 2)), ...
        DLF(msh, 3, Fcn.cst(1), d0_p, div_u, "iTrl", 2, "iTst", 1, "GInt", GInt("D3T", 0)), ...
        DLF(msh, 3, Fcn.cst(-1), div_u, d0_p, "iTrl", 1, "iTst", 2, "GInt", GInt("D3T", 0))];
    Fv = [SLF(msh, 2, p .* UNV, d0_u, "iTst", 1, "EntIdx", DirFace, "GInt", GInt("D3F", 2)), ...
        SLF(msh, 3, f, d0_p, "iTst", 2, "GInt", GInt("D3T", 1))];
    [Stiff, Load] = assemble(msh, [Uh, Ph], [Uh, Ph], Auv, Fv);
    [uh, ph] = FEF.multi([Uh, Ph], Stiff \ Load);
    divForm = @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow;
    err = [eNorm(msh, u, uh, Norm(msh, 3, d0_u, "GInt", GInt("D3T", 4))), ...
        eNorm(msh, u, uh, Norm(msh, 3, div_u, "form", divForm, "GInt", GInt("D3T", 2))), ...
        eNorm(msh, p, ph, Norm(msh, 3, d0_p, "GInt", GInt("D3T", 2)))];
end
function err = solveHDG(nSub)
    % solveHDG: solve 3D Poisson equation in mixed form by HDG with discontinuous P1 for u and p and P1 trace
    % (same scheme as Poisson_HDG_3D) and return errors [|u-uh|_L2, |p-ph|_L2, |p-lh|_L2].
    % Remark: the 2D example Posisson_HDG uses the flux u.n + tau (p - lambda) of the positive side for both elements
    % of a face, i.e. stabilization +tau on one side and -tau on the other; on 3D meshes this gave first-order
    % convergence with large errors, so the standard flux with tau > 0 on every element is used here.
    p = Fcn("D3", "sin(pi*x)*sin(pi*y)*sin(pi*z) + x*y*z");
    d0_p = [0; 0; 0];
    grad_p = cat(3, [1; 0; 0], [0; 1; 0], [0; 0; 1]);
    u = dif(p, grad_p);
    d0_u = zeros(3); div_u = eye(3);
    f = -sum(dif(u, div_u));
    UNV = MshEnt("D3F").UNV;
    msh = mshD3TS([0, 1, 0, 1, 0, 1], nSub);
    DirFace = msh.bdEnt(2, [1, 3, 5]);
    NeuFace = msh.bdEnt(2, [2, 4, 6]);
    D3TElem = MshEnt("D3T").msh;
    Ph = FES(msh, stdFE("DG1"));
    Uh = FES(msh, FE("D3T", FE.repFS("[1,x,y,z]", [3, 1]), ...
        [NdDoF("D3", D3TElem, 0, [], [0, nan, nan; 0, nan, nan; 0, nan, nan], "share", false), ...
        NdDoF("D3", D3TElem, 0, [], [nan, 0, nan; nan, 0, nan; nan, 0, nan], "share", false), ...
        NdDoF("D3", D3TElem, 0, [], [nan, nan, 0; nan, nan, 0; nan, nan, 0], "share", false)], "map", "affine"));
    Lh = FES(msh, stdFE("TP1"), BC(p, "face", DirFace));
    d0_l = [0; 0];
    tau = 1;
    trls = [Uh, Ph, Lh];
    IntFace = setdiff(1:msh.nFace, msh.bdEnt(2, 1:6));
    % Numerical flux on each element: u|_K . n_K - tau (p|_K - lambda); sign of iTrl / iTst selects the side of face.
    Auv = [DLF(msh, 3, Fcn.cst(1), d0_u, d0_u, "iTrl", 1, "iTst", 1, "GInt", GInt("D3T", 2)), ...
        DLF(msh, 3, Fcn.cst(1), d0_p, div_u, "iTrl", 2, "iTst", 1, "GInt", GInt("D3T", 1)), ...
        DLF.interface(msh, 2, -UNV, d0_l, d0_u, "iTrl", 3, "iTst", 1, "tstOpr", "jump", "GInt", GInt("D3F", 2)), ...
        DLF(msh, 3, Fcn.cst(1), d0_u, grad_p, "iTrl", 1, "iTst", 2, "GInt", GInt("D3T", 1)), ...
        DLF(msh, 2, -UNV, d0_u, d0_p, "iTrl", 1, "iTst", 2, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, UNV, d0_u, d0_p, "iTrl", -1, "iTst", -2, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(tau), d0_p, d0_p, "iTrl", 2, "iTst", 2, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(tau), d0_p, d0_p, "iTrl", -2, "iTst", -2, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(-tau), d0_l, d0_p, "iTrl", 3, "iTst", 2, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(-tau), d0_l, d0_p, "iTrl", 3, "iTst", -2, "GInt", GInt("D3F", 2)), ...
        DLF.interface(msh, 2, UNV, d0_u, d0_l, "iTrl", 1, "iTst", 3, "trlOpr", "jump", "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(-tau), d0_p, d0_l, "iTrl", 2, "iTst", 3, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(-tau), d0_p, d0_l, "iTrl", -2, "iTst", 3, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(tau), d0_l, d0_l, "iTrl", 3, "iTst", 3, "EntIdx", NeuFace, "GInt", GInt("D3F", 2)), ...
        DLF(msh, 2, Fcn.cst(2 * tau), d0_l, d0_l, "iTrl", 3, "iTst", 3, "EntIdx", IntFace, "GInt", GInt("D3F", 2))];
    Fv = [SLF(msh, 3, f, d0_p, "iTst", 2, "GInt", GInt("D3T", 2)), ...
        SLF(msh, 2, dot(u, UNV), d0_l, "iTst", 3, "EntIdx", NeuFace, "GInt", GInt("D3F", 2))];
    [Stiff, Load] = assemble(msh, trls, trls, Auv, Fv);
    [uh, ph, lh] = FEF.multi(trls, condSolve(Stiff, Load, trls, [1, 2]));
    err = [eNorm(msh, u, uh, Norm(msh, 3, d0_u, "GInt", GInt("D3T", 4))), ...
        eNorm(msh, p, ph, Norm(msh, 3, d0_p, "GInt", GInt("D3T", 4))), ...
        eNorm(msh, p, lh, Norm(msh, 2, d0_l, "GInt", GInt("D3F", 4)))];
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
