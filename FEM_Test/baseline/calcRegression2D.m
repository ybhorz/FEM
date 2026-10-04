function data = calcRegression2D(options)
    % calcRegression2D: compute 2D regression data (shared by makeBaseline2D and Regression2DTest).

    % Field   | Content                                                   | Cost
    % ------- | --------------------------------------------------------- | ----
    % small   | FE spaces, matrices, projections and errors on small mesh | Fast
    % example | error norms printed by the 9 examples in FEM_Example      | Slow

    arguments (Input)
        options.small logical = true; % Whether compute data on small mesh.
        options.example logical = false; % Whether run examples.
        options.exName (1, :) string = ["Poisson_CG", "Poisson_CR", "Poisson_DG", "Posisson_MFE", "Posisson_HDG", ...
            "Nonlinear_Poisson", "Poisson_SDG", "Stokes_SDG", "Brinkman_SDG"]; % Examples to run.
    end
    arguments (Output)
        data struct;
    end
    data = struct();
    if options.small
        data.small = calcSmall();
    end
    if options.example
        data.example = calcExample(options.exName);
    end
end
%% Local functions.
function small = calcSmall()
    % calcSmall: regression data on small meshes.
    msh = mshD2TS([0, 1, 0, 1], 4);
    BdNode = find(msh.node.type ~= 0);
    BdEdge = find(msh.edge.type ~= 0);
    IntEdge = find(msh.edge.type == 0);
    D2TElem = MshEnt("D2T").msh;
    UNV = MshEnt("D2L").UNV;
    len = MshEnt("D2L").len;
    u = Fcn("D2", "sin(pi*x)*sin(pi*y) + x*y");
    w = Fcn("D2", "[sin(pi*x)*y; cos(pi*y)*x]");
    d0 = [0; 0]; grad = cat(3, [1; 0], [0; 1]);
    d0V = [0, 0; 0, 0]; divV = [1, 0; 0, 1]; gradV = cat(3, [1, 1; 0, 0], [0, 0; 1, 1]);
    divForm = @(coef, fcn, pow) abs(sum(coef .* fcn)) .^ pow;
    matForm = @(coef, fcn, pow) coef .* sum(abs(fcn) .^ pow);

    % P1: continuous, Dirichlet on boundary nodes.
    P1_FES = FES(msh, FE("D2T", "[1,x,y]", NdDoF("D2", D2TElem, 0, [], d0)), BC(u, "node", BdNode));
    small.P1 = calcFES(msh, P1_FES, u, d0, grad);
    [small.P1.Stiff, small.P1.Load] = assemble(msh, P1_FES, P1_FES, DLF(msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 0)), ...
        [SLF(msh, 2, u, d0, "GInt", GInt("D2T", 3)), SLF(msh, 1, u, d0, "EntIdx", BdEdge(1:2:end), "GInt", GInt("D2L", 3))]);

    % P2: continuous, Dirichlet on boundary nodes and edges.
    P2_FES = FES(msh, FE("D2T", "[1,x,y,x^2,y^2,x*y]", [NdDoF("D2", D2TElem, 0, [], d0), NdDoF("D2", D2TElem, 1, 0.5, d0)]), ...
        BC(u, "node", BdNode, "edge", BdEdge));
    small.P2 = calcFES(msh, P2_FES, u, d0, grad);
    [small.P2.Stiff, small.P2.Load] = assemble(msh, P2_FES, P2_FES, DLF(msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 2)), ...
        SLF(msh, 2, u, d0, "GInt", GInt("D2T", 4)));

    % CR: Dirichlet on boundary edges.
    CR_FES = FES(msh, FE("D2T", "[1,x,y]", NdDoF("D2", D2TElem, 1, 1/2, d0)), BC(u, "edge", BdEdge));
    small.CR = calcFES(msh, CR_FES, u, d0, grad);

    % DG1: SIPG.
    DG1_FES = FES(msh, FE("D2T", "[1,x,y]", NdDoF("D2", D2TElem, 0, [], d0, "share", false)), BC(u, "node", BdNode));
    small.DG1 = calcFES(msh, DG1_FES, u, d0, grad);
    Auv = [DLF(msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 0)), ...
        DLF.interface(msh, 1, -UNV, grad, d0, "EntIdx", IntEdge, "trlOpr", "aver", "tstOpr", "jump", "GInt", GInt("D2L", 1)), ...
        DLF.interface(msh, 1, -UNV, d0, grad, "EntIdx", IntEdge, "trlOpr", "jump", "tstOpr", "aver", "GInt", GInt("D2L", 1)), ...
        DLF.interface(msh, 1, len \ 10, d0, d0, "EntIdx", IntEdge, "trlOpr", "jump", "tstOpr", "jump", "GInt", GInt("D2L", 2))];
    [small.DG1.Stiff, small.DG1.Load] = assemble(msh, DG1_FES, DG1_FES, Auv, SLF(msh, 2, u, d0, "GInt", GInt("D2T", 3)));
    small.DG1.JumpErr = eNorm(msh, u, DG1_FES.proj(u), Norm(msh, 1, d0, "EntIdx", IntEdge, "coef", len \ 1, "fcnOpr", "jump", ...
        "GInt", GInt("D2L", 4)));

    % RT0 - P0: mixed finite element.
    RT0_FES = FES(msh, FE("D2T", "[1,0; 0,1; x,y].'", MoDoF("D2", D2TElem, 1, Fcn.cst(1), d0V, "coef", UNV, "orien", true, ...
        "GInt", GInt("D2L", 4))), BC(dot(w, UNV), "edge", BdEdge(1:2:end)));
    small.RT0 = calcFES(msh, RT0_FES, w, d0V, divV, divForm);
    P0_FES = FES(msh, FE("D2T", "1", MoDoF("D2", D2TElem, 2, Fcn.cst(1), d0, "GInt", GInt("D2T", 2))));
    small.P0 = calcFES(msh, P0_FES, u, d0, grad);
    trls = [RT0_FES, P0_FES];
    Auv = [DLF(msh, 2, Fcn.cst(1), d0V, d0V, "iTrl", 1, "iTst", 1, "GInt", GInt("D2T", 2)), ...
        DLF(msh, 2, Fcn.cst(1), d0, divV, "iTrl", 2, "iTst", 1, "GInt", GInt("D2T", 0)), ...
        DLF(msh, 2, Fcn.cst(-1), divV, d0, "iTrl", 1, "iTst", 2, "GInt", GInt("D2T", 0))];
    Fv = [SLF(msh, 1, u .* UNV, d0V, "iTst", 1, "EntIdx", BdEdge(2:2:end), "GInt", GInt("D2L", 3)), ...
        SLF(msh, 2, u, d0, "iTst", 2, "GInt", GInt("D2T", 2))];
    [small.MFE.Stiff, small.MFE.Load] = assemble(msh, trls, trls, Auv, Fv);

    % BDM1.
    BDM1_FES = FES(msh, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
        [MoDoF("D2", D2TElem, 1, Fcn.cst(1), d0V, "coef", UNV, "orien", true, "GInt", GInt("D2L", 4)), ...
        MoDoF("D2", D2TElem, 1, Fcn("D2R1", "s-1/2"), d0V, "coef", UNV, "orien", true, "GInt", GInt("D2L", 4))]));
    small.BDM1 = calcFES(msh, BDM1_FES, w, d0V, divV, divForm);

    % VP1, SP1 and TP1: HDG.
    VP1_FES = FES(msh, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
        [NdDoF("D2", D2TElem, 0, [], [0, nan; 0, nan], "share", false), NdDoF("D2", D2TElem, 0, [], [nan, 0; nan, 0], "share", false)]));
    small.VP1 = calcFES(msh, VP1_FES, w, d0V, gradV, matForm);
    TP1_FES = FES(msh, FE("D2LR", "[1,s]", NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0, "share", false)), ...
        BC(u, "edge", BdEdge(1:2:end)));
    small.TP1.nGlDoF = TP1_FES.nGlDoF;
    small.TP1.Lc2Gl = TP1_FES.Lc2Gl;
    small.TP1.BCIdx = TP1_FES.BC.DoFIdx;
    small.TP1.BCVal = TP1_FES.BC.DoFVal;
    small.TP1.Proj = projDoF(TP1_FES, u);
    small.TP1.Err = eNorm(msh, u, TP1_FES.proj(u), Norm(msh, 1, 0, "GInt", GInt("D2L", 4)));
    tau = 1;
    trls = [VP1_FES, DG1_FES, TP1_FES];
    iu = 1; ip = 2; il = 3;
    Auv = [DLF(msh, 2, Fcn.cst(1), d0V, d0V, "iTrl", iu, "iTst", iu, "GInt", GInt("D2T", 2)), ...
        DLF(msh, 2, Fcn.cst(1), d0, divV, "iTrl", ip, "iTst", iu, "GInt", GInt("D2T", 1)), ...
        DLF.interface(msh, 1, -UNV, 0, d0V, "iTrl", il, "iTst", iu, "tstOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF(msh, 2, Fcn.cst(1), d0V, grad, "iTrl", iu, "iTst", ip, "GInt", GInt("D2T", 1)), ...
        DLF.interface(msh, 1, -UNV, d0V, d0, "iTrl", iu, "iTst", ip, "tstOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF.interface(msh, 1, Fcn.cst(-tau), d0, d0, "iTrl", ip, "iTst", ip, "tstOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF.interface(msh, 1, Fcn.cst(tau), 0, d0, "iTrl", il, "iTst", ip, "tstOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF.interface(msh, 1, UNV, d0V, 0, "iTrl", iu, "iTst", il, "trlOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF.interface(msh, 1, Fcn.cst(tau), d0, 0, "iTrl", ip, "iTst", il, "trlOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF(msh, 1, -Fcn.cst(tau), 0, 0, "iTrl", il, "iTst", il, "EntIdx", BdEdge(2:2:end), "GInt", GInt("D2L", 2))];
    Fv = [SLF(msh, 2, u, d0, "iTst", ip, "GInt", GInt("D2T", 2)), ...
        SLF(msh, 1, u, 0, "iTst", il, "EntIdx", BdEdge(2:2:end), "GInt", GInt("D2L", 2))];
    [small.HDG.Stiff, small.HDG.Load] = assemble(msh, trls, trls, Auv, Fv);

    % Nonlinear: linearized forms with previous solution.
    uh0 = P1_FES.proj(u);
    Awuv = LDLF(msh, 2, Fcn.cst(-3), d0, d0, d0, "GInt", GInt("D2T", 3));
    Fwv = LSLF(msh, 2, Fcn.cst(3), d0, d0, "form", @(load, pre, tst) load .* (pre).^2 .* tst, "GInt", GInt("D2T", 3));
    [small.Nonlinear.Stiff, small.Nonlinear.Load] = assemble(msh, P1_FES, P1_FES, ...
        DLF(msh, 2, Fcn.cst(1), grad, grad, "GInt", GInt("D2T", 0)), SLF(msh, 2, u, d0, "GInt", GInt("D2T", 2)), ...
        "preSol", uh0, "Awuv", Awuv, "Fwv", Fwv);

    % SDG on Alfeld split mesh.
    mshS = mshSplit(mshD2TS([0, 1, 0, 1], 2));
    PrEdge = find(ismember(mshS.edge.type, [0, 1, 2, 3, 4]));
    DlEdge = find(ismember(mshS.edge.type, 1i));
    SDG_1S_FES = FES(mshS, FE("D2T", "[1,x,y]", [MoDoF("D2", D2TElem, 1, Fcn.cst(1), d0, "EntIdx", 1, "GInt", GInt("D2L", 4)), ...
        MoDoF("D2", D2TElem, 1, Fcn("D2R1", "s-1/2"), d0, "EntIdx", 1, "GInt", GInt("D2L", 4)), ...
        MoDoF("D2", D2TElem, 2, Fcn.cst(1), d0, "GInt", GInt("D2T", 2))]));
    small.SDG_1S = calcFES(mshS, SDG_1S_FES, u, d0, grad);
    SDG_1V_FES = FES(mshS, FE("D2T", FE.repFS("[1,x,y]", [2, 1]), ...
        [MoDoF("D2", D2TElem, 1, Fcn.cst(1), d0V, "coef", UNV, "EntIdx", [2, 3], "orien", true, "GInt", GInt("D2L", 4)), ...
        MoDoF("D2", D2TElem, 1, Fcn("D2R1", "s - 1/2"), d0V, "coef", UNV, "EntIdx", [2, 3], "orien", true, "GInt", GInt("D2L", 4)), ...
        MoDoF("D2", D2TElem, 2, [Fcn.cst([1; 0]), Fcn.cst([0; 1])], d0V, "GInt", GInt("D2T", 2))]));
    small.SDG_1V = calcFES(mshS, SDG_1V_FES, w, d0V, divV, divForm);
    trls = [SDG_1V_FES, SDG_1S_FES];
    Auv = [DLF(mshS, 2, Fcn.cst(1), d0V, d0V, "iTrl", 1, "iTst", 1, "GInt", GInt("D2T", 2)), ...
        DLF(mshS, 2, Fcn.cst(1), d0, divV, "iTrl", 2, "iTst", 1, "GInt", GInt("D2T", 1)), ...
        DLF.interface(mshS, 1, -UNV, d0, d0V, "EntIdx", PrEdge, "iTrl", 2, "iTst", 1, "tstOpr", "jump", "GInt", GInt("D2L", 2)), ...
        DLF(mshS, 2, Fcn.cst(1), d0V, grad, "iTrl", 1, "iTst", 2, "GInt", GInt("D2T", 1)), ...
        DLF.interface(mshS, 1, -UNV, d0V, d0, "EntIdx", DlEdge, "iTrl", 1, "iTst", 2, "tstOpr", "jump", "GInt", GInt("D2L", 2))];
    [small.SDG.Stiff, small.SDG.Load] = assemble(mshS, trls, trls, Auv, SLF(mshS, 2, u, d0, "iTst", 2, "GInt", GInt("D2T", 2)));
end
function dat = calcFES(msh, fES, fcn, d0, d1, d1Form)
    % calcFES: regression data of FE space on element (mesh, DoF map, mass matrix, projection and errors).
    arguments
        msh Msh;
        fES FES;
        fcn Fcn;
        d0; % Order of derivative for L2 norm.
        d1; % Order of derivative for second norm.
        d1Form = @(coef, fcn, pow) sum(coef .* abs(fcn) .^ pow); % Form of second norm.
    end
    dat.nGlDoF = fES.nGlDoF;
    dat.Lc2Gl = fES.Lc2Gl;
    if ~isempty(fES.BC)
        dat.BCIdx = fES.BC.DoFIdx;
        dat.BCVal = fES.BC.DoFVal;
    end
    dat.Mass = assemble(msh, fES, fES, DLF(msh, 2, Fcn.cst(1), d0, d0, "GInt", GInt("D2T", 4)), SLF.empty, "impBC", false);
    dat.Proj = projDoF(fES, fcn);
    fEF = fES.proj(fcn);
    dat.Err = [eNorm(msh, fcn, fEF, Norm(msh, 2, d0, "GInt", GInt("D2T", 4))), ...
        eNorm(msh, fcn, fEF, Norm(msh, 2, d1, "form", d1Form, "GInt", GInt("D2T", 4)))];
end
function DoFVal = projDoF(fES, fcn)
    % projDoF: global DoF values of projection (same as FES.proj).
    DoFVal = zeros(fES.nGlDoF, 1);
    for iDoF = 1:length(fES.GlDoFs)
        val = fES.GlDoFs(iDoF).eval(fcn, "valType", "num");
        DoFVal(fES.GlDoFs.sub2ind(iDoF)) = val(:);
    end
end
function example = calcExample(exName)
    % calcExample: error norms printed by examples.
    exDir = fullfile(fileparts(fileparts(fileparts(mfilename("fullpath")))), "FEM_Example");
    example = struct();
    for iEx = 1:length(exName)
        example.(exName(iEx)) = runExample(exDir, exName(iEx));
    end
end
function err = runExample(exDir, exName)
    % runExample: run example script and extract printed error norms, e.g. "|u-uh|_L2: 1.234567e-03".
    oldDir = cd(exDir);
    cleanup = onCleanup(@() cd(oldDir));
    exOut = evalc("runScript(fullfile(exDir, exName + "".m""))");
    close all;
    tok = regexp(exOut, "\|[^|]+\|_\w+:\s*([-+.\deEinfINFnaN]+)", "tokens");
    err = cellfun(@(t) str2double(t{1}), tok);
    assert(~isempty(err), "No error norm found in output of " + exName + ".");
end
function runScript(exFile)
    % runScript: run script in a separate workspace.
    run(exFile);
end
