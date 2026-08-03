function [Stiff, Load] = assemble(msh, trls, tsts, Auvs, Fvs, options)
    % assemble: assemble stiffness matrix and load vector.
    arguments (Input)
        msh Msh; % Mesh.
        trls (1, :) FES; % Trial function spaces.
        tsts (1, :) FES; % Test fun spaces.
        Auvs (1, :) DLF; % Bilinear forms a(u, v).
        Fvs (1, :) SLF; % Linear forms f(v).
        options.matType {mustBeMember(options.matType, ["sparse", "full"])} = "sparse"; % Type of matrix.
        options.impBC logical = true; % Whether impose boundary conditions.
        options.preSol (1, :) FEF = FEF.empty; % Previous solution.
        options.Awuv (1, :) LDLF = LDLF.empty; % Linearized bilinear forms.
        options.Fwv (1, :) LSLF = LSLF.empty; % Linearized linear forms.
    end
    arguments (Output)
        Stiff (:, :); % Stiffness matrix.
        Load (:, 1); % Load vector.
    end
    assert(ismember(msh.type, "D2T"));
    if ~isempty(trls)
        assert(isequal(msh, trls.getMsh));
    end
    if ~isempty(tsts)
        assert(isequal(msh, tsts.getMsh));
    end
    if ~isempty(Auvs)
        assert(isequal(msh, Auvs.getMsh));
    end
    if ~isempty(Fvs)
        assert(isequal(msh, Fvs.getMsh));
    end
    if ~isempty(options.preSol)
        assert(isequal(msh, options.preSol.getMsh));
    end
    if ~isempty(options.Awuv)
        assert(isequal(msh, options.Awuv.getMsh));
    end
    if ~isempty(options.Fwv)
        assert(isequal(msh, options.Fwv.getMsh));
    end
    switch options.matType
        case "sparse"
            nImax = 0;
            for iAuv = 1:length(Auvs)
                nImax = nImax + numel(Auvs(iAuv).EntIdx) * trls(abs(Auvs(iAuv).iTrl)).nLcDoF * tsts(abs(Auvs(iAuv).iTst)).nLcDoF;
            end
            for iAwuv = 1:length(options.Awuv)
                nImax = nImax + numel(options.Awuv(iAwuv).EntIdx) * trls(abs(options.Awuv(iAwuv).iTrl)).nLcDoF * tsts(abs(options.Awuv(iAwuv).iTst)).nLcDoF;
            end
            Is = zeros(nImax, 1); Js = zeros(nImax, 1); Vs = zeros(nImax, 1);
            nI = 0;
        case "full"
            Stiff = zeros(tsts.cumDoF, trls.cumDoF);
    end
    Load = zeros(tsts.cumDoF, 1);
    cumTrlDoF = zeros(1, length(trls) + 1);
    for iTrl = 1:length(trls) + 1
        cumTrlDoF(iTrl) = trls.cumDoF(iTrl - 1);
    end
    cumTstDoF = zeros(1, length(tsts) + 1);
    for iTst = 1:length(tsts) + 1
        cumTstDoF(iTst) = tsts.cumDoF(iTst - 1);
    end
    for iAuv = 1:length(Auvs)
        Auv = Auvs(iAuv);
        trl = trls(abs(Auv.iTrl));
        tst = tsts(abs(Auv.iTst));
        switch Auv.EntDim
            case 2
                assert(ismember(trl.elem, "D2T"));
                assert(ismember(tst.elem, "D2T"));
                trlBaseDiv = trl.LcBase.dif(Auv.trlOrd);
                tstBaseDiv = tst.LcBase.dif(Auv.tstOrd);
                for iTrlBase = 1:trl.nLcDoF
                    for iTstBase = 1:tst.nLcDoF
                        intFcn = Auv.form(Auv.coef, trlBaseDiv(iTrlBase), tstBaseDiv(iTstBase));
                        intFun = intFcn.getFun;
                        for iElem = Auv.EntIdx
                            ElParm = msh.node.coord(:, msh.elem.node(:, iElem));
                            intVal = Auv.GInt.eval(@(x) intFun(x, ElParm), ElParm);
                            I = cumTstDoF(Auv.iTst) + abs(tst.Lc2Gl(iTstBase, iElem));
                            J = cumTrlDoF(Auv.iTrl) + abs(trl.Lc2Gl(iTrlBase, iElem));
                            V = intVal * sign(tst.Lc2Gl(iTstBase, iElem)) * sign(trl.Lc2Gl(iTrlBase, iElem));
                            asmStiff(I, J, V);
                        end
                    end
                end
            case 1
                assert(ismember(trl.elem, ["D2T", "D2LR"]));
                assert(ismember(tst.elem, ["D2T", "D2LR"]));
                EgPmSym = MshEnt("D2LR").parm;
                TrlPmSym = sym('trlParm', trl.sElParm);
                TstPmSym = sym('tstParm', tst.sElParm);
                coef = Auv.coef.tfm("D2LR");
                JNorm = Tfm("D2L").JNorm;
                switch trl.elem
                    case "D2T"
                        trlBaseDiv = trl.LcBase.dif(Auv.trlOrd).subParm(TrlPmSym).tfm("D2LR");
                    case "D2LR"
                        assert(all(Auv.trlOrd == 0));
                        trlBaseDiv = trl.LcBase.dif(Auv.trlOrd);
                end
                switch tst.elem
                    case "D2T"
                        tstBaseDiv = tst.LcBase.dif(Auv.tstOrd).subParm(TstPmSym).tfm("D2LR");
                    case "D2LR"
                        assert(all(Auv.tstOrd == 0));
                        tstBaseDiv = tst.LcBase.dif(Auv.tstOrd);
                end
                for iTrlBase = 1:trl.nLcDoF
                    for iTstBase = 1:tst.nLcDoF
                        intFcn = Auv.form(coef, trlBaseDiv(iTrlBase), tstBaseDiv(iTstBase)) .* JNorm;
                        intFun = intFcn.getFun("parm", {EgPmSym, TrlPmSym, TstPmSym});
                        for iEdge = Auv.EntIdx
                            EgParm = msh.node.coord(:, msh.edge.node(:, iEdge));
                            CnElem = msh.edge.elem(:, iEdge);
                            switch trl.elem
                                case "D2T"
                                    iTrlEnt = abs(CnElem(sign(CnElem) == sign(Auv.iTrl)));
                                case "D2LR"
                                    iTrlEnt = iEdge;
                            end
                            switch tst.elem
                                case "D2T"
                                    iTstEnt = abs(CnElem(sign(CnElem) == sign(Auv.iTst)));
                                case "D2LR"
                                    iTstEnt = iEdge;
                            end
                            if isempty(iTrlEnt) || isempty(iTstEnt)
                                continue
                            end
                            trlParm = trl.ElParm(:, :, iTrlEnt);
                            tstParm = tst.ElParm(:, :, iTstEnt);
                            intVal = Auv.GInt.eval(@(x) intFun(x, EgParm, trlParm, tstParm));
                            I = cumTstDoF(abs(Auv.iTst)) + abs(tst.Lc2Gl(iTstBase, iTstEnt));
                            J = cumTrlDoF(abs(Auv.iTrl)) + abs(trl.Lc2Gl(iTrlBase, iTrlEnt));
                            V = intVal * sign(tst.Lc2Gl(iTstBase, iTstEnt)) * sign(trl.Lc2Gl(iTrlBase, iTrlEnt));
                            asmStiff(I, J, V);
                        end
                    end
                end
        end
    end
    for iAwuv = 1:length(options.Awuv)
        Awuv = options.Awuv(iAwuv);
        preSol = options.preSol(abs(Awuv.iPre));
        trl = trls(abs(Awuv.iTrl));
        tst = tsts(abs(Awuv.iTst));
        switch Awuv.EntDim
            case 2
                assert(ismember(preSol.elem, "D2T"));
                assert(ismember(trl.elem, "D2T"));
                assert(ismember(tst.elem, "D2T"));
                preSolDiv = preSol.dif(Awuv.preOrd);
                trlBaseDiv = trl.LcBase.dif(Awuv.trlOrd);
                tstBaseDiv = tst.LcBase.dif(Awuv.tstOrd);
                for iTrlBase = 1:trl.nLcDoF
                    for iTstBase = 1:tst.nLcDoF
                        intFcn = Awuv.form(Awuv.coef, preSolDiv, trlBaseDiv(iTrlBase), tstBaseDiv(iTstBase));
                        intFun = intFcn.getFun;
                        for iElem = Awuv.EntIdx
                            ElParm = msh.node.coord(:, msh.elem.node(:, iElem));
                            preCoef = preSol.ElCoef(:, iElem);
                            intVal = Awuv.GInt.eval(@(x) intFun(x, ElParm, preCoef), ElParm);
                            I = cumTstDoF(Awuv.iTst) + abs(tst.Lc2Gl(iTstBase, iElem));
                            J = cumTrlDoF(Awuv.iTrl) + abs(trl.Lc2Gl(iTrlBase, iElem));
                            V = intVal * sign(tst.Lc2Gl(iTstBase, iElem)) * sign(trl.Lc2Gl(iTrlBase, iElem));
                            asmStiff(I, J, V);
                        end
                    end
                end
            % case 1 % TODO
        end
    end
    if isequal(options.matType, "sparse")
        Stiff = sparse(Is(1:nI), Js(1:nI), Vs(1:nI), tsts.cumDoF, trls.cumDoF);
    end
    for iFv = 1:length(Fvs)
        Fv = Fvs(iFv);
        tst = tsts(Fv.iTst);
        switch Fv.EntDim
            case 2
                assert(isequal(tst.elem, "D2T"));
                tstBaseDiv = tst.LcBase.dif(Fv.tstOrd);
                for iTstBase = 1:tst.nLcDoF
                    intFcn = Fv.form(Fv.load, tstBaseDiv(iTstBase));
                    intFun = intFcn.getFun;
                    for iElem = Fv.EntIdx
                        ElParm = msh.node.coord(:, msh.elem.node(:, iElem));
                        intVal = Fv.GInt.eval(@(x) intFun(x, ElParm), ElParm);
                        I = cumTstDoF(Fv.iTst) + abs(tst.Lc2Gl(iTstBase, iElem));
                        V = intVal * sign(tst.Lc2Gl(iTstBase, iElem));
                        asmLoad(I, V);
                    end
                end
            case 1
                assert(ismember(tst.elem, ["D2T", "D2LR"]));
                EgPmSym = MshEnt("D2LR").parm;
                TstPmSym = sym('tstParm', tst.sElParm);
                load = Fv.load.tfm("D2LR");
                JNorm = Tfm("D2L").JNorm;
                switch tst.elem
                    case "D2T"
                        tstBaseDiv = tst.LcBase.dif(Fv.tstOrd).subParm(TstPmSym).tfm("D2LR");
                    case "D2LR"
                        tstBaseDiv = tst.LcBase.dif(Fv.tstOrd);
                end
                for iTstBase = 1:tst.nLcDoF
                    intFcn = Fv.form(load, tstBaseDiv(iTstBase)) .* JNorm;
                    intFun = intFcn.getFun("parm", {EgPmSym, TstPmSym});
                    for iEdge = Fv.EntIdx
                        EgParm = msh.node.coord(:, msh.edge.node(:, iEdge));
                        CnElem = msh.edge.elem(:, iEdge);
                        switch tst.elem
                            case "D2T"
                                iTstEnt = abs(CnElem(sign(CnElem) == sign(Fv.iTst)));
                            case "D2LR"
                                iTstEnt = iEdge;
                        end
                        if isempty(iTstEnt)
                            continue
                        end
                        tstParm = tst.ElParm(:, :, iTstEnt);
                        intVal = Fv.GInt.eval(@(x) intFun(x, EgParm, tstParm));
                        I = cumTstDoF(Fv.iTst) + abs(tst.Lc2Gl(iTstBase, iTstEnt));
                        V = intVal * sign(tst.Lc2Gl(iTstBase, iTstEnt));
                        asmLoad(I, V);
                    end
                end
        end
    end
    for iFwv = 1:length(options.Fwv)
        Fwv = options.Fwv(iFwv);
        preSol = options.preSol(abs(Fwv.iPre));
        tst = tsts(Fwv.iTst);
        switch Fwv.EntDim
            case 2
                assert(isequal(preSol.elem, "D2T"));
                assert(isequal(tst.elem, "D2T"));
                preSolDiv = preSol.dif(Fwv.preOrd);
                tstBaseDiv = tst.LcBase.dif(Fwv.tstOrd);
                for iTstBase = 1:tst.nLcDoF
                    intFcn = Fwv.form(Fwv.load, preSolDiv, tstBaseDiv(iTstBase));
                    intFun = intFcn.getFun;
                    for iElem = Fwv.EntIdx
                        ElParm = msh.node.coord(:, msh.elem.node(:, iElem));
                        preCoef = preSol.ElCoef(:, iElem);
                        intVal = Fwv.GInt.eval(@(x) intFun(x, ElParm, preCoef), ElParm);
                        I = cumTstDoF(Fwv.iTst) + abs(tst.Lc2Gl(iTstBase, iElem));
                        V = intVal * sign(tst.Lc2Gl(iTstBase, iElem));
                        asmLoad(I, V);
                    end
                end
            % case 1 % TODO
        end
    end
    if options.impBC
        for iTrl = 1:length(trls)
            if ~isempty(trls(iTrl).BC)
                BdDoFIdx = trls.cumDoF(iTrl - 1) + trls(iTrl).BC.DoFIdx;
                Stiff(BdDoFIdx, :) = 0;
                Stiff(sub2ind(size(Stiff), BdDoFIdx, BdDoFIdx)) = 1;
                Load(BdDoFIdx) = trls(iTrl).BC.DoFVal;
            end
        end
    end
    % Local functions.
    function asmStiff(I, J, V)
        % asmStiff: add value V to entry (I, J) of the stiffness matrix.
        switch options.matType
            case "sparse"
                nI = nI + 1;
                if nI > nImax
                    nImax = floor(nImax * 1.5);
                    Is = [Is; zeros(nImax - length(Is), 1)];
                    Js = [Js; zeros(nImax - length(Js), 1)];
                    Vs = [Vs; zeros(nImax - length(Vs), 1)];
                end
                Is(nI) = I; Js(nI) = J; Vs(nI) = V;
            case "full"
                Stiff(I, J) = Stiff(I, J) + V;
        end
    end
    function asmLoad(I, V)
        % asmLoad: add value V to entry I of the load vector.
        Load(I) = Load(I) + V;
    end
end
