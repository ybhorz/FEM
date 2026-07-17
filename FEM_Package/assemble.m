function [Stiff, Load] = assemble(msh, trls, tsts, Auvs, Fvs, options)
    % assemble: assemble stiffness matrix and load vector.
    arguments (Input)
        msh Msh; % Mesh.
        trls (1, :) FES; % Trial finite element spaces.
        tsts (1, :) FES; % Test finite element spaces.
        Auvs (1, :) DLF; % Double linear functional.
        Fvs (1, :) SLF; % Single linear functional.
        options.matType {mustBeMember(options.matType, ["sparse", "full"])} = "sparse"; % Type of matrix.
        options.impBC logical = true; % Whether impose boundary conditions.
        options.preSol (1, :) FEF = FEF.empty; % Previous solution.
        options.Awuv (1, :) LDLF = LDLF.empty; % Linearized double linear functional.
        options.Fwv (1, :) LSLF = LSLF.empty; % Linearized single linear functional.
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
            nImax = tsts.cumDoF * sum([trls.nLcDoF]) * 10;
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
                    end
                end
            % case 1 % TODO: implement this case.
            %     assert(ismember(preSol.elem, "D2T"));
            %     assert(ismember(trl.elem, ["D2T", "D2LR"]));
            %     assert(ismember(tst.elem, ["D2T", "D2LR"]));
            %     EgPmSym = MshEnt("D2LR").parm;
            %     PrePmSym = sym('preParm', preSol.sElParm);
            %     TrlPmSym = sym('trlParm', trl.sElParm);
            %     TstPmSym = sym('tstParm', tst.sElParm);
            %     coef = Awuv.coef.tfm("D2LR");
            %     JNorm = Tfm("D2L").JNorm;
            %     switch preSol.elem
            %         case "D2T"
            %             preSolDiv = preSol.dif(Awuv.preOrd).subParm(PrePmSym).tfm("D2LR");
            %     end
            %     switch trl.elem
            %         case "D2T"
            %             trlBaseDiv = trl.LcBase.dif(Awuv.trlOrd).subParm(TrlPmSym).tfm("D2LR");
            %         case "D2LR"
            %             assert(all(Awuv.trlOrd == 0));
            %             trlBaseDiv = trl.LcBase.dif(Awuv.trlOrd);
            %     end
            %     switch tst.elem
            %         case "D2T"
            %             tstBaseDiv = tst.LcBase.dif(Awuv.tstOrd).subParm(TstPmSym).tfm("D2LR");
            %         case "D2LR"
            %             assert(all(Awuv.tstOrd == 0));
            %             tstBaseDiv = tst.LcBase.dif(Awuv.tstOrd);
            %     end
            %     for iTrlBase = 1:trl.nLcDoF
            %         for iTstBase = 1:tst.nLcDoF
            %             intFcn = Awuv.form(coef, preSolDiv, trlBaseDiv(iTrlBase), tstBaseDiv(iTstBase)) .* JNorm;
            %             intFun = intFcn.getFun("parm", {EgPmSym, PrePmSym, TrlPmSym, TstPmSym});
            %             for iEdge = Awuv.EntIdx
            %                 EgParm = msh.node.coord(:, msh.edge.node(:, iEdge));
            %                 CnElem = msh.edge.elem(:, iEdge);
            %                 switch preSol.elem
            %                     case "D2T"
            %                         iPreEnt = abs(CnElem(sign(CnElem) == sign(Awuv.iPre)));
            %                 end
            %                 switch trl.elem
            %                     case "D2T"
            %                         iTrlEnt = abs(CnElem(sign(CnElem) == sign(Awuv.iTrl)));
            %                     case "D2LR"
            %                         iTrlEnt = iEdge;
            %                 end
            %                 switch tst.elem
            %                     case "D2T"
            %                         iTstEnt = abs(CnElem(sign(CnElem) == sign(Awuv.iTst)));
            %                     case "D2LR"
            %                         iTstEnt = iEdge;
            %                 end
            %                 if isempty(iPreEnt) || isempty(iTrlEnt) || isempty(iTstEnt)
            %                     continue
            %                 end
            %                 preParm = preSol.ElParm(:, :, iPreEnt);
            %                 trlParm = trl.ElParm(:, :, iTrlEnt);
            %                 tstParm = tst.ElParm(:, :, iTstEnt);
            %                 preCoef = preSol.ElCoef(:, iPreEnt);
            %                 intVal = Awuv.GInt.eval(@(x) intFun(x, EgParm, preParm, trlParm, tstParm, preCoef));
            %                 I = cumTstDoF(abs(Awuv.iTst)) + abs(tst.Lc2Gl(iTstBase, iTstEnt));
            %                 J = cumTrlDoF(abs(Awuv.iTrl)) + abs(trl.Lc2Gl(iTrlBase, iTrlEnt));
            %                 V = intVal * sign(tst.Lc2Gl(iTstBase, iTstEnt)) * sign(trl.Lc2Gl(iTrlBase, iTrlEnt));
            %                 switch options.matType
            %                     case "sparse"
            %                         nI = nI + 1;
            %                         if nI > nImax
            %                             nImax = floor(nImax * 1.5);
            %                             Is = [Is; zeros(nImax - length(Is), 1)];
            %                             Js = [Js; zeros(nImax - length(Js), 1)];
            %                             Vs = [Vs; zeros(nImax - length(Vs), 1)];
            %                         end
            %                         Is(nI) = I; Js(nI) = J; Vs(nI) = V;
            %                     case "full"
            %                         Stiff(I, J) = Stiff(I, J) + V;
            %                 end
            %             end
            %         end
            %     end
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
                        Load(I) = Load(I) + V;
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
                        Load(I) = Load(I) + V;
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
                        Load(I) = Load(I) + V;
                    end
                end
            % case 1 % TODO: implement this case.
            %     assert(ismember(tst.elem, ["D2T", "D2LR"]));
            %     EgPmSym = MshEnt("D2LR").parm;
            %     TstPmSym = sym('tstParm', tst.sElParm);
            %     load = Fv.load.tfm("D2LR");
            %     JNorm = Tfm("D2L").JNorm;
            %     switch tst.elem
            %         case "D2T"
            %             tstBaseDiv = tst.LcBase.dif(Fv.tstOrd).subParm(TstPmSym).tfm("D2LR");
            %         case "D2LR"
            %             tstBaseDiv = tst.LcBase.dif(Fv.tstOrd);
            %     end
            %     for iTstBase = 1:tst.nLcDoF
            %         intFcn = Fv.form(load, tstBaseDiv(iTstBase)) .* JNorm;
            %         intFun = intFcn.getFun("parm", {EgPmSym, TstPmSym});
            %         for iEdge = Fv.EntIdx
            %             EgParm = msh.node.coord(:, msh.edge.node(:, iEdge));
            %             CnElem = msh.edge.elem(:, iEdge);
            %             switch tst.elem
            %                 case "D2T"
            %                     iTstEnt = abs(CnElem(sign(CnElem) == sign(Fv.iTst)));
            %                 case "D2LR"
            %                     iTstEnt = iEdge;
            %             end
            %             if isempty(iTstEnt)
            %                 continue
            %             end
            %             tstParm = tst.ElParm(:, :, iTstEnt);
            %             intVal = Fv.GInt.eval(@(x) intFun(x, EgParm, tstParm));
            %             I = cumTstDoF(Fv.iTst) + abs(tst.Lc2Gl(iTstBase, iTstEnt));
            %             V = intVal * sign(tst.Lc2Gl(iTstBase, iTstEnt));
            %             Load(I) = Load(I) + V;
            %         end
            %     end
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
end
