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
    assert(ismember(msh.type, ["D2T", "D3T"]));
    ElDomn = msh.ElDomn; % Domain of element.
    FtDomn = msh.FtDomn; % Domain of facet.
    FtRefDomn = msh.FtRefDomn; % Reference domain of facet.
    facet = msh.facet; % Facet.
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
    trcCache = containers.Map(); % Cache of derivatives of base functions (`difBase`, `trcBase`).
    cumTstDoF = zeros(1, length(tsts) + 1);
    for iTst = 1:length(tsts) + 1
        cumTstDoF(iTst) = tsts.cumDoF(iTst - 1);
    end
    for iAuv = 1:length(Auvs)
        asmDLF(Auvs(iAuv), FEF.empty);
    end
    for iAwuv = 1:length(options.Awuv)
        asmDLF(options.Awuv(iAwuv), options.preSol(abs(options.Awuv(iAwuv).iPre)));
    end
    if isequal(options.matType, "sparse")
        Stiff = sparse(Is(1:nI), Js(1:nI), Vs(1:nI), tsts.cumDoF, trls.cumDoF);
    end
    for iFv = 1:length(Fvs)
        asmSLF(Fvs(iFv), FEF.empty);
    end
    for iFwv = 1:length(options.Fwv)
        asmSLF(options.Fwv(iFwv), options.preSol(abs(options.Fwv(iFwv).iPre)));
    end
    if options.impBC
        % Rows of boundary DoFs are replaced by identity rows (by a diagonal mask, since assigning rows of a sparse
        % matrix is slow).
        isBd = false(size(Stiff, 1), 1);
        for iTrl = 1:length(trls)
            if ~isempty(trls(iTrl).BC)
                BdDoFIdx = trls.cumDoF(iTrl - 1) + trls(iTrl).BC.DoFIdx;
                isBd(BdDoFIdx) = true;
                Load(BdDoFIdx) = trls(iTrl).BC.DoFVal;
            end
        end
        if any(isBd)
            nRow = size(Stiff, 1);
            BdIdx = find(isBd);
            if issparse(Stiff)
                Stiff = spdiags(double(~isBd), 0, nRow, nRow) * Stiff + sparse(BdIdx, BdIdx, 1, nRow, size(Stiff, 2));
            else
                Stiff(BdIdx, :) = 0;
                Stiff(sub2ind(size(Stiff), BdIdx, BdIdx)) = 1;
            end
        end
    end
    % Nested functions.
    function asmDLF(Auv, preSol)
        % asmDLF: assemble double linear functional into stiffness matrix.
        % If `preSol` is non-empty, `Auv` is a linearized double linear functional (LDLF) with previous solution `preSol`.
        % On facet, sign of `iTrl`, `iTst` and `iPre` selects the connected element providing the trace.
        % Base functions (and previous solution) are evaluated numerically at Gauss points of all mesh entities; the
        % form is turned into a function handle once, with placeholders `trlVal`, `tstVal` (and `preVal`) for their values.
        trl = trls(abs(Auv.iTrl));
        tst = tsts(abs(Auv.iTst));
        isLin = ~isempty(preSol);
        trlName = sprintf("trl%d", abs(Auv.iTrl));
        tstName = sprintf("tst%d", abs(Auv.iTst));
        switch Auv.EntDim
            case msh.dim
                assert(ismember(trl.elem, ElDomn));
                assert(ismember(tst.elem, ElDomn));
                ElIdx = Auv.EntIdx;
                nEnt = length(ElIdx);
                if nEnt == 0
                    return
                end
                ElParm = reshape(msh.node.coord(:, msh.elem.node(:, ElIdx)), msh.dim, msh.elem.nNode, nEnt);
                ElPmArg = reshape(ElParm, [], 1, nEnt);
                [X, Jac] = Auv.GInt.pntVec(ElParm);
                [TrlVal, trlSz, trlDomn] = elBase(trl, Auv.trlOrd, trlName, X, ElParm);
                [TstVal, tstSz, tstDomn] = elBase(tst, Auv.tstOrd, tstName, X, ElParm);
                ArgSym = {MshEnt(ElDomn).parm};
                Args = {ElPmArg};
                if isLin
                    assert(ismember(preSol.elem, ElDomn));
                    preSolDiv = preSol.dif(Auv.preOrd);
                    PreCfArg = reshape(preSol.ElCoef(:, ElIdx), [], 1, nEnt);
                    PreVal = evalBase(preSolDiv, sprintf("el_pre%d", abs(Auv.iPre)), Auv.preOrd, X, {}, {ElPmArg}, {}, {PreCfArg});
                end
                coef = Auv.coef;
                JNorm = Fcn.cst(1);
                trlL2G = trl.Lc2Gl(:, ElIdx);
                tstL2G = tst.Lc2Gl(:, ElIdx);
            case msh.dim - 1
                assert(ismember(trl.elem, [ElDomn, FtRefDomn]));
                assert(ismember(tst.elem, [ElDomn, FtRefDomn]));
                FtPmSym = MshEnt(FtRefDomn).parm;
                TrlPmSym = sym('trlParm', trl.sElParm);
                TstPmSym = sym('tstParm', tst.sElParm);
                % Facets with connected entities on the required sides.
                FtIdx = Auv.EntIdx;
                CnElem = facet.elem(:, FtIdx);
                iTrlEnt = cnEnt(trl, CnElem, FtIdx, Auv.iTrl);
                iTstEnt = cnEnt(tst, CnElem, FtIdx, Auv.iTst);
                isVal = iTrlEnt > 0 & iTstEnt > 0;
                if isLin
                    iPreEnt = cnEnt(preSol, CnElem, FtIdx, Auv.iPre);
                    isVal = isVal & iPreEnt > 0;
                    iPreEnt = iPreEnt(isVal);
                end
                FtIdx = FtIdx(isVal); iTrlEnt = iTrlEnt(isVal); iTstEnt = iTstEnt(isVal);
                nEnt = length(FtIdx);
                if nEnt == 0
                    return
                end
                FtPmArg = reshape(msh.node.coord(:, facet.node(:, FtIdx)), [], 1, nEnt);
                TrlPmArg = reshape(trl.ElParm(:, :, iTrlEnt), [], 1, nEnt);
                TstPmArg = reshape(tst.ElParm(:, :, iTstEnt), [], 1, nEnt);
                [X, Jac] = Auv.GInt.pntVec([], nEnt);
                [TrlVal, trlSz, trlDomn] = ftBase(trl, Auv.trlOrd, trlName, X, FtPmArg, TrlPmSym, TrlPmArg, iTrlEnt);
                [TstVal, tstSz, tstDomn] = ftBase(tst, Auv.tstOrd, tstName, X, FtPmArg, TstPmSym, TstPmArg, iTstEnt);
                ArgSym = {FtPmSym};
                Args = {FtPmArg};
                if isLin
                    assert(ismember(preSol.elem, [ElDomn, FtRefDomn]));
                    PrePmSym = sym('preParm', preSol.sElParm);
                    PreCfSym = sym('preCoef', [preSol.nElCoef, 1]);
                    preSolDiv = trcPre(preSol, Auv.preOrd, PrePmSym, PreCfSym);
                    PrePmArg = reshape(preSol.ElParm(:, :, iPreEnt), [], 1, nEnt);
                    PreCfArg = reshape(preSol.ElCoef(:, iPreEnt), [], 1, nEnt);
                    PreVal = evalBase(preSolDiv, sprintf("ft_pre%d", abs(Auv.iPre)), Auv.preOrd, X, {FtPmSym, PrePmSym}, ...
                        {FtPmArg, PrePmArg}, {PreCfSym}, {PreCfArg});
                end
                coef = Auv.coef.tfm(FtRefDomn);
                JNorm = Tfm(FtDomn).JNorm;
                trlL2G = trl.Lc2Gl(:, iTrlEnt);
                tstL2G = tst.Lc2Gl(:, iTstEnt);
        end
        % Form with placeholders for values of base functions (and previous solution).
        TrlSym = sym('trlVal', trlSz);
        TstSym = sym('tstVal', tstSz);
        trlPh = Fcn(trlDomn, TrlSym);
        tstPh = Fcn(tstDomn, TstSym);
        % Coefficients and Jacobian factor are evaluated once at Gauss points and enter the form as placeholders too,
        % so that the integrand handle (called for every pair of base functions) does not repeat their expressions,
        % e.g. unit normal and area factor of facets in terms of facet vertices.
        [CfPh, CfSym, CfVal] = phCoef([coef, JNorm], X, ArgSym, Args, nEnt);
        if isLin
            PreSym = sym('preVal', size(preSolDiv(1).fun));
            intFcn = Auv.form(CfPh(1:end - 1), Fcn(trlDomn, PreSym), trlPh, tstPh) .* CfPh(end);
            intFun = intFcn.getFun("parm", ArgSym, "coef", [CfSym, {PreSym, TrlSym, TstSym}], "vec", true);
        else
            intFcn = Auv.form(CfPh(1:end - 1), trlPh, tstPh) .* CfPh(end);
            intFun = intFcn.getFun("parm", ArgSym, "coef", [CfSym, {TrlSym, TstSym}], "vec", true);
        end
        [iTstBs, iTrlBs] = ndgrid(1:tst.nLcDoF, 1:trl.nLcDoF);
        nPair = numel(iTstBs);
        nPnt = size(X, 2);
        nTst = tst.nLcDoF;
        intVal = zeros(nPair, nEnt);
        % Test base functions are batched: their values are stacked along entities (chunks of bounded size).
        for iTsts = chunkIdx(nTst, numel(TstVal) / nTst)
            m = numel(iTsts{1});
            S = reshape(permute(TstVal(:, iTsts{1}, :, :), [1, 3, 4, 2]), [], nPnt, nEnt * m);
            [XRep, ArgsRep, JacRep] = repEnt(X, Args, Jac, m);
            CfRep = cellfun(@(val) repmat(val, 1, 1, m), CfVal, "UniformOutput", false);
            if isLin
                PreRep = repmat(reshape(PreVal, [], nPnt, nEnt), 1, 1, m);
            end
            for jTrl = 1:trl.nLcDoF
                T = repmat(reshape(TrlVal(:, jTrl, :, :), [], nPnt, nEnt), 1, 1, m);
                if isLin
                    F = intFun(XRep, ArgsRep{:}, CfRep{:}, PreRep, T, S);
                else
                    F = intFun(XRep, ArgsRep{:}, CfRep{:}, T, S);
                end
                intVal(iTsts{1} + nTst * (jTrl - 1), :) = reshape(Auv.GInt.sumVec(F, JacRep), nEnt, m).';
            end
        end
        tstL2G = tstL2G(iTstBs(:), :);
        trlL2G = trlL2G(iTrlBs(:), :);
        asmStiff(cumTstDoF(abs(Auv.iTst)) + abs(tstL2G), cumTrlDoF(abs(Auv.iTrl)) + abs(trlL2G), ...
            intVal .* sign(tstL2G) .* sign(trlL2G));
    end
    function asmSLF(Fv, preSol)
        % asmSLF: assemble single linear functional into load vector.
        % If `preSol` is non-empty, `Fv` is a linearized single linear functional (LSLF) with previous solution `preSol`.
        % On facet, sign of `iTst` and `iPre` selects the connected element providing the trace.
        % Evaluation follows `asmDLF`.
        tst = tsts(abs(Fv.iTst));
        isLin = ~isempty(preSol);
        tstName = sprintf("tst%d", abs(Fv.iTst));
        switch Fv.EntDim
            case msh.dim
                assert(isequal(tst.elem, ElDomn));
                ElIdx = Fv.EntIdx;
                nEnt = length(ElIdx);
                if nEnt == 0
                    return
                end
                ElParm = reshape(msh.node.coord(:, msh.elem.node(:, ElIdx)), msh.dim, msh.elem.nNode, nEnt);
                ElPmArg = reshape(ElParm, [], 1, nEnt);
                [X, Jac] = Fv.GInt.pntVec(ElParm);
                [TstVal, tstSz, tstDomn] = elBase(tst, Fv.tstOrd, tstName, X, ElParm);
                ArgSym = {MshEnt(ElDomn).parm};
                Args = {ElPmArg};
                if isLin
                    assert(isequal(preSol.elem, ElDomn));
                    preSolDiv = preSol.dif(Fv.preOrd);
                    PreCfArg = reshape(preSol.ElCoef(:, ElIdx), [], 1, nEnt);
                    PreVal = evalBase(preSolDiv, sprintf("el_pre%d", abs(Fv.iPre)), Fv.preOrd, X, {}, {ElPmArg}, {}, {PreCfArg});
                end
                load = Fv.load;
                JNorm = Fcn.cst(1);
                tstL2G = tst.Lc2Gl(:, ElIdx);
            case msh.dim - 1
                assert(ismember(tst.elem, [ElDomn, FtRefDomn]));
                FtPmSym = MshEnt(FtRefDomn).parm;
                TstPmSym = sym('tstParm', tst.sElParm);
                FtIdx = Fv.EntIdx;
                CnElem = facet.elem(:, FtIdx);
                iTstEnt = cnEnt(tst, CnElem, FtIdx, Fv.iTst);
                isVal = iTstEnt > 0;
                if isLin
                    iPreEnt = cnEnt(preSol, CnElem, FtIdx, Fv.iPre);
                    isVal = isVal & iPreEnt > 0;
                    iPreEnt = iPreEnt(isVal);
                end
                FtIdx = FtIdx(isVal); iTstEnt = iTstEnt(isVal);
                nEnt = length(FtIdx);
                if nEnt == 0
                    return
                end
                FtPmArg = reshape(msh.node.coord(:, facet.node(:, FtIdx)), [], 1, nEnt);
                TstPmArg = reshape(tst.ElParm(:, :, iTstEnt), [], 1, nEnt);
                [X, Jac] = Fv.GInt.pntVec([], nEnt);
                [TstVal, tstSz, tstDomn] = ftBase(tst, Fv.tstOrd, tstName, X, FtPmArg, TstPmSym, TstPmArg, iTstEnt);
                ArgSym = {FtPmSym};
                Args = {FtPmArg};
                if isLin
                    assert(ismember(preSol.elem, [ElDomn, FtRefDomn]));
                    PrePmSym = sym('preParm', preSol.sElParm);
                    PreCfSym = sym('preCoef', [preSol.nElCoef, 1]);
                    preSolDiv = trcPre(preSol, Fv.preOrd, PrePmSym, PreCfSym);
                    PrePmArg = reshape(preSol.ElParm(:, :, iPreEnt), [], 1, nEnt);
                    PreCfArg = reshape(preSol.ElCoef(:, iPreEnt), [], 1, nEnt);
                    PreVal = evalBase(preSolDiv, sprintf("ft_pre%d", abs(Fv.iPre)), Fv.preOrd, X, {FtPmSym, PrePmSym}, ...
                        {FtPmArg, PrePmArg}, {PreCfSym}, {PreCfArg});
                end
                load = Fv.load.tfm(FtRefDomn);
                JNorm = Tfm(FtDomn).JNorm;
                tstL2G = tst.Lc2Gl(:, iTstEnt);
        end
        TstSym = sym('tstVal', tstSz);
        tstPh = Fcn(tstDomn, TstSym);
        if isLin
            PreSym = sym('preVal', size(preSolDiv(1).fun));
            intFcn = Fv.form(load, Fcn(tstDomn, PreSym), tstPh) .* JNorm;
            intFun = intFcn.getFun("parm", ArgSym, "coef", {PreSym, TstSym}, "vec", true);
        else
            intFcn = Fv.form(load, tstPh) .* JNorm;
            intFun = intFcn.getFun("parm", ArgSym, "coef", {TstSym}, "vec", true);
        end
        nPnt = size(X, 2);
        nTst = tst.nLcDoF;
        intVal = zeros(nTst, nEnt);
        for iTsts = chunkIdx(nTst, numel(TstVal) / nTst)
            m = numel(iTsts{1});
            S = reshape(permute(TstVal(:, iTsts{1}, :, :), [1, 3, 4, 2]), [], nPnt, nEnt * m);
            [XRep, ArgsRep, JacRep] = repEnt(X, Args, Jac, m);
            if isLin
                F = intFun(XRep, ArgsRep{:}, repmat(reshape(PreVal, [], nPnt, nEnt), 1, 1, m), S);
            else
                F = intFun(XRep, ArgsRep{:}, S);
            end
            intVal(iTsts{1}, :) = reshape(Fv.GInt.sumVec(F, JacRep), nEnt, m).';
        end
        asmLoad(cumTstDoF(abs(Fv.iTst)) + abs(tstL2G), intVal .* sign(tstL2G));
    end
    function val = evalBase(fcns, name, ord, X, ParmSym, ParmArg, CoefSym, CoefArg)
        % evalBase: values of functions `fcns` (of the same size, e.g. derivatives of base functions) at Gauss points X
        % of a batch of mesh entities: val(k, i, :, :) is k-th component of i-th function, nComp x nFcn x nPnt x nEnt.
        % `ParmSym`, `CoefSym`: symbolic parameters and coefficients of `fcns` (defaults of `Fcn.getFun` if empty);
        % `ParmArg`, `CoefArg`: their batched values. The function handle is cached by `name` and `ord`.
        key = char("fun_" + name + "_" + mat2str(size(ord)) + "_" + mat2str(ord(:)'));
        if isKey(trcCache, key)
            funH = trcCache(key);
        else
            nComp = numel(fcns(1).fun);
            comps(1:numel(fcns) * nComp) = Fcn.cst(0);
            for iFcn = 1:numel(fcns)
                for iComp = 1:nComp
                    comps((iFcn - 1) * nComp + iComp) = Fcn(fcns(iFcn).domn, fcns(iFcn).fun(iComp), fcns(iFcn).coef);
                end
            end
            funH = comps.getFun("parm", ParmSym, "coef", CoefSym, "vec", true);
            trcCache(key) = funH;
        end
        val = funH(X, ParmArg{:}, CoefArg{:});
        val = reshape(val, numel(fcns(1).fun), numel(fcns), size(val, 2), size(val, 3));
    end
    function [val, sDiv, domn] = elBase(fES, ord, name, X, ElParm)
        % elBase: values of derivative `ord` of local base functions of element space at points X of elements with
        % vertices ElParm (nComp x nBase x nPnt x nEnt), and size and domain of the derivative.
        % Mapped base functions (see `FE`) are evaluated numerically from reference base functions (`mapVal`).
        if isequal(fES.map, "none")
            baseDiv = difBase(fES, ord, name);
            val = evalBase(baseDiv, "el_" + name, ord, X, {}, {reshape(ElParm, [], 1, size(ElParm, 3))}, {}, {});
            sDiv = size(baseDiv(1).fun);
            domn = baseDiv(1).domn;
        else
            val = mapVal(fES.RefBase, fES.map, ord, X, ElParm, fES.RefKey);
            sDiv = Fcn.difSize(size(fES.LcBase(1).fun), ord);
            domn = ElDomn;
        end
    end
    function [val, sDiv, domn] = ftBase(fES, ord, name, X, FtPmArg, PmSym, PmArg, iEnt)
        % ftBase: values of derivative `ord` of local base functions at points X of reference facet, traces from
        % entities iEnt of FE space (nComp x nBase x nPnt x nEnt), and size and domain of the derivative.
        % `FtPmArg`: vertices of facets; `PmSym`, `PmArg`: symbolic and batched parameters of the entities.
        if isequal(fES.elem, ElDomn) && ~isequal(fES.map, "none")
            % Mapped base functions are evaluated at physical points of facets.
            nEnt = size(FtPmArg, 3);
            FtParm = reshape(FtPmArg, msh.dim, facet.nNode, nEnt);
            XPhy = pagemtimes(FtParm(:, 2:end, :) - FtParm(:, 1, :), X) + FtParm(:, 1, :);
            val = mapVal(fES.RefBase, fES.map, ord, XPhy, fES.ElParm(:, :, iEnt), fES.RefKey);
            sDiv = Fcn.difSize(size(fES.LcBase(1).fun), ord);
            domn = FtRefDomn;
        else
            baseDiv = trcBase(fES, ord, PmSym, name);
            val = evalBase(baseDiv, "ft_" + name, ord, X, {MshEnt(FtRefDomn).parm, PmSym}, {FtPmArg, PmArg}, {}, {});
            sDiv = size(baseDiv(1).fun);
            domn = baseDiv(1).domn;
        end
    end
    function baseDiv = difBase(fES, ord, name)
        % difBase: derivative of local base functions on element, cached by `name` (FE space) and `ord`.
        key = char("dif_" + name + "_" + mat2str(size(ord)) + "_" + mat2str(ord(:)'));
        if isKey(trcCache, key)
            baseDiv = trcCache(key);
            return;
        end
        baseDiv = fES.LcBase.dif(ord);
        trcCache(key) = baseDiv;
    end
    function baseDiv = trcBase(fES, ord, PmSym, name)
        % trcBase: derivative of local base functions as functions on reference facet.
        % Parameter of element base functions is replaced by `PmSym`; trace base functions keep facet parameter.
        % Results are cached by `name` (FE space) and `ord`: e.g. the forms from `DLF.interface` share them.
        key = char("trc_" + name + "_" + mat2str(size(ord)) + "_" + mat2str(ord(:)'));
        if isKey(trcCache, key)
            baseDiv = trcCache(key);
            return;
        end
        switch fES.elem
            case ElDomn
                baseDiv = fES.LcBase.dif(ord).subParm(PmSym).tfm(FtRefDomn);
            case FtRefDomn
                assert(all(ord == 0, "all"));
                baseDiv = fES.LcBase.dif(ord);
        end
        trcCache(key) = baseDiv;
    end
    function preDiv = trcPre(preSol, ord, PmSym, CfSym)
        % trcPre: derivative of previous solution as function on reference facet.
        % Coefficient is replaced by `CfSym`; parameter of element function is replaced by `PmSym`.
        switch preSol.elem
            case ElDomn
                preDiv = preSol.dif(ord).subParm(PmSym).subCoef(CfSym).tfm(FtRefDomn);
            case FtRefDomn
                assert(all(ord == 0, "all"));
                preDiv = preSol.dif(ord).subCoef(CfSym);
        end
    end
    function iEnt = cnEnt(fES, CnElem, FtIdx, iSide)
        % cnEnt: indices of entities of FE space (or FE function) connected to facets (by row).
        % Element space: element connected to each facet on the side given by sign of `iSide` (0 if none).
        % Trace space: the facets themselves.
        switch fES.elem
            case ElDomn
                iEnt = abs(sum(CnElem .* (sign(CnElem) == sign(iSide)), 1));
            case FtRefDomn
                iEnt = FtIdx;
        end
    end
    function asmStiff(I, J, V)
        % asmStiff: add values V to entries (I, J) of the stiffness matrix (I, J, V of the same size).
        n = numel(V);
        switch options.matType
            case "sparse"
                if nI + n > nImax
                    nImax = max(floor(nImax * 1.5), nI + n);
                    Is = [Is; zeros(nImax - length(Is), 1)];
                    Js = [Js; zeros(nImax - length(Js), 1)];
                    Vs = [Vs; zeros(nImax - length(Vs), 1)];
                end
                Is(nI + (1:n)) = I(:); Js(nI + (1:n)) = J(:); Vs(nI + (1:n)) = V(:);
                nI = nI + n;
            case "full"
                Stiff = Stiff + full(sparse(I(:), J(:), V(:), size(Stiff, 1), size(Stiff, 2)));
        end
    end
    function asmLoad(I, V)
        % asmLoad: add values V to entries I of the load vector (I, V of the same size).
        Load = Load + accumarray(I(:), V(:), size(Load));
    end
end
%% Local functions.
function [CfPh, CfSym, CfVal] = phCoef(cfs, X, ArgSym, Args, nEnt)
    % phCoef: placeholders for coefficient functions `cfs` (Fcn array), their symbols, and their values at points X of
    % a batch of mesh entities (nComp x nPnt x nEnt, or nComp x 1 x nEnt if constant).
    nCf = numel(cfs);
    CfPh(1:nCf) = Fcn.cst(0);
    CfSym = cell(1, nCf);
    CfVal = cell(1, nCf);
    for iCf = 1:nCf
        cf = cfs(iCf);
        CfSym{iCf} = sym(sprintf('cf%dVal', iCf), size(cf.fun));
        CfPh(iCf) = Fcn(cf.domn, CfSym{iCf});
        fun = cf.fun(:);
        if isempty(symvar(fun))
            CfVal{iCf} = repmat(double(fun), 1, 1, nEnt);
        else
            assert(~isempty(cf.var) && isempty(cf.coef));
            comps = repmat(Fcn.cst(0), 1, numel(fun));
            for iComp = 1:numel(fun)
                comps(iComp) = Fcn(cf.domn, fun(iComp));
            end
            funH = comps.getFun("parm", ArgSym, "vec", true);
            CfVal{iCf} = funH(X, Args{:});
        end
    end
end
function chunks = chunkIdx(n, sz)
    % chunkIdx: split 1:n into chunks (cell array, by column) such that a chunk holds about 2e7 values of size `sz` each.
    m = max(1, min(n, floor(2e7 / max(1, sz))));
    chunks = arrayfun(@(i0) i0:min(i0 + m - 1, n), 1:m:n, "UniformOutput", false);
end
function [XRep, ArgsRep, JacRep] = repEnt(X, Args, Jac, m)
    % repEnt: replicate points, batched arguments and Jacobians of mesh entities m times along entities.
    XRep = repmat(X, 1, 1, m);
    ArgsRep = cellfun(@(arg) repmat(arg, 1, 1, m), Args, "UniformOutput", false);
    JacRep = repmat(Jac, 1, 1, m);
end
