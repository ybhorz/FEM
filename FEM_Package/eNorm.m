function eNorm = eNorm(msh, reals, apprs, norms)
    % eNorm: compute error norm between real and approximate solutions.
    % Several norms (with the same power p) are combined as (sum_i |e|_i^p)^(1/p), e.g. the broken H1 norm
    % |v|_{1,h}^2 = sum_K |grad v|_K^2 + sum_e h_e^{-1} |[v]|_e^2 is given by an element norm and a jump norm.
    arguments
        msh Msh; % Mesh.
        reals (1, :) Fcn; % Real solutions.
        apprs (1, :) FEF; % Approximate solutions.
        norms (1, :) Norm; % Norms.
    end
    assert(ismember(msh.type, ["D2T", "D3T"]));
    ElDomn = msh.ElDomn; % Domain of element.
    FtDomn = msh.FtDomn; % Domain of facet.
    FtRefDomn = msh.FtRefDomn; % Reference domain of facet.
    facet = msh.facet; % Facet.
    assert(isequal(msh, apprs.getMsh));
    assert(isequal(msh, norms.getMsh));
    eNorms = zeros(1, length(norms));
    for iNorm = 1:length(norms)
        norm = norms(iNorm);
        real = reals(norm.iFcn);
        appr = apprs(norm.iFcn);
        % Approximate solution (and its derivatives) is evaluated numerically at Gauss points of all mesh entities;
        % the norm is turned into a function handle with placeholder `apprVal` for these values (cf. `assemble`).
        switch norm.EntDim
            case msh.dim
                assert(ismember(real.domn, ["VOID", msh.domn, ElDomn]));
                assert(ismember(appr.elem, ElDomn));
                ElIdx = norm.EntIdx;
                nEnt = length(ElIdx);
                ElParm = appr.ElParm(:, :, ElIdx);
                ElPmArg = reshape(ElParm, [], 1, nEnt);
                ElCfArg = reshape(appr.ElCoef(:, ElIdx), [], 1, nEnt);
                [X, Jac] = norm.GInt.pntVec(ElParm);
                apprDiv = appr.dif(norm.ord);
                ApprVal = evalVal(apprDiv, X, {}, {ElPmArg}, {}, {ElCfArg});
                % Real solution broadcast to the size of approximate solution.
                realDiv = dif(real + Fcn(msh.domn, sym(zeros(size(appr.fun)))), norm.ord);
                ApprSym = sym('apprVal', size(apprDiv.fun));
                intFcn = norm.form(norm.coef, realDiv - Fcn(ElDomn, ApprSym), norm.pow);
                intFun = intFcn.getFun("parm", {MshEnt(ElDomn).parm}, "coef", {ApprSym}, "vec", true);
                intVal = norm.GInt.sumVec(intFun(X, ElPmArg, ApprVal), Jac);
                eNorms(iNorm) = sum(intVal) ^ (1 / norm.pow);
            case msh.dim - 1
                assert(ismember(appr.elem, [FtRefDomn, ElDomn]));
                switch appr.elem
                    case FtRefDomn
                        assert(ismember(real.domn, ["VOID", msh.domn, FtDomn]));
                        assert(all(norm.ord == 0));
                        assert(isequal(norm.fcnOpr, "none"));
                        FtIdx = norm.EntIdx;
                        nEnt = length(FtIdx);
                        FtPmArg = reshape(appr.ElParm(:, :, FtIdx), [], 1, nEnt);
                        FtCfArg = reshape(appr.ElCoef(:, FtIdx), [], 1, nEnt);
                        [X, Jac] = norm.GInt.pntVec([], nEnt);
                        apprDiv = appr.dif(norm.ord);
                        ApprVal = evalVal(apprDiv, X, {}, {FtPmArg}, {}, {FtCfArg});
                        realDiv = real.tfm(FtRefDomn) + Fcn(FtRefDomn, sym(zeros(size(appr.fun))));
                        ApprSym = sym('apprVal', size(apprDiv.fun));
                        intFcn = norm.form(norm.coef.tfm(FtRefDomn), realDiv - Fcn(FtRefDomn, ApprSym), norm.pow) .* Tfm(FtDomn).JNorm;
                        intFun = intFcn.getFun("parm", {MshEnt(FtRefDomn).parm}, "coef", {ApprSym}, "vec", true);
                        intVal = norm.GInt.sumVec(intFun(X, FtPmArg, ApprVal), Jac);
                        eNorms(iNorm) = sum(intVal) ^ (1 / norm.pow);
                    case ElDomn
                        assert(ismember(real.domn, ["VOID", msh.domn]));
                        assert(all(norm.ord == 0, "all"));
                        assert(isequal(norm.fcnOpr, "jump"));
                        % Interior facets only.
                        FtIdx = norm.EntIdx;
                        CnElem = abs(facet.elem(:, FtIdx));
                        isVal = all(CnElem ~= 0, 1);
                        FtIdx = FtIdx(isVal); iElem1 = CnElem(1, isVal); iElem2 = CnElem(2, isVal);
                        nEnt = length(FtIdx);
                        if nEnt == 0
                            continue
                        end
                        FtParm = reshape(msh.node.coord(:, facet.node(:, FtIdx)), msh.dim, facet.nNode, nEnt);
                        [X, Jac] = norm.GInt.pntVec(FtParm);
                        % Values of approximate solution from both connected elements at Gauss points of facets.
                        apprDiv = appr.dif(norm.ord);
                        El1Val = evalVal(apprDiv, X, {}, {reshape(appr.ElParm(:, :, iElem1), [], 1, nEnt)}, {}, ...
                            {reshape(appr.ElCoef(:, iElem1), [], 1, nEnt)});
                        El2Val = evalVal(apprDiv, X, {}, {reshape(appr.ElParm(:, :, iElem2), [], 1, nEnt)}, {}, ...
                            {reshape(appr.ElCoef(:, iElem2), [], 1, nEnt)});
                        JumpSym = sym('jumpVal', size(apprDiv.fun));
                        intFcn = norm.form(norm.coef, Fcn(FtDomn, JumpSym), norm.pow);
                        intFun = intFcn.getFun("parm", {MshEnt(FtDomn).parm}, "coef", {JumpSym}, "vec", true);
                        intVal = norm.GInt.sumVec(intFun(X, reshape(FtParm, [], 1, nEnt), El1Val - El2Val), Jac);
                        eNorms(iNorm) = sum(intVal) ^ (1 / norm.pow);
                end
        end
    end
    pow = [norms.pow];
    assert(all(pow == pow(1)));
    eNorm = sum(eNorms .^ pow(1)) ^ (1 / pow(1));
end
%% Local functions.
function val = evalVal(fcn, X, ParmSym, ParmArg, CoefSym, CoefArg)
    % evalVal: values of (vector- or matrix-valued) function at points X of a batch of mesh entities,
    % val(k, :, :) is k-th component, nComp x nPnt x nEnt (see `Fcn.getFun` with "vec").
    nComp = numel(fcn.fun);
    comps(1:nComp) = Fcn.cst(0);
    for iComp = 1:nComp
        comps(iComp) = Fcn(fcn.domn, fcn.fun(iComp), fcn.coef);
    end
    funH = comps.getFun("parm", ParmSym, "coef", CoefSym, "vec", true);
    val = funH(X, ParmArg{:}, CoefArg{:});
end
