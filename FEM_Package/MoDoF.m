classdef MoDoF < DoF
    % MoDoF: moment degree of freedom.
    % Integration of `form(coef, fcn.dif(ord))` multiplied by specific test function.

    % Valid properties:
    % domn   | msh.type | EntDim| tst.domn  | coef.domn | fcn.domn
    % -------|----------|-------|-----------|-----------|---------
    % D2     | D2T      | 1     | D2R1 D2LR | D2 D2L    | D2 D2L
    % D2     | D2T      | 2     | D2R  D2TR | D2 D2T    | D2 D2T
    % D2R1   | D1       | 1     | D2R1      | D2R1      | D2R1
    % D2R1   | D2T      | 1     | D2R1      | D2R1      | D2 D2L

    properties
        tst (1, :) Fcn; % Test functions.
        % Remark: when `domn` is "D2" and `EntDim` is 1, test functions must be symmetric (either odd or even) on the reference line [0, 1].
        ord (:, :); % Order of derivative.
        % ord(i,j) = k: k-th derivative of j-th function component w.r.t i-th variable.
        % Remark: when `domn` is "D2R1", only zero order of derivative is supported.
        coef (1, :) Fcn; % Coefficient function.
        form % Form of DoF: function handle of `coef`, `fcn`, and `tst`.
        GInt GInt; % Gaussian integration.
    end
    properties (Dependent)
        nEnt; % Number of mesh entities.
        nTst; % Number of test functions.
        nSamp; % Number of samplers.
        nDoF; % Number of degree of freedoms.
        sDoF; % Size of DoF group: [nEnt, nTst].
        LFun; % Length of function value.
    end
    methods
        %% Constructor.
        function MoDoF = MoDoF(domn, msh, EntDim, tst, ord, options)
            arguments
                domn {mustBeMember(domn, ["D2", "D2R1"])};
                msh Msh;
                EntDim {mustBeMember(EntDim, [1, 2])};
                tst (1, :) Fcn;
                ord (:, :);
                options.EntIdx (1, :) = [];
                options.share logical = [];
                options.orien logical = false;
                options.coef (1, :) Fcn = Fcn.cst(1);
                options.form = @(coef, fcn, tst) sum(coef .* fcn .* tst);
                options.GInt GInt = GInt.empty;
            end
            isTstOdd = checkProp(domn, msh.type, EntDim, tst, ord, options.share, options.orien, options.coef, options.GInt);
            MoDoF.domn = domn;
            MoDoF.msh = msh;
            MoDoF.EntDim = EntDim;
            MoDoF.tst = tst;
            MoDoF.ord = ord;
            if ~isempty(options.EntIdx)
                MoDoF.EntIdx = options.EntIdx;
            else
                MoDoF.EntIdx = 1:msh.nEnt(EntDim);
            end
            if ~isempty(options.share)
                MoDoF.share = options.share;
            else
                if isequal(domn, "D2") && ismember(EntDim, 1)
                    MoDoF.share = true;
                else
                    MoDoF.share = false;
                end
            end
            MoDoF.orien = options.orien;
            if isTstOdd
                % Note: if test function is odd, `orien` is reversed.
                MoDoF.orien = ~MoDoF.orien;
            end
            MoDoF.coef = options.coef;
            MoDoF.form = options.form;
            if ~isempty(options.GInt)
                MoDoF.GInt = options.GInt;
            else
                warning("No GInt specified.");
                switch EntDim
                    case 2
                        MoDoF.GInt = GInt("D2T", 1);
                    case 1
                        MoDoF.GInt = GInt("D2L", 1);
                end
            end
        end
        % Get functions.
        function nEnt = get.nEnt(MoDoF)
            nEnt = length(MoDoF.EntIdx);
        end
        function nTst = get.nTst(MoDoF)
            nTst = length(MoDoF.tst);
        end
        function nSamp = get.nSamp(MoDoF)
            nSamp = MoDoF.nTst;
        end
        function nDoF = get.nDoF(MoDoF)
            nDoF = MoDoF.nEnt * MoDoF.nTst;
        end
        function sDoF = get.sDoF(MoDoF)
            sDoF = [MoDoF.nEnt, MoDoF.nTst];
        end
        function LFun = get.LFun(MoDoF)
            LFun = size(MoDoF.ord, 2);
        end
        %% Public functions.
        function val = eval(moDoF, fcn, options)
            % MoDoF.eval: evaluate moment DoF on function.
            % val(i,j): function's moment DoF value w.r.t. j-th test function on i-th mesh entity.
            arguments
                moDoF MoDoF;
                fcn Fcn;
                options.valType {mustBeMember(options.valType, ["sym", "num"])} = "sym"; % Value type.
                options.rawEval logical = false; % Raw evaluation.
                % Sample `fcn` directly, discarding the DoF's `ord`, `coef`, and `form`.
                % Only supported for scalar-valued function.
            end
            if options.rawEval
                % assert(isscalar(fcn.fun));
                switch moDoF.domn
                    case "D2"
                        moDoF.ord = zeros(2, fcn.nFun);
                    case "D2R1"
                        moDoF.ord = zeros(1, fcn.nFun);
                end
                moDoF.coef = Fcn.cst(1);
                moDoF.form = @(coef, fcn, tst) sum(coef .* fcn .* tst);
            end
            checkFcn(moDoF.domn, moDoF.msh.type, moDoF.EntDim, moDoF.ord, fcn);
            switch options.valType
                case "sym"
                    val = sym(zeros(moDoF.nEnt, moDoF.nTst));
                    switch moDoF.domn
                        case "D2"
                            switch moDoF.EntDim
                                case 1
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef.tfm("D2LR"), fcn.dif(moDoF.ord).tfm("D2LR"), moDoF.tst(iTst)) .* Tfm("D2L").JNorm;
                                        intVal = intFcn.int("D2LR").getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            iEdge = moDoF.EntIdx(iEnt);
                                            EgParm = moDoF.msh.node.coord(:, moDoF.msh.edge.node(:, iEdge));
                                            val(iEnt, iTst) = intVal([0; 0], EgParm);
                                        end
                                    end
                                case 2
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef.tfm("D2TR"), fcn.dif(moDoF.ord).tfm("D2TR"), moDoF.tst(iTst)) .* Tfm("D2T").JDet;
                                        intVal = intFcn.int("D2TR").getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            iElem = moDoF.EntIdx(iEnt);
                                            ElParm = moDoF.msh.node.coord(:, moDoF.msh.elem.node(:, iElem));
                                            val(iEnt, iTst) = intVal([0; 0], ElParm);
                                        end
                                    end
                            end
                        case "D2R1"
                            switch moDoF.msh.type
                                case "D1"
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef, fcn.dif(moDoF.ord), moDoF.tst(iTst));
                                        intVal = intFcn.int("D2LR").getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            val(iEnt, iTst) = intVal(0);
                                        end
                                    end
                                case "D2T"
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef, fcn.tfm("D2LR").dif(moDoF.ord), moDoF.tst(iTst));
                                        intVal = intFcn.int("D2LR").getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            iEdge = moDoF.EntIdx(iEnt);
                                            EgParm = moDoF.msh.node.coord(:, moDoF.msh.edge.node(:, iEdge));
                                            val(iEnt, iTst) = intVal(0, EgParm);
                                        end
                                    end
                            end
                    end
                    val = simplify(sym(val));
                case "num"
                    val = zeros(moDoF.nEnt, moDoF.nTst);
                    switch moDoF.domn
                        case "D2"
                            switch moDoF.EntDim
                                case 1
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef.tfm("D2LR"), fcn.dif(moDoF.ord).tfm("D2LR"), moDoF.tst(iTst)) .* Tfm("D2L").JNorm;
                                        intFun = intFcn.getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            iEdge = moDoF.EntIdx(iEnt);
                                            EgParm = moDoF.msh.node.coord(:, moDoF.msh.edge.node(:, iEdge));
                                            val(iEnt, iTst) = moDoF.GInt.eval(@(x) intFun(x, EgParm));
                                        end
                                    end
                                case 2
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef.tfm("D2TR"), fcn.dif(moDoF.ord).tfm("D2TR"), moDoF.tst(iTst)) .* Tfm("D2T").JDet;
                                        intFun = intFcn.getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            iElem = moDoF.EntIdx(iEnt);
                                            ElParm = moDoF.msh.node.coord(:, moDoF.msh.elem.node(:, iElem));
                                            val(iEnt, iTst) = moDoF.GInt.eval(@(x) intFun(x, ElParm));
                                        end
                                    end
                            end
                        case "D2R1"
                            switch moDoF.msh.type
                                case "D1"
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef, fcn.dif(moDoF.ord), moDoF.tst(iTst));
                                        intFun = intFcn.getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            val(iEnt, iTst) = moDoF.GInt.eval(@(x) intFun(x));
                                        end
                                    end
                                case "D2T"
                                    for iTst = 1:moDoF.nTst
                                        intFcn = moDoF.form(moDoF.coef, fcn.tfm("D2LR").dif(moDoF.ord), moDoF.tst(iTst));
                                        intFun = intFcn.getFun;
                                        for iEnt = 1:moDoF.nEnt
                                            iEdge = moDoF.EntIdx(iEnt);
                                            EgParm = moDoF.msh.node.coord(:, moDoF.msh.edge.node(:, iEdge));
                                            val(iEnt, iTst) = moDoF.GInt.eval(@(x) intFun(x, EgParm));
                                        end
                                    end
                            end
                    end
            end
        end
    end
end
%% Local functions.
function isTstOdd = checkProp(domn, mshType, EntDim, tst, ord, share, orien, coef, GInt)
    % checkProp: check validity of properties.
    isTstOdd = false;
    switch domn
        case "D2"
            assert(ismember(mshType, "D2T"));
        case "D2R1"
            assert(ismember(mshType, ["D1", "D2T"]));
    end
    switch domn
        case "D2"
            assert(ismember(EntDim, [1, 2]));
        case "D2R1"
            assert(ismember(EntDim, 1));
    end
    switch domn
        case "D2"
            switch EntDim
                case 1
                    assert(all(ismember([tst.domn], ["VOID", "D2R1", "D2LR"])));
                    isSym = zeros(1, length(tst));
                    for iTst = 1:length(tst)
                        tstFun = tst(iTst).fun;
                        tstVar = tst(iTst).var;
                        tstFunInv = subs(tstFun, tstVar, 1 - tstVar);
                        switch tstFun
                            case tstFunInv
                                isSym(iTst) = 1;
                            case -tstFunInv
                                isSym(iTst) = -1;
                        end
                        assert(all(isSym == 1) || all(isSym == -1));
                        if all(isSym == -1)
                            isTstOdd = true;
                        end
                    end
                case 2
                    assert(all(ismember([tst.domn], ["VOID", "D2R", "D2TR"])));
            end
        case "D2R1"
            assert(all(ismember([tst.domn], ["VOID", "D2R1"])));
    end
    switch domn
        case "D2"
            assert(size(ord, 1) == 2);
        case "D2R1"
            assert(size(ord, 1) == 1);
            assert(all(ord == 0));
    end
    if ~isempty(share)
        switch domn
            case "D2"
                switch EntDim
                    case 1
                        assert(ismember(share, [true, false]));
                    case 2
                        assert(ismember(share, false));
                end
            case "D2R1"
                assert(ismember(share, false));
        end
    end
    switch domn
        case "D2"
            switch EntDim
                case 1
                    assert(ismember(orien, [true, false]));
                case 2
                    assert(ismember(orien, false));
            end
        case "D2R1"
            assert(ismember(orien, false));
    end
    switch domn
        case "D2"
            switch EntDim
                case 1
                    assert(ismember(coef.getDomn, ["VOID", "D2", "D2L"]));
                case 2
                    assert(ismember(coef.getDomn, ["VOID", "D2", "D2T"]));
            end
        case "D2R1"
            assert(ismember(coef.getDomn, ["VOID", "D2R1"]));
    end
    if ~isempty(GInt)
        switch EntDim
            case 2
                assert(isequal(GInt.domn, "D2T"));
            case 1
                assert(isequal(GInt.domn, "D2L"));
        end
    end
end
function checkFcn(domn, mshType, EntDim, ord, fcn)
    % checkFcn: check validity of function.
    switch domn
        case "D2"
            switch EntDim
                case 1
                    assert(ismember(fcn.domn, ["D2", "D2L"]));
                case 2
                    assert(ismember(fcn.domn, ["D2", "D2T"]));
            end
        case "D2R1"
            switch mshType
                case "D1"
                    assert(ismember(fcn.domn, "D2R1"));
                case "D2T"
                    assert(ismember(fcn.domn, ["D2", "D2L"]));
            end
    end
    assert(fcn.nFun == size(ord, 2));
end
