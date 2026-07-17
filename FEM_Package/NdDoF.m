classdef NdDoF < DoF
    % NdDoF: nodal degree of freedom.
    % Function value of `form(coef, fcn.dif(ord))` at specific node.

    % Valid properties:
    % domn   | msh.type | EntDim | coef.domn | fcn.domn
    %-----------------------------------------------------
    % D2     | D2T      | 0      | D2        | D2
    % D2     | D2T      | 1      | D2 D2L    | D2 D2L
    % D2     | D2T      | 2      | D2 D2T    | D2 D2T
    % D2R1   | D1       | 1      | D2R1      | D2R1
    % D2R1   | D2T      | 1      | D2R1      | D2 D2L

    properties
        coord (:, :) sym; % Relative coordinates of sample node in mesh entity (by column).
        % EntDim = 0: no relative coordinate is required.
        % EntDim = 1: use 1D coordinate "a" where the node is computed as (1 - a) * vertex_1 + a * vertex_2.
        % EntDim = 2: use 2D coordinates (a, b) where the node is computed as (1 - a - b) * vertex_1 + a * vertex_2 + b * vertex_3.
        % When `domn` is "D2", `msh.type` is "D2T" and `EntDim` is 1, node are required to distribute symmetrically, e.g. coord = [1/3, 1/2 ,2/3].
        ord (:, :); % Order of derivative.
        % ord(i,j) = k: k-th derivative of j-th function component with respect to i-th variable.
        % When `domn` is "D2R1", order of derivative is required to be all zero.
        coef (1, :) Fcn; % Coefficient function.
        form; % Form of DoF: function handle of `coef` and `fcn`.
    end
    properties (Dependent)
        nEnt; % Number of mesh entities.
        nNode; % Number of sample node.
        nSamp; % Number of samplers.
        nDoF; % Number of degree of freedoms.
        sDoF; % Size of DoF group: [nEnt, nNode].
        loc (2, :); % Geometry location of sample node in mesh entity (by column).
        % Only supported for `EntDim` = 2.
        % loc(1,j) = 0, loc(2,j) = k: j-th node is at k-th vertex of mesh entity.
        % loc(1,j) = 1, loc(2,j) = k: j-th node is at k-th edge of mesh entity.
        LFun; % Length of function value.
    end
    methods
        %% Constructor.
        function NdDoF = NdDoF(domn, msh, EntDim, coord, ord, options)
            arguments
                domn {mustBeMember(domn, ["D2", "D2R1"])};
                msh Msh;
                EntDim {mustBeMember(EntDim, [0, 1, 2])};
                coord (:, :) {mustBeA(coord, ["sym", "string", "double"])};
                ord (:, :);
                options.EntIdx (1, :) = [];
                options.share logical = [];
                options.orien logical = false;
                options.coef (1, :) Fcn = Fcn.cst(1);
                options.form = @(coef, fcn) sum(coef .* fcn);
            end
            checkProp(domn, msh.type, EntDim, coord, ord, options.share, options.orien, options.coef);
            NdDoF.domn = domn;
            NdDoF.msh = msh;
            NdDoF.EntDim = EntDim;
            if isstring(coord)
                coord = str2sym(coord);
            end
            if isnumeric(coord)
                coord = sym(coord);
            end
            NdDoF.coord = coord;
            NdDoF.ord = ord;
            if ~isempty(options.EntIdx)
                NdDoF.EntIdx = options.EntIdx;
            else
                NdDoF.EntIdx = 1:msh.nEnt(EntDim);
            end
            if ~isempty(options.share)
                NdDoF.share = options.share;
            else
                if isequal(domn, "D2") && ismember(EntDim, [0, 1])
                    NdDoF.share = true;
                else
                    NdDoF.share = false;
                end
            end
            NdDoF.orien = options.orien;
            NdDoF.coef = options.coef;
            NdDoF.form = options.form;
        end
        %% Get functions.
        function nEnt = get.nEnt(NdDoF)
            nEnt = length(NdDoF.EntIdx);
        end
        function nNode = get.nNode(NdDoF)
            if NdDoF.EntDim == 0
                nNode = 1;
            else
                nNode = size(NdDoF.coord, 2);
            end
        end
        function nSamp = get.nSamp(NdDoF)
            nSamp = NdDoF.nNode;
        end
        function nDoF = get.nDoF(NdDoF)
            nDoF = NdDoF.nEnt * NdDoF.nNode;
        end
        function sDoF = get.sDoF(NdDoF)
            sDoF = [NdDoF.nEnt, NdDoF.nNode];
        end
        function loc = get.loc(NdDoF)
            assert(NdDoF.EntDim == 2);
            loc = zeros(2, NdDoF.nNode);
            loc(1, :) = 2;
            loc(2, :) = 1;
            for iNode = 1:NdDoF.nNode
                if isequal(NdDoF.coord(:, iNode), sym([0; 0]))
                    loc(1, iNode) = 0;
                    loc(2, iNode) = 1;
                elseif isequal(NdDoF.coord(:, iNode), sym([1; 0]))
                    loc(1, iNode) = 0;
                    loc(2, iNode) = 2;
                elseif isequal(NdDoF.coord(:, iNode), sym([0; 1]))
                    loc(1, iNode) = 0;
                    loc(2, iNode) = 3;
                elseif isequal(NdDoF.coord(2, iNode), sym(0))
                    loc(1, iNode) = 1;
                    loc(2, iNode) = 1;
                elseif isequal(sum(NdDoF.coord(:, iNode)), sym(1))
                    loc(1, iNode) = 1;
                    loc(2, iNode) = 2;
                elseif isequal(NdDoF.coord(1, iNode), sym(0))
                    loc(1, iNode) = 1;
                    loc(2, iNode) = 3;
                end
            end
        end
        function LFun = get.LFun(NdDoF)
            LFun = size(NdDoF.ord, 2);
        end
        %% Public functions.
        function val = eval(ndDoF, fcn, options)
            % NdDoF.eval: eval nodal DoF on function.
            % val(i,j): function's nodal DoF value at j-th sample node of i-th mesh entity.
            arguments
                ndDoF NdDoF;
                fcn Fcn;
                options.valType {mustBeMember(options.valType, ["sym", "num"])} = "sym"; % Value type.
                options.rawEval logical = false; % Whether evaluate raw function directly.
            end
            if options.rawEval
                switch ndDoF.domn
                    case "D2"
                        ndDoF.ord = zeros(2, fcn.nFun);
                    case "D2R1"
                        ndDoF.ord = zeros(1, fcn.nFun);
                end
                ndDoF.coef = Fcn.cst(1);
                ndDoF.form = @(coef, fcn) sum(coef .* fcn);
            end
            checkFcn(ndDoF.domn, ndDoF.msh.type, ndDoF.EntDim, ndDoF.ord, fcn);
            switch options.valType
                case "sym"
                    DoFCrd = ndDoF.coord;
                    val = sym(zeros(ndDoF.nEnt, ndDoF.nNode));
                case "num"
                    DoFCrd = double(ndDoF.coord);
                    val = zeros(ndDoF.nEnt, ndDoF.nNode);
            end
            switch ndDoF.domn
                case "D2"
                    switch ndDoF.EntDim
                        case 0
                            fun = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord)).getFun;
                            node = ndDoF.msh.node.coord(:, ndDoF.EntIdx);
                            val(:) = fun(node);
                        case 1
                            fcn = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord));
                            fun = fcn.getFun;
                            for iEnt = 1:ndDoF.nEnt
                                iEdge = ndDoF.EntIdx(iEnt);
                                EgNd1 = ndDoF.msh.node.coord(:, ndDoF.msh.edge.node(1, iEdge));
                                EgNd2 = ndDoF.msh.node.coord(:, ndDoF.msh.edge.node(2, iEdge));
                                a = DoFCrd;
                                node = (1 - a) .* EgNd1 + a .* EgNd2;
                                switch fcn.domn
                                    case "D2"
                                        val(iEnt, :) = fun(node);
                                    case "D2L"
                                        EgParm = [EgNd1, EgNd2];
                                        val(iEnt, :) = fun(node, EgParm);
                                end
                            end
                        case 2
                            fcn = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord));
                            fun = fcn.getFun;
                            for iEnt = 1:ndDoF.nEnt
                                iElem = ndDoF.EntIdx(iEnt);
                                ElNd1 = ndDoF.msh.node.coord(:, ndDoF.msh.elem.node(1, iElem));
                                ElNd2 = ndDoF.msh.node.coord(:, ndDoF.msh.elem.node(2, iElem));
                                ElNd3 = ndDoF.msh.node.coord(:, ndDoF.msh.elem.node(3, iElem));
                                a = DoFCrd(1, :);
                                b = DoFCrd(2, :);
                                node = (1 - a - b) .* ElNd1 + a .* ElNd2 + b .* ElNd3;
                                switch fcn.domn
                                    case "D2"
                                        val(iEnt, :) = fun(node);
                                    case "D2T"
                                        ElParm = [ElNd1, ElNd2, ElNd3];
                                        val(iEnt, :) = fun(node, ElParm);
                                end
                            end
                    end
                case "D2R1"
                    switch ndDoF.msh.type
                        case "D1"
                            fun = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord)).getFun;
                            for iEnt = 1:ndDoF.nEnt
                                iElem = ndDoF.EntIdx(iEnt);
                                ElNd1 = ndDoF.msh.node.coord(:, ndDoF.msh.elem.node(1, iElem));
                                ElNd2 = ndDoF.msh.node.coord(:, ndDoF.msh.elem.node(2, iElem));
                                a = DoFCrd;
                                node = (1 - a) .* ElNd1 + a .* ElNd2;
                                val(iEnt, :) = fun(node);
                            end
                        case "D2T"
                            fun = ndDoF.form(ndDoF.coef, fcn.tfm("D2LR").dif(ndDoF.ord)).getFun;
                            ElNd1 = MshEnt("D2LR").node.coord(:, 1);
                            ElNd2 = MshEnt("D2LR").node.coord(:, 2);
                            a = DoFCrd;
                            node = (1 - a) .* ElNd1 + a .* ElNd2;
                            for iEnt = 1:ndDoF.nEnt
                                iEdge = ndDoF.EntIdx(iEnt);
                                EgParm = ndDoF.msh.node.coord(:, ndDoF.msh.edge.node(:, iEdge));
                                val(iEnt, :) = fun(node, EgParm);
                            end
                    end
            end
            if isequal(options.valType, "sym")
                val = simplify(sym(val));
            end
        end
    end
end
%% Local functions.
function checkProp(domn, mshType, EntDim, coord, ord, share, orien, coef)
    % checkProp: check validity of properties.
    switch domn
        case "D2"
            assert(ismember(mshType, "D2T"));
        case "D2R1"
            assert(ismember(mshType, ["D1", "D2T"]));
    end
    switch domn
        case "D2"
            assert(ismember(EntDim, [0, 1, 2]));
        case "D2R1"
            assert(ismember(EntDim, 1));
    end
    assert(size(coord, 1) == EntDim);
    switch EntDim
        case 1
            assert(all(coord >= 0) && all(coord <= 1));
        case 2
            assert(all(coord(1, :) >= 0) && all(coord(2, :) >= 0) && all(sum(coord, 1) <= 1));
    end
    if isequal(domn, "D2") && ismember(mshType, "D2T") && EntDim == 1
        assert(all(abs(coord - (1 - coord(end:-1:1))) < 10 * eps));
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
                    case {0, 1}
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
                case {0, 2}
                    assert(ismember(orien, false));
                case 1
                    assert(ismember(orien, [true, false]));
            end
        case "D2R1"
            assert(ismember(orien, false));
    end
    switch domn
        case "D2"
            switch EntDim
                case 0
                    assert(ismember(coef.getDomn, ["VOID", "D2"]));
                case 1
                    assert(ismember(coef.getDomn, ["VOID", "D2", "D2L"]));
                case 2
                    assert(ismember(coef.getDomn, ["VOID", "D2", "D2T"]));
            end
        case "D2R1"
            assert(ismember(coef.getDomn, ["VOID", "D2R1"]));
    end
end
function checkFcn(domn, mshType, EntDim, ord, fcn)
    % checkFcn: check validity of function.
    switch domn
        case "D2"
            switch EntDim
                case 0
                    assert(ismember(fcn.domn, "D2"));
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
