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
    % D3     | D3T      | 0      | D3        | D3
    % D3     | D3T      | 1      | D3 D3L    | D3 D3L
    % D3     | D3T      | 2      | D3 D3F    | D3 D3F
    % D3     | D3T      | 3      | D3 D3T    | D3 D3T
    % D3R2   | D2T      | 2      | D3R2      | D3R2
    % D3R2   | D3T      | 2      | D3R2      | D3 D3F

    properties
        coord (:, :) sym; % Relative coordinates of sample node in mesh entity (by column).
        % EntDim = 0: no relative coordinate is required.
        % EntDim = 1: use 1D coordinate `a` where the node is computed as (1 - a) * vertex_1 + a * vertex_2.
        % EntDim = 2: use 2D coordinates (a; b) where the node is computed as (1 - a - b) * vertex_1 + a * vertex_2 + b * vertex_3.
        % EntDim = 3: use 3D coordinates (a; b; c) where the node is computed as
        % (1 - a - b - c) * vertex_1 + a * vertex_2 + b * vertex_3 + c * vertex_4.
        % Remark: when `domn` is "D2" or "D3" and `EntDim` is 1, nodes must distribute symmetrically, e.g. coord = [1/3, 1/2 ,2/3].
        ord (:, :); % Order of derivative.
        % ord(i,j) = k: k-th derivative of j-th function component w.r.t. i-th variable.
        % Remark: when `domn` is "D2R1" or "D3R2", only zero order of derivative is supported.
        coef (1, :) Fcn; % Coefficient function.
        form; % Form of DoF: function handle of `coef` and `fcn`.
    end
    properties (Dependent)
        nEnt; % Number of mesh entities.
        nNode; % Number of sample nodes.
        nSamp; % Number of samplers.
        nDoF; % Number of degrees of freedom.
        sDoF; % Size of DoF group: [nEnt, nNode].
        loc (2, :); % Geometry location of sample node in mesh entity (by column).
        % Only supported for `domn` = "D2" and `EntDim` = 2.
        % - loc(1,j) = 0, loc(2,j) = k: j-th node is at k-th vertex of mesh entity.
        % - loc(1,j) = 1, loc(2,j) = k: j-th node is at k-th edge of mesh entity.
        LFun; % Length of function value.
    end
    methods
        %% Constructor.
        function NdDoF = NdDoF(domn, msh, EntDim, coord, ord, options)
            arguments
                domn {mustBeMember(domn, ["D2", "D2R1", "D3", "D3R2"])};
                msh Msh;
                EntDim {mustBeMember(EntDim, [0, 1, 2, 3])};
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
                if ismember(domn, ["D2", "D3"]) && EntDim < msh.dim
                    % DoFs on vertices, edges (and faces in 3D) are shared by default.
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
            assert(isequal(NdDoF.domn, "D2") && NdDoF.EntDim == 2);
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
            % NdDoF.eval: evaluate nodal DoF on function.
            % val(i,j): function's nodal DoF value at j-th sample node of i-th mesh entity.
            arguments
                ndDoF NdDoF;
                fcn Fcn;
                options.valType {mustBeMember(options.valType, ["sym", "num"])} = "sym"; % Value type.
                options.rawEval logical = false; % Raw evaluation.
                % Sample `fcn` directly, discarding the DoF's `ord`, `coef`, and `form`.
                % Only supported for scalar-valued function.
            end
            if options.rawEval
                assert(isscalar(fcn.fun));
                switch ndDoF.domn
                    case "D2"
                        ndDoF.ord = zeros(2, fcn.nFun);
                    case "D2R1"
                        ndDoF.ord = zeros(1, fcn.nFun);
                    case "D3"
                        ndDoF.ord = zeros(3, fcn.nFun);
                    case "D3R2"
                        ndDoF.ord = zeros(2, fcn.nFun);
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
                case {"D2", "D3"}
                    switch ndDoF.EntDim
                        case 0
                            fun = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord)).getFun;
                            node = ndDoF.msh.node.coord(:, ndDoF.EntIdx);
                            val(:) = fun(node);
                        otherwise
                            fcn = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord));
                            EntNode = ndDoF.msh.ent(ndDoF.EntDim).node;
                            if isequal(options.valType, "num")
                                % Evaluate on all mesh entities at once.
                                nEnt = ndDoF.nEnt;
                                if nEnt == 0
                                    return;
                                end
                                EntParm = reshape(ndDoF.msh.node.coord(:, EntNode(:, ndDoF.EntIdx)), ndDoF.msh.dim, size(EntNode, 1), nEnt);
                                node = pagemtimes(EntParm, [1 - sum(DoFCrd, 1); DoFCrd]);
                                fun = fcn.getFun("vec", true);
                                if ismember(fcn.domn, ["D2", "D3"])
                                    F = fun(node);
                                else
                                    F = fun(node, reshape(EntParm, [], 1, nEnt));
                                end
                                val = reshape(F + zeros(1, ndDoF.nNode, nEnt), ndDoF.nNode, nEnt).';
                                return;
                            end
                            fun = fcn.getFun;
                            for iEnt = 1:ndDoF.nEnt
                                EntParm = ndDoF.msh.node.coord(:, EntNode(:, ndDoF.EntIdx(iEnt)));
                                node = barNode(EntParm, DoFCrd);
                                if ismember(fcn.domn, ["D2", "D3"])
                                    val(iEnt, :) = fun(node);
                                else
                                    % Function defined on mesh entity, e.g. "D2L", "D3F".
                                    val(iEnt, :) = fun(node, EntParm);
                                end
                            end
                    end
                case {"D2R1", "D3R2"}
                    % Trace function on reference facet, e.g. "D2LR", "D3FR".
                    RefDomn = trcRefDomn(ndDoF.domn);
                    EntNode = ndDoF.msh.ent(ndDoF.EntDim).node;
                    if isequal(options.valType, "num")
                        % Evaluate on all mesh entities at once.
                        nEnt = ndDoF.nEnt;
                        if nEnt == 0
                            return;
                        end
                        EntParm = reshape(double(ndDoF.msh.node.coord(:, EntNode(:, ndDoF.EntIdx))), ndDoF.msh.dim, size(EntNode, 1), nEnt);
                        if ndDoF.msh.dim == MshEnt(RefDomn).dim
                            fun = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord)).getFun("vec", true);
                            F = fun(pagemtimes(EntParm, [1 - sum(DoFCrd, 1); DoFCrd]));
                        else
                            fun = ndDoF.form(ndDoF.coef, fcn.tfm(RefDomn).dif(ndDoF.ord)).getFun("vec", true);
                            node = barNode(double(MshEnt(RefDomn).node.coord), DoFCrd);
                            F = fun(repmat(node, 1, 1, nEnt), reshape(EntParm, [], 1, nEnt));
                        end
                        val = reshape(F + zeros(1, ndDoF.nNode, nEnt), ndDoF.nNode, nEnt).';
                        return;
                    end
                    if ndDoF.msh.dim == MshEnt(RefDomn).dim
                        % Mesh of reference facet itself.
                        fun = ndDoF.form(ndDoF.coef, fcn.dif(ndDoF.ord)).getFun;
                        for iEnt = 1:ndDoF.nEnt
                            node = barNode(ndDoF.msh.node.coord(:, EntNode(:, ndDoF.EntIdx(iEnt))), DoFCrd);
                            val(iEnt, :) = fun(node);
                        end
                    else
                        % Facets of mesh: sample function transformed to reference facet, parameterized by facet vertices.
                        fun = ndDoF.form(ndDoF.coef, fcn.tfm(RefDomn).dif(ndDoF.ord)).getFun;
                        node = barNode(MshEnt(RefDomn).node.coord, DoFCrd);
                        for iEnt = 1:ndDoF.nEnt
                            FtParm = ndDoF.msh.node.coord(:, EntNode(:, ndDoF.EntIdx(iEnt)));
                            val(iEnt, :) = fun(node, FtParm);
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
        case "D3"
            assert(ismember(mshType, "D3T"));
        case "D3R2"
            assert(ismember(mshType, ["D2T", "D3T"]));
    end
    switch domn
        case "D2"
            assert(ismember(EntDim, [0, 1, 2]));
        case "D2R1"
            assert(ismember(EntDim, 1));
        case "D3"
            assert(ismember(EntDim, [0, 1, 2, 3]));
        case "D3R2"
            assert(ismember(EntDim, 2));
    end
    assert(size(coord, 1) == EntDim);
    if EntDim >= 1
        assert(all(coord(:) >= 0) && all(sum(coord, 1) <= 1));
    end
    if ismember(domn, ["D2", "D3"]) && EntDim == 1
        assert(all(abs(coord - (1 - coord(end:-1:1))) < 10 * eps));
    end
    switch domn
        case "D2"
            assert(size(ord, 1) == 2);
        case "D2R1"
            assert(size(ord, 1) == 1);
            assert(all(ord == 0));
        case "D3"
            assert(size(ord, 1) == 3);
        case "D3R2"
            assert(size(ord, 1) == 2);
            assert(all(ord == 0, "all"));
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
            case "D3"
                switch EntDim
                    case {0, 1, 2}
                        assert(ismember(share, [true, false]));
                    case 3
                        assert(ismember(share, false));
                end
            case "D3R2"
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
        case "D3"
            switch EntDim
                case {0, 3}
                    assert(ismember(orien, false));
                case {1, 2}
                    assert(ismember(orien, [true, false]));
            end
        case "D3R2"
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
        case "D3"
            switch EntDim
                case 0
                    assert(ismember(coef.getDomn, ["VOID", "D3"]));
                case 1
                    assert(ismember(coef.getDomn, ["VOID", "D3", "D3L"]));
                case 2
                    assert(ismember(coef.getDomn, ["VOID", "D3", "D3F"]));
                case 3
                    assert(ismember(coef.getDomn, ["VOID", "D3", "D3T"]));
            end
        case "D3R2"
            assert(ismember(coef.getDomn, ["VOID", "D3R2"]));
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
        case "D3"
            switch EntDim
                case 0
                    assert(ismember(fcn.domn, "D3"));
                case 1
                    assert(ismember(fcn.domn, ["D3", "D3L"]));
                case 2
                    assert(ismember(fcn.domn, ["D3", "D3F"]));
                case 3
                    assert(ismember(fcn.domn, ["D3", "D3T"]));
            end
        case "D3R2"
            switch mshType
                case "D2T"
                    assert(ismember(fcn.domn, "D3R2"));
                case "D3T"
                    assert(ismember(fcn.domn, ["D3", "D3F"]));
            end
    end
    assert(fcn.nFun == size(ord, 2));
end
function RefDomn = trcRefDomn(domn)
    % trcRefDomn: reference facet of trace function domain.
    switch domn
        case "D2R1"
            RefDomn = "D2LR";
        case "D3R2"
            RefDomn = "D3FR";
    end
end
function node = barNode(EntParm, DoFCrd)
    % barNode: compute sample nodes from vertices of mesh entity and relative coordinates (by column).
    % node(:,j) = (1 - sum(DoFCrd(:,j))) * EntParm(:,1) + EntParm(:,2:end) * DoFCrd(:,j).
    node = EntParm * [1 - sum(DoFCrd, 1); DoFCrd];
end
