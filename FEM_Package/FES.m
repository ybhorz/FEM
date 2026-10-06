classdef FES
    % FES: Finite Element space.
    properties
        msh Msh; % Mesh.
        elem {mustBeMember(elem, ["VOID", "D2T", "D2LR", "D3T", "D3FR"])} = "VOID"; % Element type.
        LcBase (1, :) Fcn; % Local base function.
        map {mustBeMember(map, ["none", "affine", "piolaDiv", "piolaCurl"])} = "none"; % Mapping of base functions (see `FE`).
        RefBase (1, :) Fcn; % Base function on reference element before mapping (empty if map is "none").
        RefKey (1, 1) string = ""; % Key of reference base functions (see `FE`).
        ElParm (:, :, :); % Parameter of base function on each element.
        % ElParm(:, :, i): parameter for i-th element.
        GlDoFs (1, :) DoF = NdDoF.empty; % Global degrees of freedom.
        Lc2Gl (:, :); % Mapping from local DoFs to global DoFs.
        % Lc2Gl(i, j): global DoF index of i-th local DoF on j-th element.
        % Negative value indicates corresponding base function takes negative sign.
        BC BC; % Boundary condition.
    end
    properties (Dependent)
        nLcDoF; % Number of local degrees of freedom.
        sElParm (1, 2); % Size of element parameter.
        nGlDoF; % Number of global degrees of freedom.
    end
    methods
        %% Constructor.
        function FES = FES(msh, fE, bC)
            arguments
                msh Msh;
                fE FE;
                bC BC = BC.empty;
            end
            switch msh.type
                case "D2T"
                    assert(ismember(fE.elem, ["D2T", "D2LR"]));
                case "D3T"
                    assert(ismember(fE.elem, ["D3T", "D3FR"]));
                otherwise
                    error("Unsupported mesh type.");
            end
            FES.msh = msh;
            FES.elem = fE.elem;
            FES.LcBase = fE.base;
            FES.map = fE.map;
            FES.RefBase = fE.RefBase;
            FES.RefKey = fE.RefKey;
            switch FES.elem
                case {"D2T", "D3T"}
                    % FES.ElParm = zeros(msh.dim, msh.elem.nNode, msh.nElem);
                    % for iElem = 1:msh.nElem
                    %     FES.ElParm(:, :, iElem) = msh.node.coord(:, msh.elem.node(:, iElem));
                    % end
                    FES.ElParm = reshape(msh.node.coord(:, msh.elem.node), msh.dim, msh.elem.nNode, msh.nElem);
                case "D2LR"
                    % FES.ElParm = zeros(msh.dim, msh.edge.nNode, msh.nEdge);
                    % for iEdge = 1:msh.nEdge
                    %     FES.ElParm(:, :, iEdge) = msh.node.coord(:, msh.edge.node(:, iEdge));
                    % end
                    FES.ElParm = reshape(msh.node.coord(:, msh.edge.node), msh.dim, msh.edge.nNode, msh.nEdge);
                case "D3FR"
                    % Trace space: "elements" are faces, parameterized by global face vertices.
                    FES.ElParm = reshape(msh.node.coord(:, msh.face.node), msh.dim, msh.face.nNode, msh.nFace);
            end
            [FES.GlDoFs, FES.Lc2Gl] = asmDoF(msh, fE.DoFs);
            if ~isempty(bC)
                FES.BC = bC.setDoF(FES.GlDoFs);
            else
                FES.BC = bC;
            end
        end
        % Get functions.
        function nLcDoF = get.nLcDoF(FES)
            nLcDoF = length(FES.LcBase);
        end
        function sElParm = get.sElParm(FES)
            sElParm = size(FES.ElParm);
            sElParm = sElParm(1:2);
        end
        function nGlDoF = get.nGlDoF(FES)
            nGlDoF = cumDoF(FES.GlDoFs);
        end
        function msh = getMsh(FESs)
            msh = FESs(1).msh;
            for iFES = 2:length(FESs)
                assert(isequal(FESs(iFES).msh, msh));
            end
        end
        %% Public functions.
        function cumDoF = cumDoF(FESs, iFES)
            % FES.cumDoF: cumulative number of DoFs up to the i-th FE space.
            if nargin == 1
                cumDoF = sum([FESs.nGlDoF]);
            elseif nargin == 2
                if iFES <= 0
                    cumDoF = 0;
                else
                    cumDoF = sum([FESs(1:iFES).nGlDoF]);
                end
            end
        end
        function fEF = proj(fES, fcn)
            % FES.proj: project function to finite element space.
            arguments (Input)
                fES FES;
                fcn Fcn;
            end
            arguments (Output)
                fEF FEF;
            end
            DoFVal = zeros(fES.nGlDoF, 1);
            for iDoF = 1:length(fES.GlDoFs)
                val = fES.GlDoFs(iDoF).eval(fcn, "valType", "num");
                DoFVal(fES.GlDoFs.cumDoF(iDoF - 1) + 1:fES.GlDoFs.cumDoF(iDoF)) = val(:);
            end
            fEF = FEF(fES, DoFVal);
        end
    end
    %% Static functions.
    methods (Static)
        function [permTab, sgnTab] = facePerm(LcDoF)
            % FES.facePerm: permutation of samples of DoF on face.
            % A face with permutation code c in an element has local nodes L = G(P(c, :)), where G are global face nodes
            % (see `Elem.facePerm`). Local sample i then corresponds to global sample permTab(i, c) with sign sgnTab(i, c):
            % - NdDoF: sample node with local barycentric coordinates (1 - a - b, a, b) w.r.t. L has global barycentric
            %   coordinates lambda_G(P(c, k)) = lambda_L(k); the set of sample nodes must be invariant under permutation.
            % - MoDoF: test function q_i in local coordinates equals sgnTab(i, c) * q_j in global coordinates, j = permTab(i, c);
            %   the set of test functions must be invariant (up to sign) under permutation.
            arguments (Input)
                LcDoF DoF;
            end
            arguments (Output)
                permTab (:, 6); % permTab(i, c): global sample corresponding to local sample i under code c.
                sgnTab (:, 6); % sgnTab(i, c): sign (+1 or -1).
            end
            assert(isequal(LcDoF.domn, "D3") && LcDoF.EntDim == 2);
            P = [1, 2, 3; 2, 3, 1; 3, 1, 2; 1, 3, 2; 3, 2, 1; 2, 1, 3];
            nSamp = LcDoF.nSamp;
            permTab = zeros(nSamp, 6);
            sgnTab = ones(nSamp, 6);
            switch class(LcDoF)
                case "NdDoF"
                    crd = double(LcDoF.coord);
                    for c = 1:6
                        barG = zeros(3, nSamp);
                        barG(P(c, :), :) = [1 - sum(crd, 1); crd];
                        for iSamp = 1:nSamp
                            jSamp = find(all(abs(crd - barG(2:3, iSamp)) < 1e-12, 1));
                            % Sample nodes on face must be invariant under permutation.
                            assert(isscalar(jSamp));
                            permTab(iSamp, c) = jSamp;
                        end
                    end
                case "MoDoF"
                    % Test functions are compared numerically at points of the reference face (no symmetry).
                    var = MshEnt.getInfo("D3R2").var;
                    pnt = [0.21, 0.13, 0.52, 0.34; 0.17, 0.61, 0.25, 0.08];
                    nPnt = size(pnt, 2);
                    barG = [1 - sum(pnt, 1); pnt];
                    tstFun = cell(1, nSamp);
                    tstVal = zeros(nSamp, nPnt);
                    for iSamp = 1:nSamp
                        assert(ismember(LcDoF.tst(iSamp).domn, ["VOID", "D3R2", "D3FR"]));
                        tstFun{iSamp} = matlabFunction(LcDoF.tst(iSamp).fun, "Vars", {var});
                        for iPnt = 1:nPnt
                            tstVal(iSamp, iPnt) = tstFun{iSamp}(pnt(:, iPnt));
                        end
                    end
                    tol = 1e-10 * max(1, max(abs(tstVal), [], "all"));
                    for c = 1:6
                        % Local coordinates (lambda_L(2), lambda_L(3)) = (lambda_G(P(c, 2)), lambda_G(P(c, 3))).
                        pntL = barG(P(c, 2:3), :);
                        for iSamp = 1:nSamp
                            tstL = zeros(1, nPnt);
                            for iPnt = 1:nPnt
                                tstL(iPnt) = tstFun{iSamp}(pntL(:, iPnt));
                            end
                            isFound = false;
                            for jSamp = 1:nSamp
                                for s = [1, -1]
                                    if all(abs(tstL - s * tstVal(jSamp, :)) < tol)
                                        permTab(iSamp, c) = jSamp;
                                        sgnTab(iSamp, c) = s;
                                        isFound = true;
                                        break;
                                    end
                                end
                                if isFound
                                    break;
                                end
                            end
                            % Test functions on face must be invariant (up to sign) under permutation.
                            assert(isFound);
                        end
                    end
            end
        end
    end
end
%% Local functions.
function [GlDoFs, Lc2Gl] = asmDoF(msh, LcDoFs)
    % asmDoF: assemble local DoFs to global DoFs.
    assert(ismember(msh.type, ["D2T", "D3T"]));
    GlDoFs = LcDoFs;
    if ismember(LcDoFs.getDomn, ["D2", "D3"])
        Lc2Gl = zeros(LcDoFs.cumDoF, msh.nElem);
    elseif ismember(LcDoFs.getDomn, ["D2R1", "D3R2"])
        % Trace space: "elements" are facets.
        Lc2Gl = zeros(LcDoFs.cumDoF, msh.nEnt(msh.dim - 1));
    else
        error("Unsupported DoF domain type.");
    end
    for iDoF = 1:length(LcDoFs)
        LcDoF = LcDoFs(iDoF); GlDoF = GlDoFs(iDoF);
        GlDoF.msh = msh;
        sgn = [];
        switch LcDoF.domn
            case {"D2", "D3"}
                switch LcDoF.EntDim
                    case 0
                        El_EntIdx = msh.elem.node(LcDoF.EntIdx, :);
                    case msh.dim
                        % Element.
                        El_EntIdx = 1:msh.nElem;
                    case 1
                        El_EntIdx = abs(msh.elem.edge(LcDoF.EntIdx, :));
                        sgn = sign(msh.elem.edge(LcDoF.EntIdx, :));
                    case 2
                        % Face (3D).
                        El_EntIdx = abs(msh.elem.face(LcDoF.EntIdx, :));
                        sgn = sign(msh.elem.face(LcDoF.EntIdx, :));
                        code = msh.elem.facePerm(LcDoF.EntIdx, :);
                end
                if LcDoF.share
                    GlDoF.EntIdx = unique(El_EntIdx)';
                    map = zeros(msh.nEnt(LcDoF.EntDim), 1); map(GlDoF.EntIdx) = 1:GlDoF.nEnt;
                    El_iEnt = reshape(map(El_EntIdx), size(El_EntIdx));
                else
                    GlDoF.EntIdx = El_EntIdx(:)';
                    El_iEnt = reshape(1:GlDoF.nEnt, size(El_EntIdx));
                end
            case {"D2R1", "D3R2"}
                assert(isequal(LcDoF.EntDim, msh.dim - 1));
                GlDoF.EntIdx = 1:msh.nEnt(msh.dim - 1);
                El_iEnt = 1:msh.nEnt(msh.dim - 1);
        end
        GlDoFs(iDoF) = GlDoF;
        if ismember(LcDoF.domn, ["D2", "D3"]) && LcDoF.EntDim == 1
            % Samples on edge are reversed if edge orientation in element is opposite to global orientation.
            for iSamp = 1:LcDoF.nSamp
                Lc2Gl(LcDoFs.sub2ind(iDoF, 1:LcDoF.nEnt, iSamp), :) = (GlDoFs.sub2ind(iDoF, El_iEnt, iSamp) .* (1 + sgn) / 2 + GlDoFs.sub2ind(iDoF, El_iEnt, GlDoF.nSamp - iSamp + 1) .* (1 - sgn) / 2) .* (sgn * LcDoF.orien + 1 * (1 - LcDoF.orien));
            end
        elseif isequal(LcDoF.domn, "D3") && LcDoF.EntDim == 2
            % Samples on face are permuted according to permutation code of face in element (see `FES.facePerm`);
            % sign of base function follows orientation of face in element.
            [permTab, sgnTab] = FES.facePerm(LcDoF);
            for iSamp = 1:LcDoF.nSamp
                GlSamp = reshape(permTab(iSamp, code(:)), size(code));
                GlSgn = reshape(sgnTab(iSamp, code(:)), size(code));
                Lc2Gl(LcDoFs.sub2ind(iDoF, 1:LcDoF.nEnt, iSamp), :) = GlDoFs.sub2ind(iDoF, El_iEnt, GlSamp) .* GlSgn .* (sgn * LcDoF.orien + 1 * (1 - LcDoF.orien));
            end
        else
            for iSamp = 1:LcDoF.nSamp
                Lc2Gl(LcDoFs.sub2ind(iDoF, 1:LcDoF.nEnt, iSamp), :) = GlDoFs.sub2ind(iDoF, El_iEnt, iSamp);
            end
        end
    end
end
