classdef FES
    % FES: Finite Element space.
    properties
        msh Msh; % Mesh.
        elem {mustBeMember(elem, ["VOID", "D2T", "D2LR"])} = "VOID"; % Element type.
        LcBase (1, :) Fcn; % Local base function.
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
            assert(ismember(msh.type, "D2T"));
            assert(ismember(fE.elem, ["D2T", "D2LR"]));
            FES.msh = msh;
            FES.elem = fE.elem;
            FES.LcBase = fE.base;
            switch FES.elem
                case "D2T"
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
end
%% Local functions.
function [GlDoFs, Lc2Gl] = asmDoF(msh, LcDoFs)
    % asmDoF: assemble local DoFs to global DoFs.
    assert(ismember(msh.type, "D2T"));
    GlDoFs = LcDoFs;
    if ismember(LcDoFs.getDomn, "D2")
        Lc2Gl = zeros(LcDoFs.cumDoF, msh.nElem);
    elseif ismember(LcDoFs.getDomn, "D2R1")
        Lc2Gl = zeros(LcDoFs.cumDoF, msh.nEdge);
    else
        error("Unsupported DoF domain type.");
    end
    for iDoF = 1:length(LcDoFs)
        LcDoF = LcDoFs(iDoF); GlDoF = GlDoFs(iDoF);
        GlDoF.msh = msh;
        switch LcDoF.domn
            case "D2"
                switch LcDoF.EntDim
                    case 0
                        switch LcDoF.share
                            case true
                                El_EntIdx = msh.elem.node(LcDoF.EntIdx, :);
                                GlDoF.EntIdx = unique(El_EntIdx);
                                map = zeros(msh.nNode, 1); map(GlDoF.EntIdx) = 1:GlDoF.nEnt;
                                El_iEnt = reshape(map(El_EntIdx), size(El_EntIdx));
                            case false
                                El_EntIdx = msh.elem.node(LcDoF.EntIdx, :);
                                GlDoF.EntIdx = El_EntIdx(:);
                                El_iEnt = reshape(1:GlDoF.nEnt, size(El_EntIdx));
                        end
                    case 1
                        switch LcDoF.share
                            case true
                                El_EntIdx = abs(msh.elem.edge(LcDoF.EntIdx, :));
                                sgn = sign(msh.elem.edge(LcDoF.EntIdx, :));
                                GlDoF.EntIdx = unique(El_EntIdx);
                                map = zeros(msh.nEdge, 1); map(GlDoF.EntIdx) = 1:GlDoF.nEnt;
                                El_iEnt = reshape(map(El_EntIdx), size(El_EntIdx));
                            case false
                                El_EntIdx = abs(msh.elem.edge(LcDoF.EntIdx, :));
                                sgn = sign(msh.elem.edge(LcDoF.EntIdx, :));
                                GlDoF.EntIdx = El_EntIdx(:);
                                El_iEnt = reshape(1:GlDoF.nEnt, size(El_EntIdx));
                        end
                    case 2
                        GlDoF.EntIdx = 1:msh.nElem;
                        El_iEnt = 1:msh.nElem;
                end
            case "D2R1"
                assert(isequal(LcDoF.EntDim, 1));
                GlDoF.EntIdx = 1:msh.nEdge;
                El_iEnt = 1:msh.nEdge;
        end
        GlDoFs(iDoF) = GlDoF;
        if isequal(LcDoF.domn, "D2") && LcDoF.EntDim == 1
            for iSamp = 1:LcDoF.nSamp
                Lc2Gl(LcDoFs.sub2ind(iDoF, 1:LcDoF.nEnt, iSamp), :) = (GlDoFs.sub2ind(iDoF, El_iEnt, iSamp) .* (1 + sgn) / 2 + GlDoFs.sub2ind(iDoF, El_iEnt, GlDoF.nSamp - iSamp + 1) .* (1 - sgn) / 2) .* (sgn * LcDoF.orien + 1 * (1 - LcDoF.orien));
            end
        else
            for iSamp = 1:LcDoF.nSamp
                Lc2Gl(LcDoFs.sub2ind(iDoF, 1:LcDoF.nEnt, iSamp), :) = GlDoFs.sub2ind(iDoF, El_iEnt, iSamp);
            end
        end
    end
end
