classdef BC
    % BC: boundary condition.
    properties
        fcn (1, :) Fcn; % Boundary condition function.
                        % Degree of freedom takes specified value at the boundary.
        % Mesh.
        node (1, :); % Index of boundary node.
        edge (1, :); % Index of boundary edge.
        face (1, :); % Index of boundary face (3D).
        DoF (1, :); % Index of single DoFs (in the FE space) taking specified values, e.g. to fix pressure at one DoF.
        % Degree of freedom.
        DoFIdx (1, :); % Index of boundary DoF.
        DoFVal (1, :); % Value of boundary DoF.
    end
    properties (Dependent)
        nDoF; % Number of boundary DoF.
    end
    methods
        % Constructor.
        function BC = BC(fcn, options)
            arguments
                fcn (1, :) Fcn;
                options.node (1, :) = [];
                options.edge (1, :) = [];
                options.face (1, :) = [];
                options.DoF (1, :) = [];
            end
            BC.fcn = fcn;
            BC.node = options.node;
            BC.edge = options.edge;
            BC.face = options.face;
            BC.DoF = options.DoF;
        end
        % Get functions.
        function nDoF = get.nDoF(BC)
            assert(length(BC.DoFIdx) == length(BC.DoFVal));
            nDoF = length(BC.DoFIdx);
        end
        % Public functions.
        function BC = setDoF(BC, DoFs)
            % BC.setDoF: set degree of freedom for boundary condition.
            if isscalar(BC.fcn)
                BC.fcn = repmat(BC.fcn, 1, length(DoFs));
            else
                assert(length(BC.fcn) == length(DoFs));
            end
            assert(ismember(DoFs.getDomn, ["D2", "D2R1", "D3", "D3R2"]));
            assert(ismember(DoFs.getMsh.type, ["D2T", "D3T"]));
            BC.DoFIdx = zeros(1, DoFs.cumDoF);
            BC.DoFVal = zeros(1, DoFs.cumDoF);
            for iDoF = 1:length(DoFs)
                DoF = DoFs(iDoF);
                assert(isscalar(BC.fcn(iDoF)));
                val = DoF.eval(BC.fcn(iDoF), "valType", "num", "rawEval", true);
                % Single DoFs of this group (`val(:)` is ordered as `DoFs.sub2ind`).
                idx = DoFs.sub2ind(iDoF);
                isSel = ismember(idx, BC.DoF);
                BC.DoFIdx(idx(isSel)) = 1;
                BC.DoFVal(idx(isSel)) = val(isSel);
                if ~isempty(BC.DoF) && isempty([BC.node, BC.edge, BC.face])
                    continue;
                end
                % Boundary entities of this group (rows of `val` and of `idx`).
                switch DoF.EntDim
                    case 0
                        assert(~isempty(BC.node));
                        BdEnt = BC.node;
                    case 1
                        assert(~isempty(BC.edge));
                        BdEnt = BC.edge;
                    case 2
                        % Face in 3D (element in 2D: no boundary condition).
                        if DoF.msh.dim ~= 3
                            continue;
                        end
                        assert(~isempty(BC.face));
                        BdEnt = BC.face;
                    otherwise
                        continue;
                end
                isEnt = ismember(DoF.EntIdx(:), BdEnt);
                idx = reshape(idx, DoF.sDoF);
                BC.DoFIdx(idx(isEnt, :)) = 1;
                BC.DoFVal(idx(isEnt, :)) = val(isEnt, :);
            end
            BC.DoFIdx = find(BC.DoFIdx);
            BC.DoFVal = BC.DoFVal(BC.DoFIdx);
        end
    end
end
