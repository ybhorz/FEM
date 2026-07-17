classdef BC
    % BC: boundary condition.
    properties
        fcn (1, :) Fcn; % Boundary function.
        % Mesh.
        node (1, :); % Index of boundary node.
        edge (1, :); % Index of boundary edge.
        % Degree of freedom.
        DoFIdx (1, :); % Index of boundary degree of freedom.
        DoFVal (1, :); % Value of boundary degree of freedom.
    end
    methods
        % Constructor.
        function BC = BC(fcn, options)
            arguments
                fcn (1, :) Fcn;
                options.node (1, :) = [];
                options.edge (1, :) = [];
            end
            BC.fcn = fcn;
            BC.node = options.node;
            BC.edge = options.edge;
        end
        % Public functions.
        function BC = setDoF(BC, DoFs)
            % BC.setDoF: set degree of freedom for boundary condition.
            if isscalar(BC.fcn)
                BC.fcn = repmat(BC.fcn, 1, length(DoFs));
            else
                assert(length(BC.fcn) == length(DoFs));
            end
            assert(ismember(DoFs.getDomn, ["D2", "D2R1"]));
            assert(ismember(DoFs.getMsh.type, "D2T"));
            BC.DoFIdx = zeros(1, DoFs.cumDoF);
            BC.DoFVal = zeros(1, DoFs.cumDoF);
            for iDoF = 1:length(DoFs)
                DoF = DoFs(iDoF);
                val = DoF.eval(BC.fcn(iDoF), "valType", "num", "rawEval", true);
                switch DoF.EntDim
                    case 0
                        assert(~isempty(BC.node));
                        for iEnt = 1:length(DoF.EntIdx)
                            iNode = DoF.EntIdx(iEnt);
                            if ismember(iNode, BC.node)
                                idx = DoFs.sub2ind(iDoF, iEnt, 1:DoF.nSamp);
                                BC.DoFIdx(idx) = 1;
                                BC.DoFVal(idx) = val(iEnt, :);
                            end
                        end
                    case 1
                        assert(~isempty(BC.edge));
                        for iEnt = 1:length(DoF.EntIdx)
                            iEdge = DoF.EntIdx(iEnt);
                            if ismember(iEdge, BC.edge)
                                idx = DoFs.sub2ind(iDoF, iEnt, 1:DoF.nSamp);
                                BC.DoFIdx(idx) = 1;
                                BC.DoFVal(idx) = val(iEnt, :);
                            end
                        end
                end
            end
            BC.DoFIdx = find(BC.DoFIdx);
            BC.DoFVal = BC.DoFVal(BC.DoFIdx);
        end
    end
end
