classdef DoF
    % DoF: degree of freedom.
    % Degree of freedom is a functional that samples function on specific mesh entity.
    % One `DoF` represent a group of DoFs including `nSamp` samplers on each of `nEnt` mesh entities with dimension `EntDim`.
    properties
        domn {mustBeMember(domn, ["VOID", "D2", "D2R1"])} = "VOID"; % Domain of function.
        msh Msh; % Mesh.
        EntDim {mustBeMember(EntDim, [0, 1, 2])}; % Dimension of mesh entity.
        % EntDim = 0: point.
        % EntDim = 1: line.
        % EntDim = 2: face.
        EntIdx (1, :); % Indices of mesh entities.
        share logical; % Whether DoF is shared among mesh entities.
        orien logical; % Whether orientation of mesh entity is considered.
    end
    properties (Abstract)
        nEnt % Number of mesh entities.
        nSamp % Number of samplers
        nDoF % Number of degree of freedoms.
        sDoF (1, 2) % Size of DoF group: [nEnt, nSamp].
    end
    methods (Abstract)
        val = eval(DoF, fcn, options)
        % DoF.eval: evaluate DoF on function.
        % val(i,j): function's DoF value with respect to j-th sampler on i-th mesh entity.
    end
    methods
        %% Get functions.
        function domn = getDomn(DoFs)
            domn = DoFs(1).domn;
            for iDoF = 2:length(DoFs)
                assert(isequal(DoFs(iDoF).domn, domn));
            end
        end
        function msh = getMsh(DoFs)
            msh = DoFs(1).msh;
            for iDoF = 2:length(DoFs)
                assert(isequal(DoFs(iDoF).msh, msh));
            end
        end
        %% Public functions.
        function cumDoF = cumDoF(DoFs, iDoF)
            % DoF.cumDoF: cumulative number of DoFs to the i-th DoF group.
            if nargin == 1
                cumDoF = sum([DoFs.nDoF]);
            elseif nargin == 2
                if iDoF <= 0
                    cumDoF = 0;
                else
                    cumDoF = sum([DoFs(1:iDoF).nDoF]);
                end
            end
        end
        function DoFIdx = sub2ind(DoFs, iDoF, iEnt, iSamp)
            % DoF.sub2ind: convert subscript of DoF to index.
            if nargin == 2
                DoFIdx = cumDoF(DoFs, iDoF - 1) + 1:cumDoF(DoFs, iDoF);
            elseif nargin == 4
                DoFIdx = cumDoF(DoFs, iDoF - 1) + sub2ind(DoFs(iDoF).sDoF, iEnt, iSamp);
            end
        end
        function [iDoF, iEnt, iSamp] = ind2sub(DoFs, DoFIdx)
            % DoF.ind2sub: convert index of DoF to subscript.
            for iDoF = 1:length(DoFs)
                if DoFIdx <= cumDoF(DoFs, iDoF)
                    break;
                end
            end
            DoFIdx = DoFIdx - cumDoF(DoFs, iDoF - 1);
            [iEnt, iSamp] = ind2sub(DoFs(iDoF).sDoF, DoFIdx);
        end
    end
end
