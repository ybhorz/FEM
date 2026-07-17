classdef FE
    % FE: Finite Element.

    % Valid properties:
    % elem | DoFs.domn | DoFs.msh
    % -----|-----------|------------------
    % D2T  | D2        | MshEnt("D2T").msh
    % D2LR | D2R1      | MshEnt("D2LR").msh

    properties
        elem {mustBeMember(elem, ["VOID", "D2T", "D2LR"])} = "VOID"; % Element.
        FS (:, :, :) sym; % Function space.
        % If dim(FS) = 2, function space consists of scalar/vector-valued functions spanned by FS(:,i).
        % If dim(FS) = 3, function space consists of matrix-valued functions spanned by FS(:,:,i).
        DoFs (1, :) DoF = NdDoF.empty; % Degree of freedoms.
        base (1, :) Fcn; % Base function.
    end
    properties (Dependent)
        nDoF; % Number of degree of freedoms.
    end
    methods
        % Constructor.
        function FE = FE(elem, FS, DoFs)
            arguments
                elem {mustBeMember(elem, ["D2T", "D2LR"])};
                FS (:, :, :) {mustBeA(FS, ["sym", "string"])};
                DoFs (1, :) DoF;
            end
            if isstring(FS)
                FS = str2sym(FS);
            end
            checkProp(elem, FS, DoFs);
            FE.elem = elem;
            FE.FS = FS;
            FE.DoFs = DoFs;
            FE.base = genBase(FE);
        end
        % Get functions.
        function nDoF = get.nDoF(FE)
            nDoF = cumDoF(FE.DoFs);
        end
    end
    % Static functions.
    methods (Static)
        function FS = repFS(FS, sz)
            % repFS: replicate function space.
            arguments
                FS (:, :) {mustBeA(FS, ["sym", "string"])};
                sz (1, 2);
            end
            if isstring(FS)
                FS = str2sym(FS);
            end
            assert(isrow(FS));
            nFS = length(FS);
            temp = FS;
            if sz(2) == 1
                FS = sym(zeros(sz(1), nFS * sz(1)));
                for i = 1:sz(1)
                    FS(i, (i - 1) * nFS + (1:nFS)) = temp;
                end
            else
                FS = sym(zeros(sz(1), sz(2), nFS * sz(1) * sz(2)));
                for j = 1:sz(2)
                    for i = 1:sz(1)
                        FS(i, j, (j - 1) * nFS * sz(1) + (i - 1) * nFS + (1:nFS)) = temp;
                    end
                end
            end
        end
    end
    % Private functions.
    methods (Access = private)
        function base = genBase(FE)
            % FE.genBase: generate base functions.
            FSFcn(1:FE.nDoF) = Fcn.cst(0);
            for iFcn = 1:FE.nDoF
                switch FE.elem
                    case "D2T"
                        if ndims(FE.FS) == 2
                            FSFcn(iFcn) = Fcn("D2", FE.FS(:, iFcn));
                        else
                            FSFcn(iFcn) = Fcn("D2", FE.FS(:, :, iFcn));
                        end
                    case "D2LR"
                        if ndims(FE.FS) == 2
                            FSFcn(iFcn) = Fcn("D2R1", FE.FS(:, iFcn));
                        else
                            FSFcn(iFcn) = Fcn("D2R1", FE.FS(:, :, iFcn));
                        end
                end
            end
            FSDoF = sym(zeros(FE.nDoF));
            for iDoF = 1:length(FE.DoFs)
                for iFcn = 1:FE.nDoF
                    DoFVal = FE.DoFs(iDoF).eval(FSFcn(iFcn));
                    FSDoF(FE.DoFs.sub2ind(iDoF), iFcn) = DoFVal(:);
                end
            end
            I = eye(FE.nDoF);
            BsCoef = FSDoF \ I;
            base(1:FE.nDoF) = Fcn.cst(0);
            for iBase = 1:FE.nDoF
                if ndims(FE.FS) == 2
                    base(iBase) = Fcn(FE.elem, FE.FS * BsCoef(:, iBase)).simplify;
                else
                    base(iBase) = Fcn(FE.elem, sum(FE.FS .* reshape(repmat(BsCoef(:, iBase).', [size(FE.FS, 1) * size(FE.FS, 2), 1]), size(FE.FS)), 3)).simplify;
                end
            end
        end
    end
end
% Local functions.
function checkProp(elem, FS, DoFs)
    % checkProp: check validity of properties.
    assert(size(FS, ndims(FS)) == cumDoF(DoFs));
    if ndims(FS) == 2
        assert(all(size(FS, 1) == [DoFs.LFun]));
    else
        assert(all(size(FS, 1) * size(FS, 2) == [DoFs.LFun]));
    end
    switch elem
        case "D2T"
            assert(isequal(DoFs.getDomn, "D2"));
            assert(isequal(DoFs.getMsh, MshEnt("D2T").msh));
        case "D2LR"
            assert(isequal(DoFs.getDomn, "D2R1"));
            assert(isequal(DoFs.getMsh, MshEnt("D2LR").msh));
    end
end
