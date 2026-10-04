classdef FE
    % FE: Finite Element.

    % Valid properties:
    % elem | DoFs.domn | DoFs.msh
    % -----|-----------|------------------
    % D2T  | D2        | MshEnt("D2T").msh
    % D2LR | D2R1      | MshEnt("D2LR").msh
    % D3T  | D3        | MshEnt("D3T").msh
    % D3FR | D3R2      | MshEnt("D3FR").msh

    % Mapping of base functions:
    % map      | elem     | DoFs
    % ---------|----------|-------------------------------------------------------------------------------------------
    % none     | all      | all
    % affine   | D2T D3T  | NdDoF with zero-order derivatives and constant coefficients
    % piolaDiv | D2T D3T  | MoDoF on facets of vector-valued function with zero-order derivatives and normal coefficient
    % piolaCurl| D2T D3T  | MoDoF on edges of vector-valued function with zero-order derivatives and tangent coefficient

    properties
        elem {mustBeMember(elem, ["VOID", "D2T", "D2LR", "D3T", "D3FR"])} = "VOID"; % Element.
        FS (:, :, :) sym; % Function space.
        % If dim(FS) = 2, function space consists of scalar/vector-valued functions spanned by FS(:,i).
        % If dim(FS) = 3, function space consists of matrix-valued functions spanned by FS(:,:,i).
        DoFs (1, :) DoF = NdDoF.empty; % Degrees of freedom.
        base (1, :) Fcn; % Base function.
        map {mustBeMember(map, ["none", "affine", "piolaDiv", "piolaCurl"])} = "none"; % Mapping of base functions.
        % none: invert DoF matrix on general element symbolically.
        % affine: invert DoF matrix on reference element, then compose base functions with transformation to reference element.
        % Remark: "affine" avoids symbolic inversion with element parameters, which is very slow in 3D
        % (e.g. 3D P2 element takes about 20 minutes with "none"). It is valid only for affine-equivalent elements,
        % i.e. DoFs are function values at nodes defined by relative coordinates (Lagrange, CR, DG, etc.).
        % piolaDiv: as "affine", but vector-valued base functions are mapped by contravariant Piola transformation
        % v(x) = B * v_ref(toRef(x)) / det(B), where x = B * x_ref + b. It preserves normal moments on facets
        % (int_F v.n q ds = int_F_ref v_ref.n_ref q ds_ref), hence is valid for H(div) elements (RT, BDM) whose DoFs are
        % such moments (`MoDoF` with `coef` = UNV). Point values of normal component (`NdDoF`) are not preserved.
        % piolaCurl: as "affine", but vector-valued base functions are mapped by covariant Piola transformation
        % v(x) = B^{-T} * v_ref(toRef(x)). It preserves tangential moments on edges
        % (int_e v.t q ds = int_e_ref v_ref.t_ref q ds_ref), hence is valid for H(curl) elements (Nedelec) whose DoFs are
        % such moments (`MoDoF` on edges with `coef` = UTV).
    end
    properties (Dependent)
        nDoF; % Number of degrees of freedom.
    end
    methods
        % Constructor.
        function FE = FE(elem, FS, DoFs, options)
            arguments
                elem {mustBeMember(elem, ["D2T", "D2LR", "D3T", "D3FR"])};
                FS (:, :, :) {mustBeA(FS, ["sym", "string"])};
                DoFs (1, :) DoF;
                options.map {mustBeMember(options.map, ["none", "affine", "piolaDiv", "piolaCurl"])} = "none";
            end
            if isstring(FS)
                FS = str2sym(FS);
            end
            checkProp(elem, FS, DoFs, options.map);
            FE.elem = elem;
            FE.FS = FS;
            FE.DoFs = DoFs;
            FE.map = options.map;
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
                    case "D3T"
                        if ndims(FE.FS) == 2
                            FSFcn(iFcn) = Fcn("D3", FE.FS(:, iFcn));
                        else
                            FSFcn(iFcn) = Fcn("D3", FE.FS(:, :, iFcn));
                        end
                    case "D3FR"
                        if ndims(FE.FS) == 2
                            FSFcn(iFcn) = Fcn("D3R2", FE.FS(:, iFcn));
                        else
                            FSFcn(iFcn) = Fcn("D3R2", FE.FS(:, :, iFcn));
                        end
                end
            end
            % Domain of base functions before mapping, function space and DoFs on it.
            DoFs = FE.DoFs;
            switch FE.map
                case "none"
                    BsDomn = FE.elem;
                    BsFS = FE.FS;
                case {"affine", "piolaDiv", "piolaCurl"}
                    % Reference element, e.g. "D2TR", "D3TR".
                    BsDomn = MshEnt.getInfo(FE.elem).dual;
                    BsFS = subs(FE.FS, MshEnt.getInfo(FE.elem).var, MshEnt.getInfo(BsDomn).var);
                    for iDoF = 1:length(DoFs)
                        DoFs(iDoF).msh = MshEnt(BsDomn).msh;
                    end
            end
            FSDoF = sym(zeros(FE.nDoF));
            for iDoF = 1:length(DoFs)
                for iFcn = 1:FE.nDoF
                    DoFVal = DoFs(iDoF).eval(FSFcn(iFcn));
                    FSDoF(DoFs.sub2ind(iDoF), iFcn) = DoFVal(:);
                end
            end
            I = eye(FE.nDoF);
            BsCoef = FSDoF \ I;
            base(1:FE.nDoF) = Fcn.cst(0);
            for iBase = 1:FE.nDoF
                if ndims(BsFS) == 2
                    base(iBase) = Fcn(BsDomn, BsFS * BsCoef(:, iBase)).simplify;
                else
                    base(iBase) = Fcn(BsDomn, sum(BsFS .* reshape(repmat(BsCoef(:, iBase).', [size(BsFS, 1) * size(BsFS, 2), 1]), size(BsFS)), 3)).simplify;
                end
            end
            switch FE.map
                case "affine"
                    % Compose with transformation to reference element.
                    base = base.tfm(FE.elem);
                case "piolaDiv"
                    % Contravariant Piola transformation: v = B * (v_ref o toRef) / det(B).
                    elTfm = Tfm(FE.elem);
                    B = jacobian(elTfm.toOrg, elTfm.refVar);
                    base = base.tfm(FE.elem);
                    for iBase = 1:FE.nDoF
                        base(iBase) = Fcn(FE.elem, B * base(iBase).fun / det(B));
                    end
                case "piolaCurl"
                    % Covariant Piola transformation: v = B^{-T} * (v_ref o toRef).
                    elTfm = Tfm(FE.elem);
                    BInvT = inv(jacobian(elTfm.toOrg, elTfm.refVar)).';
                    base = base.tfm(FE.elem);
                    for iBase = 1:FE.nDoF
                        base(iBase) = Fcn(FE.elem, BInvT * base(iBase).fun);
                    end
            end
        end
    end
end
% Local functions.
function checkProp(elem, FS, DoFs, map)
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
        case "D3T"
            assert(isequal(DoFs.getDomn, "D3"));
            assert(isequal(DoFs.getMsh, MshEnt("D3T").msh));
        case "D3FR"
            assert(isequal(DoFs.getDomn, "D3R2"));
            assert(isequal(DoFs.getMsh, MshEnt("D3FR").msh));
    end
    switch map
        case "affine"
            assert(ismember(elem, ["D2T", "D3T"]));
            for iDoF = 1:length(DoFs)
                assert(isa(DoFs(iDoF), "NdDoF"));
                assert(all(~(DoFs(iDoF).ord(:) > 0)));
                assert(isequal(DoFs(iDoF).coef.getDomn, "VOID"));
            end
        case "piolaDiv"
            assert(ismember(elem, ["D2T", "D3T"]));
            dim = MshEnt(elem).dim;
            assert(ismatrix(FS) && size(FS, 1) == dim);
            for iDoF = 1:length(DoFs)
                assert(isa(DoFs(iDoF), "MoDoF"));
                assert(DoFs(iDoF).EntDim == dim - 1);
                assert(all(~(DoFs(iDoF).ord(:) > 0)));
                assert(isequal(DoFs(iDoF).coef.getDomn, DoFs(iDoF).msh.FtDomn));
            end
        case "piolaCurl"
            assert(ismember(elem, ["D2T", "D3T"]));
            dim = MshEnt(elem).dim;
            assert(ismatrix(FS) && size(FS, 1) == dim);
            % Domain of edge, e.g. "D2L", "D3L".
            EgDomn = replace(string(elem), "T", "L");
            for iDoF = 1:length(DoFs)
                assert(isa(DoFs(iDoF), "MoDoF"));
                assert(DoFs(iDoF).EntDim == 1);
                assert(all(~(DoFs(iDoF).ord(:) > 0)));
                assert(isequal(DoFs(iDoF).coef.getDomn, EgDomn));
            end
    end
end
