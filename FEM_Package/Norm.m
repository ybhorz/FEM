classdef Norm
    % Norm: norm of function.
    % Integration of `form(coef, fcn.dif(ord), pow)` over specific mesh entity.
    properties
        %% Integrated domain.
        msh Msh; % Mesh.
        EntDim {mustBeMember(EntDim, [1, 2])}; % Dimension of mesh entity.
        % EntDim = 1: line.
        % EntDim = 2: face.
        EntIdx (1, :); % Indices of mesh entities.
        %% Integrant.
        ord (:, :, :) % Order of derivative.
        pow; % Power of norm.
        iFcn; % Index of function.
        coef (1, :) Fcn; % Coefficient function.
        form; % Form of norm: function handle of `coef`, `fcn` and `pow`.
        fcnOpr {mustBeMember(fcnOpr, ["none", "jump"])} = "none"; % Operation on function.
        %% Integral method.
        GInt GInt; % Gauss integration.
    end
    methods
        % Constructor.
        function Norm = Norm(msh, EntDim, ord, options)
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [1, 2])};
                ord (:, :, :);
                options.EntIdx (1, :) = [];
                options.pow = 2;
                options.iFcn = 1;
                options.coef (1, :) Fcn = Fcn.cst(1);
                options.form = @(coef, fcn, pow) sum(coef .* abs(fcn) .^ pow);
                options.fcnOpr {mustBeMember(options.fcnOpr, ["none", "jump"])} = "none"
                options.GInt GInt = GInt.empty;
            end
            checkProp(msh.type, EntDim, options.coef, options.fcnOpr, options.GInt);
            Norm.msh = msh;
            Norm.EntDim = EntDim;
            Norm.ord = ord;
            if ~isempty(options.EntIdx)
                Norm.EntIdx = options.EntIdx;
            else
                Norm.EntIdx = 1:msh.nEnt(EntDim);
            end
            Norm.pow = options.pow;
            Norm.iFcn = options.iFcn;
            Norm.coef = options.coef;
            Norm.form = options.form;
            Norm.fcnOpr = options.fcnOpr;
            if ~isempty(options.GInt)
                Norm.GInt = options.GInt;
            else
                warning("No GInt specified.");
                switch EntDim
                    case 2
                        Norm.GInt = GInt("D2T", 1);
                    case 1
                        Norm.GInt = GInt("D2L", 1);
                end
            end
        end
        % Get functions.
        function msh = getMsh(Norms)
            msh = Norms(1).msh;
            for iNorm = 2:length(Norms)
                assert(isequal(Norms(iNorm).msh, msh));
            end
        end
    end
end
% Local functions.
function checkProp(mshType, EntDim, coef, fcnOpr, GInt)
    % checkProp: check validity of properties.
    assert(ismember(mshType, "D2T"));
    switch EntDim
        case 2
            assert(ismember(coef.getDomn, ["VOID", "D2", "D2T"]));
        case 1
            assert(ismember(coef.getDomn, ["VOID", "D2", "D2L"]));
    end
    switch EntDim
        case 2
            assert(ismember(fcnOpr, "none"));
        case 1
            assert(ismember(fcnOpr, ["none", "jump"]));
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
