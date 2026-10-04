classdef SLF
    % SLF: single linear functional.
    % Integration of `form(load, tst.dif(tstOrd))` over specific mesh entity.
    properties
        %% Integrated domain.
        msh Msh; % Mesh.
        EntDim {mustBeMember(EntDim, [1, 2, 3])}; % Dimension of mesh entity.
        % EntDim = 1: line.
        % EntDim = 2: face.
        % EntDim = 3: volume.
        % `EntDim` is either `msh.dim` (element) or `msh.dim - 1` (facet).
        EntIdx (1, :); % Indices of mesh entities.
        %% Integrant.
        load (1, :) Fcn; % Load function.
        tstOrd (:, :, :) % Order of derivative of test function.
        iTst; % Index of test function.
        % When `EntDim` = `msh.dim - 1`, positive/negative `iTst` indicates taking function trace from positive/negative connected element of facet.
        form; % Form of SLF: function handle of `load` and `tst`.
        %% Integral method.
        GInt GInt; % Gauss integration.
    end
    methods
        % Constructor.
        function SLF = SLF(msh, EntDim, load, tstOrd, options)
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [1, 2, 3])};
                load (1, :) Fcn;
                tstOrd (:, :, :);
                options.EntIdx (1, :) = [];
                options.iTst = 1;
                options.form = @(load, tst) sum(load .* tst);
                options.GInt GInt = GInt.empty;
            end
            checkProp(msh, EntDim, load, options.GInt);
            SLF.msh = msh;
            SLF.EntDim = EntDim;
            SLF.load = load;
            SLF.tstOrd = tstOrd;
            if ~isempty(options.EntIdx)
                SLF.EntIdx = options.EntIdx;
            else
                SLF.EntIdx = 1:msh.nEnt(EntDim);
            end
            SLF.iTst = options.iTst;
            SLF.form = options.form;
            if ~isempty(options.GInt)
                SLF.GInt = options.GInt;
            else
                warning("No GInt specified.");
                switch EntDim
                    case msh.dim
                        SLF.GInt = GInt(msh.ElDomn, 1);
                    case msh.dim - 1
                        SLF.GInt = GInt(msh.FtDomn, 1);
                end
            end
        end
        % Get functions.
        function msh = getMsh(SLFs)
            msh = SLFs(1).msh;
            for iSLF = 2:length(SLFs)
                assert(isequal(SLFs(iSLF).msh, msh));
            end
        end
    end
end
% Local functions.
function checkProp(msh, EntDim, load, GInt)
    % checkProp: check validity of properties.
    assert(ismember(msh.type, ["D2T", "D3T"]));
    assert(ismember(EntDim, [msh.dim - 1, msh.dim]));
    switch EntDim
        case msh.dim
            assert(ismember(load.getDomn, ["VOID", msh.domn, msh.ElDomn]));
        case msh.dim - 1
            assert(ismember(load.getDomn, ["VOID", msh.domn, msh.FtDomn]));
    end
    if ~isempty(GInt)
        switch EntDim
            case msh.dim
                assert(isequal(GInt.domn, msh.ElDomn));
            case msh.dim - 1
                assert(isequal(GInt.domn, msh.FtDomn));
        end
    end
end
