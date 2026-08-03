classdef DLF
    % DLF: double linear functional.
    % Integration of `form(coef, trl.dif(trlOrd), tst.dif(tstOrd))` over specific mesh entity.
    properties
        %% Integrated domain.
        msh Msh; % Mesh.
        EntDim {mustBeMember(EntDim, [1, 2])}; % Dimension of mesh entity.
        % EntDim = 1: line.
        % EntDim = 2: face.
        EntIdx (1, :); % Indices of mesh entities.
        %% Integrant.
        coef (1, :) Fcn; % Coefficient function.
        trlOrd (:, :, :) % Order of derivative of trial function.
        tstOrd (:, :, :) % Order of derivative of test function.
        iTrl; % Index of trial function.
        iTst; % Index of test function.
        % When `EntDim` = 1, positive/negative `iTrl/iTst` indicates taking function trace from positive/negative connected element of edge.
        form; % Form of DLF: function handle of `coef`, `trl`, and `tst`.
        %% Integral method.
        GInt GInt; % Gauss integration.
    end
    methods
        % Constructor.
        function DLF = DLF(msh, EntDim, coef, trlOrd, tstOrd, options)
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [1, 2])};
                coef (1, :) Fcn;
                trlOrd (:, :, :);
                tstOrd (:, :, :);
                options.EntIdx (1, :) = [];
                options.iTrl = 1;
                options.iTst = 1;
                options.form = @(coef, trl, tst) sum(coef .* trl .* tst);
                options.GInt GInt = GInt.empty;
            end
            checkProp(msh.type, EntDim, coef, options.GInt);
            DLF.msh = msh;
            DLF.EntDim = EntDim;
            DLF.coef = coef;
            DLF.trlOrd = trlOrd;
            DLF.tstOrd = tstOrd;
            if ~isempty(options.EntIdx)
                DLF.EntIdx = options.EntIdx;
            else
                DLF.EntIdx = 1:msh.nEnt(EntDim);
            end
            DLF.iTrl = options.iTrl;
            DLF.iTst = options.iTst;
            DLF.form = options.form;
            if ~isempty(options.GInt)
                DLF.GInt = options.GInt;
            else
                warning("No GInt specified.");
                switch EntDim
                    case 2
                        DLF.GInt = GInt("D2T", 1);
                    case 1
                        DLF.GInt = GInt("D2L", 1);
                end
            end
        end
        % Get functions.
        function msh = getMsh(DLFs)
            msh = DLFs(1).msh;
            for iDLF = 2:length(DLFs)
                assert(isequal(DLFs(iDLF).msh, msh));
            end
        end
    end
    % Static functions.
    methods (Static)
        function DLFs = interface(msh, EntDim, coef, trlOrd, tstOrd, options)
            % DLF.interface: double linear functional on interface.
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, 1)};
                coef (1, :) Fcn;
                trlOrd (:, :, :);
                tstOrd (:, :, :);
                options.EntIdx (1, :) = [];
                options.iTrl = 1;
                options.iTst = 1;
                options.trlOpr {mustBeMember(options.trlOpr, ["none", "aver", "jump"])} = "none"; % Operation on trial function.
                options.tstOpr {mustBeMember(options.tstOpr, ["none", "aver", "jump"])} = "none"; % Operation on test function.
                options.form = @(coef, trl, tst) sum(coef .* trl .* tst);
                options.GInt GInt = GInt.empty;
            end
            checkProp(msh.type, EntDim, coef, options.GInt);
            switch options.trlOpr
                case "none"
                    nRow = 1; trlSgn = 1; trlWgt = 1;
                case "aver"
                    nRow = 2; trlSgn = [1, -1]; trlWgt = [1, 1] / 2;
                case "jump"
                    nRow = 2; trlSgn = [1, -1]; trlWgt = [1, -1];
            end
            switch options.tstOpr
                case "none"
                    nCol = 1; tstSgn = 1; tstWgt = 1;
                case "aver"
                    nCol = 2; tstSgn = [1, -1]; tstWgt = [1, 1] / 2;
                case "jump"
                    nCol = 2; tstSgn = [1, -1]; tstWgt = [1, -1];
            end
            DLFs(1:nRow, 1:nCol) = DLF(msh, EntDim, coef, trlOrd, tstOrd, "EntIdx", options.EntIdx, "iTrl", options.iTrl, "iTst", options.iTst, "form", options.form, "GInt", options.GInt);
            baseForm = options.form;
            for iRow = 1:nRow
                for iCol = 1:nCol
                    DLFs(iRow, iCol).iTrl = DLFs(iRow, iCol).iTrl * trlSgn(iRow);
                    DLFs(iRow, iCol).iTst = DLFs(iRow, iCol).iTst * tstSgn(iCol);
                    DLFs(iRow, iCol).form = @(varargin) baseForm(varargin{:}) * trlWgt(iRow) * tstWgt(iCol);
                end
            end
            DLFs = DLFs(:)';
        end
    end
end
% Local functions.
function checkProp(mshType, EntDim, coef, GInt)
    % checkProp: check validity of properties.
    assert(ismember(mshType, "D2T"));
    switch EntDim
        case 2
            assert(ismember(coef.getDomn, ["VOID", "D2", "D2T"]));
        case 1
            assert(ismember(coef.getDomn, ["VOID", "D2", "D2L"]));
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
