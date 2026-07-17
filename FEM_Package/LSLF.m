classdef LSLF < SLF
    % LSLF: linearized single linear functional.
    % Integration of `form(load, pre.dif(preOrd), tst.dif(tstOrd))` over specific mesh entity.
    properties
        preOrd (:, :, :); % Order of derivative of previous solution.
        iPre; % Index of previous solution.
    end
    methods
        % Constructor.
        function LSLF = LSLF(msh, EntDim, load, preOrd, tstOrd, options)
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [1, 2])};
                load (1, :) Fcn;
                preOrd (:, :, :);
                tstOrd (:, :, :);
                options.EntIdx (1, :) = [];
                options.iPre = 1;
                options.iTst = 1;
                options.form = @(load, pre, tst) sum(load .* pre .* tst);
                options.GInt GInt = GInt.empty;
            end
            LSLF = LSLF@SLF(msh, EntDim, load, tstOrd, "EntIdx", options.EntIdx, "iTst", options.iTst, "form", options.form, "GInt", options.GInt);
            LSLF.preOrd = preOrd;
            LSLF.iPre = options.iPre;
        end
    end
end