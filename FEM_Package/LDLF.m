classdef LDLF < DLF
    % LDLF: linearized double linear functional.
    % Integration of `form(coef, pre.dif(preOrd), trl.dif(trlOrd), tst.dif(tstOrd))` over specific mesh entity.
    properties
        preOrd (:, :, :); % Order of derivative of previous solution.
        iPre; % Index of previous solution.
    end
    methods
        % Constructor.
        function LDLF = LDLF(msh, EntDim, coef, preOrd, trlOrd, tstOrd, options)
            arguments
                msh Msh;
                EntDim {mustBeMember(EntDim, [1, 2])};
                coef (1, :) Fcn;
                preOrd (:, :, :);
                trlOrd (:, :, :);
                tstOrd (:, :, :);
                options.EntIdx (1, :) = [];
                options.iPre = 1;
                options.iTrl = 1;
                options.iTst = 1;
                options.form = @(coef, pre, trl, tst) sum(coef .* pre .* trl .* tst);
                options.GInt GInt = GInt.empty;
            end
            LDLF = LDLF@DLF(msh, EntDim, coef, trlOrd, tstOrd, "EntIdx", options.EntIdx, "iTrl", options.iTrl, "iTst", options.iTst, "form", options.form, "GInt", options.GInt);
            LDLF.preOrd = preOrd;
            LDLF.iPre = options.iPre;
        end
    end
end
