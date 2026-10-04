function val = mapVal(RefBase, map, ord, X, ElParm, RefKey)
    % mapVal: values of (derivatives of) mapped base functions at physical points, evaluated numerically.
    % Base functions of an FE with map "affine", "piolaDiv" or "piolaCurl" (see `FE`) are v = M * (v_ref o toRef),
    % where x = B * x_ref + b maps the reference element (vertices 0, e_1, ..., e_d) to the element and M = I, B / det(B)
    % or B^{-T}. Their derivatives follow from the chain rule, d/dx_j = sum_k (B^{-1})_{kj} d/dx_ref_k, so only the
    % (small) reference base functions are differentiated and evaluated, instead of the mapped functions whose symbolic
    % expressions contain all element parameters.
    % val(k, i, p, e): k-th component of derivative `ord` (see `Fcn.dif`) of i-th base function at X(:, p, e).
    arguments (Input)
        RefBase (1, :) Fcn; % Base functions on reference element.
        map {mustBeMember(map, ["affine", "piolaDiv", "piolaCurl"])};
        ord (:, :, :); % Order of derivative.
        X (:, :, :); % Physical points, d x nPnt x nEnt.
        ElParm (:, :, :); % Vertices of elements, d x (d + 1) x nEnt.
        RefKey (1, 1) string = ""; % Key of `RefBase` (see `FE`): function handles of derivatives of reference base
        % functions are cached by key across calls (not cached if empty).
    end
    arguments (Output)
        val (:, :, :, :); % nOut x nBase x nPnt x nEnt.
    end
    [dim, nPnt, nEnt] = size(X);
    nBase = length(RefBase);
    sFun = size(RefBase(1).fun);
    nFun = prod(sFun);
    assert(size(ord, 1) == dim && size(ord, 2) == nFun);
    assert(size(ElParm, 3) == nEnt);
    % Element map and reference points.
    b = ElParm(:, 1, :);
    B = ElParm(:, 2:end, :) - b;
    BInv = pageinv(B);
    XRef = pagemtimes(BInv, X - b);
    switch map
        case "affine"
            M = repmat(eye(nFun), 1, 1, nEnt);
        case "piolaDiv"
            assert(nFun == dim);
            M = B ./ detVec(B);
        case "piolaCurl"
            assert(nFun == dim);
            M = pagetranspose(BInv);
    end
    % Output entries: (iFun, iOut) with derivative ord(:, iFun, iOut) of component iFun (see `Fcn.dif`).
    nOutSlc = size(ord, 3);
    if nOutSlc > 1
        assert(nFun == 1 || sFun(2) == 1);
    end
    val = zeros(nFun * nOutSlc, nBase, nPnt, nEnt);
    refVal = containers.Map(); % Values of derivatives of reference base functions at XRef, by multi-index.
    for iOut = 1:nOutSlc
        for iFun = 1:nFun
            alpha = ord(:, iFun, iOut);
            if ~all(alpha >= 0)
                continue;
            end
            % Directions of derivative, e.g. alpha = [1; 0; 2] -> [1, 3, 3].
            dirs = repelem(1:dim, alpha);
            n = length(dirs);
            % D(:, i, p, e): derivative of (all components of) v_ref o toRef, summed over sequences of reference
            % directions l_1, ..., l_n with coefficient prod_k (B^{-1})_{l_k, dirs(k)}.
            D = zeros(nFun, nBase, nPnt, nEnt);
            for iSeq = 1:dim ^ n
                seq = seqIdx(iSeq, dim, n);
                coef = ones(1, 1, 1, nEnt);
                for k = 1:n
                    coef = coef .* reshape(BInv(seq(k), dirs(k), :), 1, 1, 1, nEnt);
                end
                beta = accumarray(seq(:), 1, [dim, 1]);
                D = D + coef .* refDif(beta);
            end
            % Mapped component iFun: sum_i M(iFun, i) * D(i, ...).
            Mrow = reshape(M(iFun, :, :), nFun, 1, 1, nEnt);
            val(iFun + (iOut - 1) * nFun, :, :, :) = sum(Mrow .* D, 1);
        end
    end
    % Nested functions.
    function V = refDif(beta)
        % refDif: derivative beta of all components of reference base functions at XRef, nFun x nBase x nPnt x nEnt.
        key = char(mat2str(beta'));
        if isKey(refVal, key)
            V = refVal(key);
            return;
        end
        funH = refHandle(RefBase, beta, RefKey);
        out = cell(1, nFun * nBase);
        [out{:}] = funH(reshape(XRef, dim, nPnt * nEnt));
        V = zeros(nFun * nBase, nPnt * nEnt);
        for jOut = 1:nFun * nBase
            V(jOut, :) = out{jOut};
        end
        V = reshape(V, nFun, nBase, nPnt, nEnt);
        refVal(key) = V;
    end
end
%% Local functions.
function funH = refHandle(RefBase, beta, RefKey)
    % refHandle: function handle of derivative beta of all components of reference base functions (one output per
    % component and base function), cached by `RefKey` and `beta` across calls.
    persistent cache
    if isempty(cache)
        cache = containers.Map();
    end
    key = char(RefKey + "|" + mat2str(beta'));
    if strlength(RefKey) > 0 && isKey(cache, key)
        funH = cache(key);
        return;
    end
    var = RefBase(1).var;
    nBase = length(RefBase);
    nFun = numel(RefBase(1).fun);
    funs = cell(1, nFun * nBase);
    for iBase = 1:nBase
        fun = RefBase(iBase).fun;
        for jFun = 1:nFun
            f = fun(jFun);
            for iVar = 1:length(var)
                f = diff(f, var(iVar), beta(iVar));
            end
            funs{jFun + (iBase - 1) * nFun} = f;
        end
    end
    funH = matlabFunction(funs{:}, "Vars", {var});
    if strlength(RefKey) > 0
        cache(key) = funH;
    end
end
function seq = seqIdx(iSeq, dim, n)
    % seqIdx: iSeq-th sequence of n directions in 1:dim (lexicographic).
    seq = zeros(1, n);
    r = iSeq - 1;
    for k = n:-1:1
        seq(k) = mod(r, dim) + 1;
        r = floor(r / dim);
    end
end
function d = detVec(B)
    % detVec: determinants of a batch of 2 x 2 or 3 x 3 matrices, 1 x 1 x nEnt.
    switch size(B, 1)
        case 2
            d = B(1, 1, :) .* B(2, 2, :) - B(1, 2, :) .* B(2, 1, :);
        case 3
            d = B(1, 1, :) .* (B(2, 2, :) .* B(3, 3, :) - B(2, 3, :) .* B(3, 2, :)) ...
                - B(1, 2, :) .* (B(2, 1, :) .* B(3, 3, :) - B(2, 3, :) .* B(3, 1, :)) ...
                + B(1, 3, :) .* (B(2, 1, :) .* B(3, 2, :) - B(2, 2, :) .* B(3, 1, :));
    end
end
