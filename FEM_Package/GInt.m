classdef GInt
    % GInt: Gauss integration.
    properties
        domn {mustBeMember(domn, ["VOID", "D2T", "D2L", "D3T", "D3F", "D3L"])} = "VOID"; % Domain of integration.
        ord {mustBeMember(ord, [0, 1, 2, 3, 4, 5, 6, 7])}; % Order of integration.
        % - algebraic precision: the rule integrates polynomials up to this degree exactly.
        pnt (:, :); % Coordinates of Gauss points in reference domain (by column).
        wgt (1, :); % Weights of Gauss points in reference domain.
    end
    properties (Access = private)
        orgTfm; % Transformation to original domain.
        Jac; % Jacobian of transformation.
        sParm; % Size of parameter of integrated domain.
    end
    methods
        % Constructor.
        function GInt = GInt(domn, ord)
            arguments
                domn {mustBeMember(domn, ["D2T", "D2L", "D3T", "D3F", "D3L"])};
                ord {mustBeMember(ord, [0, 1, 2, 3, 4, 5, 6, 7])};
            end
            GInt.domn = domn;
            GInt.ord = ord;
            GInt.sParm = MshEnt(domn).sParm;
            switch domn
                case {"D2T", "D3F"}
                    assert(isequal(MshEnt(domn).dual.node.coord, sym([0, 1, 0; 0, 0, 1])));
                    switch ord
                        case {0, 1}
                            GInt.pnt = [1/3; 1/3];
                            GInt.wgt = 1/2;
                        case 2
                            GInt.pnt = [1/6, 2/3, 1/6;
                                        1/6, 1/6, 2/3];
                            GInt.wgt = [1/6, 1/6, 1/6];
                        case 3
                            GInt.pnt = [1/3, 3/5, 1/5, 1/5;
                                        1/3, 1/5, 3/5, 1/5];
                            GInt.wgt = [-9/32, 25/96, 25/96, 25/96];
                        case 4
                            GInt.pnt = [0.4459484909, 0.4459484909, 0.1081030182, 0.0915762135, 0.0915762135, 0.8168475729;
                                        0.4459484909, 0.1081030182, 0.4459484909, 0.0915762135, 0.8168475729, 0.0915762135];
                            GInt.wgt = [0.1116907948390055, 0.1116907948390055, 0.1116907948390055, ...
                                        0.054975871827661, 0.054975871827661, 0.054975871827661];
                        otherwise
                            error("Not supported integration order.");
                    end
                    GInt.orgTfm = Tfm(GInt.domn).orgTfm.getFun;
                    switch domn
                        case "D2T"
                            GInt.Jac = Tfm(GInt.domn).JDet.getFun;
                        case "D3F"
                            GInt.Jac = Tfm(GInt.domn).JNorm.getFun;
                    end
                case {"D2L", "D3L"}
                    assert(isequal(MshEnt(domn).dual.node.coord, sym([0, 1])));
                    switch ord
                        case {0, 1}
                            GInt.pnt = 1/2;
                            GInt.wgt = 1;
                        case {2, 3}
                            GInt.pnt = [(1 - sqrt(3) / 3) / 2, (1 + sqrt(3) / 3) / 2];
                            GInt.wgt = [1/2, 1/2];
                        case {4, 5}
                            GInt.pnt = [(1 - sqrt(15) / 5) / 2, 1/2, (1 + sqrt(15) / 5) / 2];
                            GInt.wgt = [5/18, 4/9, 5/18];
                        case {6, 7}
                            a = sqrt((3 + 2 * sqrt(6/5)) / 7);
                            b = sqrt((3 - 2 * sqrt(6/5)) / 7);
                            GInt.pnt = [(1 - a) / 2, (1 - b) / 2, (1 + b) / 2, (1 + a) / 2];
                            GInt.wgt = [(18 - sqrt(30)) / 72, (18 + sqrt(30)) / 72, (18 + sqrt(30)) / 72, (18 - sqrt(30)) / 72];
                        otherwise
                            error("Not supported integration order.");
                    end
                    GInt.orgTfm = Tfm(GInt.domn).orgTfm.getFun;
                    GInt.Jac = Tfm(GInt.domn).JNorm.getFun;
                case "D3T"
                    assert(isequal(MshEnt(domn).dual.node.coord, sym([0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1])));
                    switch ord
                        case {0, 1}
                            GInt.pnt = [1/4; 1/4; 1/4];
                            GInt.wgt = 1/6;
                        case 2
                            a = (5 - sqrt(5)) / 20;
                            b = (5 + 3 * sqrt(5)) / 20;
                            GInt.pnt = [a, b, a, a;
                                        a, a, b, a;
                                        a, a, a, b];
                            GInt.wgt = [1/24, 1/24, 1/24, 1/24];
                        case 3
                            GInt.pnt = [1/4, 1/6, 1/2, 1/6, 1/6;
                                        1/4, 1/6, 1/6, 1/2, 1/6;
                                        1/4, 1/6, 1/6, 1/6, 1/2];
                            GInt.wgt = [-2/15, 3/40, 3/40, 3/40, 3/40];
                        case {4, 5, 6, 7}
                            % Collapsed (Duffy) product of n-point Gauss-Legendre rules, algebraic precision 2n - 3:
                            % l = u, m = (1 - u) v, n = (1 - u) (1 - v) w, Jacobian (1 - u)^2 (1 - v).
                            [pnt1D, wgt1D] = gaussLeg(ceil((ord + 3) / 2));
                            [u, v, w] = ndgrid(pnt1D, pnt1D, pnt1D);
                            [wu, wv, ww] = ndgrid(wgt1D, wgt1D, wgt1D);
                            u = u(:)'; v = v(:)'; w = w(:)';
                            GInt.pnt = [u; (1 - u) .* v; (1 - u) .* (1 - v) .* w];
                            GInt.wgt = wu(:)' .* wv(:)' .* ww(:)' .* (1 - u).^2 .* (1 - v);
                        otherwise
                            error("Not supported integration order.");
                    end
                    GInt.orgTfm = Tfm(GInt.domn).orgTfm.getFun;
                    GInt.Jac = Tfm(GInt.domn).JDet.getFun;
            end
        end
        % Public functions.
        function val = eval(gInt, fun, parm)
            % GInt.eval: evaluate integral of function.
            arguments
                gInt GInt;
                fun;
                parm (:, :) = []; % Parameter of integrated domain.
                % If empty, function is integrated over reference domain.
            end
            if ~isempty(parm)
                assert(isequal(size(parm), gInt.sParm));
                % Jacobian of affine transformation is constant: evaluate at origin of reference domain.
                val = sum(fun(gInt.orgTfm(gInt.pnt, parm)) .* gInt.wgt) .* gInt.Jac(zeros(size(gInt.pnt, 1), 1), parm);
            else
                val = sum(fun(gInt.pnt) .* gInt.wgt);
            end
        end
        function val = evalVec(gInt, fun, parm, nEnt)
            % GInt.evalVec: evaluate integrals over a batch of mesh entities.
            % fun(X): integrands at points X (d x nPnt x nEnt), returns nFun x nPnt x nEnt (or broadcastable) values.
            % val(k, i): integral of k-th integrand over i-th mesh entity.
            arguments
                gInt GInt;
                fun;
                parm (:, :, :) = []; % Parameters of integrated domains: parm(:, :, i) for i-th mesh entity.
                % If empty, function is integrated over reference domain for `nEnt` mesh entities.
                nEnt = size(parm, 3);
            end
            [X, Jac] = gInt.pntVec(parm, nEnt);
            val = gInt.sumVec(fun(X), Jac);
        end
        function [X, Jac] = pntVec(gInt, parm, nEnt)
            % GInt.pntVec: Gauss points X (d x nPnt x nEnt) and Jacobians Jac (1 x 1 x nEnt) of a batch of mesh entities.
            % See `GInt.evalVec` for arguments.
            arguments
                gInt GInt;
                parm (:, :, :) = [];
                nEnt = size(parm, 3);
            end
            if ~isempty(parm)
                assert(isequal([size(parm, 1), size(parm, 2)], gInt.sParm));
                % Affine transformation from reference domain with vertices 0, e_1, ..., e_k (asserted in constructor).
                B = parm(:, 2:end, :) - parm(:, 1, :);
                X = pagemtimes(B, gInt.pnt) + parm(:, 1, :);
                Jac = jacVec(B);
            else
                X = repmat(gInt.pnt, 1, 1, nEnt);
                Jac = ones(1, 1, nEnt);
            end
        end
        function val = sumVec(gInt, F, Jac)
            % GInt.sumVec: integrals val(k, i) from values F (nFun x nPnt x nEnt, or broadcastable) at Gauss points X
            % and Jacobians Jac given by `GInt.pntVec`.
            nEnt = size(Jac, 3);
            F = F + zeros(1, size(gInt.pnt, 2), nEnt);
            val = reshape(sum(F .* gInt.wgt, 2) .* Jac, size(F, 1), nEnt);
        end
    end
end
% Local functions.
function Jac = jacVec(B)
    % jacVec: Jacobian of affine transformations x = B(:, :, i) * p + b (1 x 1 x n), consistent with `Tfm`.
    % Element (B square): signed determinant `JDet`; line or face: square root of Gram determinant `JNorm`.
    [d, k, ~] = size(B);
    if d == k
        switch d
            case 2
                Jac = B(1, 1, :) .* B(2, 2, :) - B(1, 2, :) .* B(2, 1, :);
            case 3
                Jac = B(1, 1, :) .* (B(2, 2, :) .* B(3, 3, :) - B(2, 3, :) .* B(3, 2, :)) ...
                    - B(1, 2, :) .* (B(2, 1, :) .* B(3, 3, :) - B(2, 3, :) .* B(3, 1, :)) ...
                    + B(1, 3, :) .* (B(2, 1, :) .* B(3, 2, :) - B(2, 2, :) .* B(3, 1, :));
        end
    else
        G = pagemtimes(B, "transpose", B, "none");
        switch k
            case 1
                Jac = sqrt(G);
            case 2
                Jac = sqrt(G(1, 1, :) .* G(2, 2, :) - G(1, 2, :) .* G(2, 1, :));
        end
    end
end
function [pnt, wgt] = gaussLeg(n)
    % gaussLeg: n-point Gauss-Legendre rule on [0, 1] by Golub-Welsch algorithm (by row).
    k = 1:n - 1;
    beta = k ./ sqrt(4 * k.^2 - 1);
    [V, D] = eig(diag(beta, 1) + diag(beta, -1));
    [pnt, idx] = sort(diag(D));
    pnt = (pnt' + 1) / 2;
    wgt = V(1, idx).^2;
end
