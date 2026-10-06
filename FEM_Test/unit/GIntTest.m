classdef GIntTest < matlab.unittest.TestCase
    % GIntTest: unit tests of GInt.

    % Algebraic precision of each rule:
    % domn | ord           | precision
    % ---- | ------------- | ---------
    % D2T  | 0, 1, 2, 3, 4 | 1, 1, 2, 3, 4
    % D2L  | 0, 1, ..., 7  | 1, 1, 3, 3, 5, 5, 7, 7
    % D3T  | 0, 1, ..., 7  | 1, 1, 2, 3, 5, 5, 7, 7
    % D3F  | same as D2T
    % D3L  | same as D2L

    properties
        P2T = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1]; % Vertices of a counter-clockwise triangle (by column).
        P2L = [0.2, 1.1; -0.3, 0.9]; % Vertices of a line (by column).
        P3T = [0.1, 1.2, 0.3, 0.2; -0.1, 0.2, 1.1, 0.3; 0.05, 0.1, 0.2, 1.3]; % Vertices of a positively oriented tetrahedron.
        P3F = [0.1, 1.2, 0.3; -0.1, 0.2, 1.1; 0.05, 0.1, 0.9]; % Vertices of a triangle in 3D.
    end
    methods (Test, TestTags = {'Fast'})
        function D2TRef(tc)
            % Monomial l^a m^b on reference triangle: a! b! / (a + b + 2)!.
            for ord = 0:4
                gInt = GInt("D2T", ord);
                tc.verifyEqual(sum(gInt.wgt), 1/2, "AbsTol", 1e-10);
                prec = max(ord, 1);
                for deg = 0:prec + 1
                    err = zeros(1, deg + 1);
                    for a = 0:deg
                        b = deg - a;
                        val = gInt.eval(@(p) p(1, :).^a .* p(2, :).^b);
                        err(a + 1) = abs(val - factorial(a) * factorial(b) / factorial(deg + 2));
                    end
                    if deg <= prec
                        tc.verifyLessThan(max(err), 1e-10, sprintf("GInt(D2T, %d) is not exact for degree %d.", ord, deg));
                    else
                        tc.verifyGreaterThan(max(err), 1e-6, sprintf("GInt(D2T, %d) is exact for degree %d.", ord, deg));
                    end
                end
            end
        end
        function D2TOrg(tc)
            % Compare with symbolic integration on original triangle.
            funs = ["x + 1", "x^2 + x*y + 1", "x^3 + x^2*y + y + 1", "x^4 + x*y^3 + x^2 + 1"];
            for ord = 1:4
                fcn = Fcn("D2", funs(ord));
                intVal = fcn.int("D2T").getFun;
                tc.verifyEqual(GInt("D2T", ord).eval(fcn.getFun, tc.P2T), intVal([0; 0], tc.P2T), "AbsTol", 1e-10);
            end
        end
        function parmSize(tc)
            % Parameter of integrated domain must match the domain.
            tc.verifyError(@() GInt("D2T", 1).eval(@(p) p(1, :), [0, 1; 0, 1]), "MATLAB:assertion:failed");
            tc.verifyError(@() GInt("D2L", 1).eval(@(p) p(1, :), [0, 1, 0; 0, 0, 1]), "MATLAB:assertion:failed");
            tc.verifyError(@() GInt("D3T", 1).eval(@(p) p(1, :), [0, 1, 0; 0, 0, 1; 0, 0, 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() GInt("D3F", 1).eval(@(p) p(1, :), [0, 1, 0; 0, 0, 1]), "MATLAB:assertion:failed");
        end
        function D3TRef(tc)
            % Monomial l^a m^b n^c on reference tetrahedron: a! b! c! / (a + b + c + 3)!.
            precs = [1, 1, 2, 3, 5, 5, 7, 7];
            for ord = 0:7
                gInt = GInt("D3T", ord);
                tc.verifyEqual(sum(gInt.wgt), 1/6, "AbsTol", 1e-14);
                % Gauss points are inside the tetrahedron.
                tc.verifyTrue(all(gInt.pnt >= 0, "all") && all(sum(gInt.pnt, 1) <= 1));
                prec = precs(ord + 1);
                for deg = 0:prec + 1
                    err = 0;
                    for a = 0:deg
                        for b = 0:deg - a
                            c = deg - a - b;
                            val = gInt.eval(@(p) p(1, :).^a .* p(2, :).^b .* p(3, :).^c);
                            err = max(err, abs(val - factorial(a) * factorial(b) * factorial(c) / factorial(deg + 3)));
                        end
                    end
                    if deg <= prec
                        tc.verifyLessThan(err, 1e-12, sprintf("GInt(D3T, %d) is not exact for degree %d.", ord, deg));
                    else
                        tc.verifyGreaterThan(err, 1e-9, sprintf("GInt(D3T, %d) is exact for degree %d.", ord, deg));
                    end
                end
            end
        end
        function D3TOrg(tc)
            P = tc.P3T;
            vol = det(P(:, 2:4) - P(:, 1)) / 6;
            x = P(1, :); y = P(2, :);
            % int_K lambda_i lambda_j = vol (1 + delta_ij) / 20.
            tc.verifyEqual(GInt("D3T", 1).eval(@(p) ones(1, size(p, 2)), P), vol, "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3T", 1).eval(@(p) p(1, :), P), vol * mean(x), "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3T", 2).eval(@(p) p(1, :).^2, P), vol / 20 * (sum(x.^2) + sum(x)^2), "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3T", 2).eval(@(p) p(1, :) .* p(2, :), P), vol / 20 * (sum(x .* y) + sum(x) * sum(y)), "AbsTol", 1e-13);
            % Rules of different orders agree on polynomial of degree 5.
            fun = @(p) p(1, :).^3 .* p(2, :) .* p(3, :) + p(2, :).^4 - p(3, :).^2 + 1;
            tc.verifyEqual(GInt("D3T", 4).eval(fun, P), GInt("D3T", 7).eval(fun, P), "AbsTol", 1e-13);
        end
        function D3F(tc)
            P = tc.P3F;
            area = norm(cross(P(:, 2) - P(:, 1), P(:, 3) - P(:, 1))) / 2;
            x = P(1, :);
            % int_F lambda_i lambda_j = area (1 + delta_ij) / 12.
            tc.verifyEqual(GInt("D3F", 1).eval(@(p) ones(1, size(p, 2)), P), area, "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3F", 1).eval(@(p) p(1, :) + p(3, :), P), area * (mean(x) + mean(P(3, :))), "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3F", 2).eval(@(p) p(1, :).^2, P), area / 12 * (sum(x.^2) + sum(x)^2), "AbsTol", 1e-13);
            % In plane z = 0, same as D2T.
            P2 = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1];
            fun3 = @(p) p(1, :).^3 + p(1, :) .* p(2, :).^2 + 1;
            fun2 = @(p) p(1, :).^3 + p(1, :) .* p(2, :).^2 + 1;
            tc.verifyEqual(GInt("D3F", 3).eval(fun3, [P2; 0, 0, 0]), GInt("D2T", 3).eval(fun2, P2), "AbsTol", 1e-13);
        end
        function D3L(tc)
            P = [0, 1; 0, 2; 0, 2];
            tc.verifyEqual(GInt("D3L", 1).eval(@(p) ones(1, size(p, 2)), P), 3, "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3L", 1).eval(@(p) p(1, :), P), 1.5, "AbsTol", 1e-13);
            tc.verifyEqual(GInt("D3L", 3).eval(@(p) p(1, :).^2 + p(2, :) .* p(3, :), P), 3 * (1/3 + 4/3), "AbsTol", 1e-13);
        end
        function evalVec(tc)
            % Batch integration agrees with integration on each domain.
            rng(3);
            nEnt = 5;
            domns = ["D2T", "D2L", "D3T", "D3F", "D3L"];
            sParms = {[2, 3], [2, 2], [3, 4], [3, 3], [3, 2]};
            ords = [4, 7, 7, 4, 5];
            for iDomn = 1:length(domns)
                gInt = GInt(domns(iDomn), ords(iDomn));
                parm = rand([sParms{iDomn}, nEnt]);
                fun = @(x) x(1, :, :).^2 + sin(x(end, :, :)) + 1;
                val = gInt.evalVec(fun, parm);
                tc.verifyEqual(size(val), [1, nEnt]);
                for iEnt = 1:nEnt
                    tc.verifyEqual(val(iEnt), gInt.eval(fun, parm(:, :, iEnt)), "AbsTol", 1e-13, domns(iDomn));
                end
                % Reference domain.
                val = gInt.evalVec(@(x) x(1, :, :) + 1, [], nEnt);
                tc.verifyEqual(val, repmat(gInt.eval(@(x) x(1, :) + 1), 1, nEnt), "AbsTol", 1e-14, domns(iDomn));
                % Several integrands at once: one row per integrand.
                val = gInt.evalVec(@(x) [x(1, :, :); x(end, :, :).^2], parm);
                tc.verifyEqual(size(val), [2, nEnt]);
                tc.verifyEqual(val(2, :), arrayfun(@(i) gInt.eval(@(x) x(end, :).^2, parm(:, :, i)), 1:nEnt), "AbsTol", 1e-13);
                % Constant integrand: measure of domain.
                val = gInt.evalVec(@(x) 1, parm);
                tc.verifyEqual(val, arrayfun(@(i) gInt.eval(@(x) ones(1, size(x, 2)), parm(:, :, i)), 1:nEnt), "AbsTol", 1e-14);
            end
        end
        function D2LRef(tc)
            % Monomial s^k on reference line: 1 / (k + 1).
            for ord = 0:7
                gInt = GInt("D2L", ord);
                tc.verifyEqual(sum(gInt.wgt), 1, "AbsTol", 1e-14);
                prec = max(ord + mod(ord + 1, 2), 1);
                for deg = 0:prec + 1
                    err = abs(gInt.eval(@(p) p.^deg) - 1 / (deg + 1));
                    if deg <= prec
                        tc.verifyLessThan(err, 1e-13, sprintf("GInt(D2L, %d) is not exact for degree %d.", ord, deg));
                    else
                        tc.verifyGreaterThan(err, 1e-6, sprintf("GInt(D2L, %d) is exact for degree %d.", ord, deg));
                    end
                end
            end
        end
        function D2LOrg(tc)
            % Compare with symbolic integration on original line.
            funs = ["x + y + 1", "x^3 + x*y^2 + 1", "x^5 + y^4 + x", "x^7 + x^3*y^3 + 1"];
            ords = [1, 3, 5, 7];
            for iFun = 1:length(funs)
                fcn = Fcn("D2", funs(iFun));
                intVal = fcn.int("D2L").getFun;
                tc.verifyEqual(GInt("D2L", ords(iFun)).eval(fcn.getFun, tc.P2L), intVal([0; 0], tc.P2L), "AbsTol", 1e-12);
            end
        end
    end
end
