classdef FcnTest < matlab.unittest.TestCase
    % FcnTest: unit tests of Fcn.
    properties
        P2T = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1]; % Vertices of a counter-clockwise triangle (by column).
        P3T = [0.1, 1.2, 0.3, 0.2; -0.1, 0.2, 1.1, 0.3; 0.05, 0.1, 0.2, 1.3]; % Vertices of a positively oriented tetrahedron.
        P3F = [0.1, 1.2, 0.3; -0.1, 0.2, 1.1; 0.05, 0.1, 0.9]; % Vertices of a triangle in 3D.
    end
    methods (Test, TestTags = {'Fast'})
        function constructor(tc)
            fcn = Fcn("D2T", "[c1*x; c2*y]", "[c1; c2]");
            tc.verifyEqual(fcn.domn, "D2T");
            tc.verifyTrue(symEq(fcn.fun, str2sym("[c1*x; c2*y]")));
            tc.verifyTrue(symEq(fcn.coef, str2sym("[c1; c2]")));
            tc.verifyTrue(symEq(fcn.var, str2sym("[x; y]")));
            tc.verifyTrue(symEq(fcn.parm, str2sym("[x1, x2, x3; y1, y2, y3]")));
            tc.verifyEqual(fcn.nVar, 2);
            tc.verifyEqual(fcn.nFun, 2);
            tc.verifyEqual(fcn.nCoef, 2);
            tc.verifyEqual(fcn.sParm, [2, 3]);
            tc.verifyEqual(Fcn("D2T", "[x, y; y, x]").nFun, 4);
            tc.verifyEqual(Fcn.cst(3).domn, "VOID");
            tc.verifyTrue(symEq(Fcn.cst(3).fun, sym(3)));
        end
        function domnMerge(tc)
            tc.verifyEqual(getDomn([Fcn("D2", "1"), Fcn("D2T", "1")]), "D2T");
            tc.verifyEqual(getDomn([Fcn("D2L", "1"), Fcn("D2T", "1")]), "D2T");
            tc.verifyEqual(getDomn([Fcn.cst(1), Fcn("D2R1", "s")]), "D2R1");
            tc.verifyError(@() getDomn([Fcn("D2T", "1"), Fcn("D2TR", "1")]), "MATLAB:assertion:failed");
        end
        function funHandle(tc)
            fun = Fcn("D2", "x + 2*y").getFun;
            tc.verifyEqual(fun([1, 2; 3, 4]), [7, 10], "AbsTol", 1e-14);
            fun = Fcn("D2", "1").getFun;
            tc.verifyEqual(fun([1, 2, 3; 1, 2, 3]), [1, 1, 1]);
            fun = Fcn("D2T", "c1*x + c2*y + x1 + x2", "[c1; c2]").getFun;
            tc.verifyEqual(fun([1; 1], [1, 2, 3; 4, 5, 6], [10; 20]), 33, "AbsTol", 1e-14);
            fun = Fcn("D2", "c1*x + x1").getFun("parm", {"[x1, x2; y1, y2]"}, "coef", {"c1"});
            tc.verifyEqual(fun([2; 0], [1, 0; 0, 0], 3), 7, "AbsTol", 1e-14);
        end
        function funHandleVec(tc)
            % Batch evaluation agrees with point-wise evaluation.
            rng(2);
            nPnt = 3; nEnt = 4;
            fcn = Fcn("D3T", "x*y1 + z^2*x4 - y*z3 + 1") .* Fcn("D3T", "x2 - x");
            fun = fcn.getFun;
            vecFun = fcn.getFun("vec", true);
            X = rand(3, nPnt, nEnt);
            P = rand(3, 4, nEnt);
            val = vecFun(X, reshape(P, [], 1, nEnt));
            tc.verifyEqual(size(val), [1, nPnt, nEnt]);
            for iEnt = 1:nEnt
                tc.verifyEqual(val(1, :, iEnt), fun(X(:, :, iEnt), P(:, :, iEnt)), "AbsTol", 1e-14);
            end
            % Manually set parameter and coefficient.
            fcn = Fcn("D2", "x*p1_1 + y*p2_2 - p1_2*c1 + c2^2");
            vecFun = fcn.getFun("parm", {sym("p", [2, 2])}, "coef", {sym("c", [2, 1])}, "vec", true);
            X = rand(2, nPnt, nEnt); P = rand(2, 2, nEnt); C = rand(2, 1, nEnt);
            val = vecFun(X, reshape(P, [], 1, nEnt), C);
            ref = X(1, :, :) .* P(1, 1, :) + X(2, :, :) .* P(2, 2, :) - P(1, 2, :) .* C(1, 1, :) + C(2, 1, :).^2;
            tc.verifyEqual(val, ref, "AbsTol", 1e-14);
            % Function independent of points: output is broadcastable.
            vecFun = Fcn("D2T", "x1 + 2").getFun("vec", true);
            P = rand(2, 3, nEnt);
            val = vecFun(rand(2, nPnt, nEnt), reshape(P, [], 1, nEnt)) + zeros(1, nPnt, nEnt);
            tc.verifyEqual(val, repmat(P(1, 1, :) + 2, 1, nPnt, 1), "AbsTol", 1e-14);
            % Array of functions with the same arguments: one output row per function (constant included).
            fcns = [Fcn("D2T", "x*y1 - y"), Fcn("D2T", "3"), Fcn("D2T", "x^2 + x2")];
            vecFun = fcns.getFun("vec", true);
            X = rand(2, nPnt, nEnt);
            val = vecFun(X, reshape(P, [], 1, nEnt));
            tc.verifyEqual(size(val), [3, nPnt, nEnt]);
            for iFcn = 1:3
                fun = fcns(iFcn).getFun;
                for iEnt = 1:nEnt
                    tc.verifyEqual(val(iFcn, :, iEnt), fun(X(:, :, iEnt), P(:, :, iEnt)), "AbsTol", 1e-14);
                end
            end
            % Functions with different arguments cannot be combined.
            tc.verifyError(@() getFun([Fcn("D2T", "x"), Fcn("D2", "x")], "vec", true), "MATLAB:assertion:failed");
            % Vector-valued function is not supported.
            tc.verifyError(@() Fcn("D3", "[x; y]").getFun("vec", true), "MATLAB:assertion:failed");
        end
        function operator(tc)
            fcn = Fcn("D2", "x") + Fcn("D2T", "y");
            tc.verifyEqual(fcn.domn, "D2T");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x + y")));
            fcn = Fcn("D2T", "c1*x", "c1") + Fcn("D2T", "c2*y", "c2");
            tc.verifyTrue(symEq(fcn.coef, str2sym("[c1; c2]")));
            fcn = Fcn("D2", "x") - Fcn("D2", "y");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x - y")));
            fcn = -Fcn("D2", "x");
            tc.verifyTrue(symEq(fcn.fun, str2sym("-x")));
            fcn = Fcn("D2", "x") * 2;
            tc.verifyEqual(fcn.domn, "D2");
            tc.verifyTrue(symEq(fcn.fun, str2sym("2*x")));
            fcn = Fcn("D2", "x") \ 3;
            tc.verifyTrue(symEq(fcn.fun, str2sym("3/x")));
            fcn = Fcn("D2", "[x; y]") .* Fcn("D2", "[y; x]");
            tc.verifyTrue(symEq(fcn.fun, str2sym("[x*y; x*y]")));
            fcn = Fcn("D2", "x") .^ 2;
            tc.verifyTrue(symEq(fcn.fun, str2sym("x^2")));
            fcn = Fcn("D2", "[x, y]").';
            tc.verifyTrue(symEq(fcn.fun, str2sym("[x; y]")));
            fcn = dot(Fcn("D2", "[x; y]"), Fcn("D2", "[1; 2]"));
            tc.verifyTrue(symEq(fcn.fun, str2sym("x + 2*y")));
            fcn = sum(Fcn("D2", "[x, y; 1, 2]"));
            tc.verifyTrue(symEq(fcn.fun, str2sym("x + y + 3")));
            fcn = sum(Fcn("D2", "[x, y; 1, 2]"), 2);
            tc.verifyTrue(symEq(fcn.fun, str2sym("[x + y; 3]")));
        end
        function evalComb(tc)
            tc.verifyTrue(symEq(Fcn("D2", "x + y").eval([1; 2]), sym(3)));
            fcn = comb([Fcn("D2", "x"), Fcn("D2", "y")], [2; 3]);
            tc.verifyEqual(fcn.domn, "D2");
            tc.verifyTrue(symEq(fcn.fun, str2sym("2*x + 3*y")));
            tc.verifyEmpty(fcn.coef);
            fcn = comb([Fcn("D2", "x"), Fcn("D2T", "y")], "[c1; c2]");
            tc.verifyEqual(fcn.domn, "D2T");
            tc.verifyTrue(symEq(fcn.fun, str2sym("c1*x + c2*y")));
            tc.verifyTrue(symEq(fcn.coef, str2sym("[c1; c2]")));
        end
        function difOrd(tc)
            fcn = Fcn("D2", "x^2*y").dif([1; 0]);
            tc.verifyTrue(symEq(fcn.fun, str2sym("2*x*y")));
            fcn = Fcn("D2", "x^2*y").dif([1; 1]);
            tc.verifyTrue(symEq(fcn.fun, str2sym("2*x")));
            fcn = Fcn("D2", "x^2*y").dif(cat(3, [1; 0], [0; 1]));
            tc.verifyTrue(symEq(fcn.fun, str2sym("[2*x*y; x^2]")));
            fcn = Fcn("D2", "[x^2; x*y]").dif([1, 0; 0, 1]);
            tc.verifyTrue(symEq(fcn.fun, str2sym("[2*x; x]")));
            fcn = Fcn("D2", "[x^2; x*y]").dif([1, nan; 0, nan]);
            tc.verifyTrue(symEq(fcn.fun, str2sym("[2*x; 0]")));
            fcn = Fcn("D2", "[x^2; x*y]").dif(cat(3, [1, 1; 0, 0], [0, 0; 1, 1]));
            tc.verifyTrue(symEq(fcn.fun, str2sym("[2*x, 0; y, x]")));
        end
        function tfmDomn(tc)
            fcn = Fcn("D2", "x").tfm("D2TR");
            tc.verifyEqual(fcn.domn, "D2TR");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x1 + l*(x2 - x1) + m*(x3 - x1)")));
            fcn = Fcn("D2R", "l").tfm("D2T");
            tc.verifyEqual(fcn.domn, "D2T");
            tc.verifyTrue(symEq(fcn.eval(str2sym("[x2; y2]")), sym(1)));
            tc.verifyTrue(symEq(fcn.eval(str2sym("[x3; y3]")), sym(0)));
            fcn = Fcn("D2", "x").tfm("D2LR");
            tc.verifyEqual(fcn.domn, "D2LR");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x1 + s*(x2 - x1)")));
            % Round trip.
            fcn = Fcn("D2", "x^2 + y").tfm("D2TR").tfm("D2T");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x^2 + y")));
            % Array operation.
            fcns = tfm([Fcn("D2", "x"), Fcn("D2", "y")], "D2TR");
            tc.verifyTrue(symEq(fcns(2).fun, str2sym("y1 + l*(y2 - y1) + m*(y3 - y1)")));
        end
        function intDomn(tc)
            P = tc.P2T;
            area = ((P(1, 2) - P(1, 1)) * (P(2, 3) - P(2, 1)) - (P(1, 3) - P(1, 1)) * (P(2, 2) - P(2, 1))) / 2;
            % Reference domain.
            tc.verifyTrue(symEq(Fcn("D2R", "1").int("D2TR").fun, sym(1/2)));
            tc.verifyTrue(symEq(Fcn("D2R", "l").int("D2TR").fun, sym(1/6)));
            tc.verifyTrue(symEq(Fcn("D2R", "l*m").int("D2TR").fun, sym(1/24)));
            tc.verifyTrue(symEq(Fcn("D2R1", "s^2").int("D2LR").fun, sym(1/3)));
            % Original triangle.
            intVal = Fcn("D2", "1").int("D2T");
            tc.verifyEqual(intVal.domn, "D2T");
            tc.verifyTrue(symEq(intVal.fun, str2sym("((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))/2")));
            tc.verifyEqual(numSym(Fcn("D2", "x").int("D2T").fun, "D2T", P), area * mean(P(1, :)), "AbsTol", 1e-12);
            % Original line.
            tc.verifyEqual(numSym(Fcn("D2", "x").int("D2L").fun, "D2L", [0, 3; 0, 4]), 7.5, "AbsTol", 1e-12);
            % Edge of original triangle.
            tc.verifyEqual(numSym(Fcn("D2T", "1").int("D2L", 2).fun, "D2T", P), norm(P(:, 3) - P(:, 2)), "AbsTol", 1e-12);
            tc.verifyEqual(numSym(Fcn("D2T", "x").int("D2L", 3).fun, "D2T", P), ...
                norm(P(:, 1) - P(:, 3)) * (P(1, 1) + P(1, 3)) / 2, "AbsTol", 1e-12);
            % Edge of reference triangle.
            tc.verifyTrue(symEq(Fcn("D2R", "l").int("D2LR", 1).fun, sym(1/2)));
            tc.verifyTrue(symEq(Fcn("D2R", "l").int("D2LR", 2).fun, sqrt(sym(2)) / 2));
            tc.verifyTrue(symEq(Fcn("D2R", "l").int("D2LR", 3).fun, sym(0)));
        end
        function tfmDomn3D(tc)
            P = tc.P3T;
            fcn = Fcn("D3", "x").tfm("D3TR");
            tc.verifyEqual(fcn.domn, "D3TR");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x1 + l*(x2 - x1) + m*(x3 - x1) + n*(x4 - x1)")));
            fcn = Fcn("D3R", "l").tfm("D3T");
            tc.verifyEqual(fcn.domn, "D3T");
            for iNode = 1:4
                tc.verifyEqual(numSym(fcn.eval(P(:, iNode)), "D3T", P), double(iNode == 2), "AbsTol", 1e-12);
            end
            fcn = Fcn("D3", "x").tfm("D3FR");
            tc.verifyEqual(fcn.domn, "D3FR");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x1 + s*(x2 - x1) + t*(x3 - x1)")));
            fcn = Fcn("D3", "z").tfm("D3LR");
            tc.verifyEqual(fcn.domn, "D3LR");
            tc.verifyTrue(symEq(fcn.fun, str2sym("z1 + s*(z2 - z1)")));
            % Round trip.
            fcn = Fcn("D3", "x^2 + y*z").tfm("D3TR").tfm("D3T");
            pnt = [0.3; 0.2; 0.4];
            tc.verifyEqual(numSym(fcn.eval(pnt), "D3T", P), pnt(1)^2 + pnt(2) * pnt(3), "AbsTol", 1e-12);
            % 2D and 3D domains are not compatible.
            tc.verifyError(@() Fcn("D2", "x").tfm("D3TR"), "MATLAB:assertion:failed");
        end
        function intDomn3D(tc)
            P = tc.P3T;
            vol = det(P(:, 2:4) - P(:, 1)) / 6;
            % Reference domain.
            tc.verifyTrue(symEq(Fcn("D3R", "1").int("D3TR").fun, sym(1/6)));
            tc.verifyTrue(symEq(Fcn("D3R", "l").int("D3TR").fun, sym(1/24)));
            tc.verifyTrue(symEq(Fcn("D3R", "l*m*n").int("D3TR").fun, sym(1/720)));
            tc.verifyTrue(symEq(Fcn("D3R2", "s*t").int("D3FR").fun, sym(1/24)));
            tc.verifyTrue(symEq(Fcn("D3R1", "s^2").int("D3LR").fun, sym(1/3)));
            % Original tetrahedron, face and line.
            intVal = Fcn("D3", "x").int("D3T");
            tc.verifyEqual(intVal.domn, "D3T");
            tc.verifyEqual(numSym(intVal.fun, "D3T", P), vol * mean(P(1, :)), "AbsTol", 1e-12);
            area = norm(cross(tc.P3F(:, 2) - tc.P3F(:, 1), tc.P3F(:, 3) - tc.P3F(:, 1))) / 2;
            tc.verifyEqual(numSym(Fcn("D3", "1").int("D3F").fun, "D3F", tc.P3F), area, "AbsTol", 1e-12);
            tc.verifyEqual(numSym(Fcn("D3", "x").int("D3L").fun, "D3L", [0, 1; 0, 2; 0, 2]), 1.5, "AbsTol", 1e-12);
            % Faces of original tetrahedron: areas, and divergence theorem for [x*y; 0; z] (divergence y + 1).
            mshEnt = MshEnt("D3T");
            UNV = mshEnt.UNV;
            flux = 0;
            for iFace = 1:4
                FcNd = P(:, mshEnt.face.node(:, iFace));
                tc.verifyEqual(numSym(Fcn("D3T", "1").int("D3F", iFace).fun, "D3T", P), ...
                    norm(cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1))) / 2, "AbsTol", 1e-12);
                flux = flux + numSym(int(dot(Fcn("D3T", "[x*y; 0; z]"), UNV(iFace)), "D3F", iFace).fun, "D3T", P);
            end
            tc.verifyEqual(flux, vol * (mean(P(2, :)) + 1), "AbsTol", 1e-12);
            % Edge of original tetrahedron.
            tc.verifyEqual(numSym(Fcn("D3T", "1").int("D3L", 6).fun, "D3T", P), norm(P(:, 4) - P(:, 3)), "AbsTol", 1e-12);
            % Faces and edges of reference tetrahedron.
            tc.verifyTrue(symEq(Fcn("D3R", "1").int("D3FR", 1).fun, sqrt(sym(3)) / 2));
            tc.verifyTrue(symEq(Fcn("D3R", "1").int("D3FR", 2).fun, sym(1/2)));
            tc.verifyTrue(symEq(Fcn("D3R", "l").int("D3FR", 4).fun, sym(1/6)));
            tc.verifyTrue(symEq(Fcn("D3R", "l").int("D3FR", 2).fun, sym(0)));
            tc.verifyTrue(symEq(Fcn("D3R", "l").int("D3LR", 1).fun, sym(1/2)));
            tc.verifyTrue(symEq(Fcn("D3R", "1").int("D3LR", 6).fun, sqrt(sym(2))));
        end
        function parmCoef(tc)
            tc.verifyEqual(Fcn("D2", "x").clrParm.domn, "D2");
            tc.verifyEqual(Fcn("D2T", "x").clrParm.domn, "D2");
            tc.verifyEqual(Fcn("D2L", "x").clrParm.domn, "D2");
            tc.verifyEqual(Fcn("D2TR", "l").clrParm.domn, "D2R");
            tc.verifyEqual(Fcn("D2LR", "s").clrParm.domn, "D2R1");
            tc.verifyEqual(Fcn("D3T", "x").clrParm.domn, "D3");
            tc.verifyEqual(Fcn("D3FR", "s").clrParm.domn, "D3R2");
            tc.verifyEqual(Fcn("D3LR", "s").clrParm.domn, "D3R1");
            fcn = Fcn("D2T", "x + x2").subParm(sym([0, 1, 0; 0, 0, 1]));
            tc.verifyEqual(fcn.domn, "D2");
            tc.verifyTrue(symEq(fcn.fun, str2sym("x + 1")));
            fcn = Fcn("D2T", "c1*x", "c1").subCoef(sym(2));
            tc.verifyTrue(symEq(fcn.fun, str2sym("2*x")));
            tc.verifyEmpty(fcn.coef);
            fcn = Fcn("D2T", "c1*x", "c1").clrCoef;
            tc.verifyEmpty(fcn.coef);
        end
    end
end
