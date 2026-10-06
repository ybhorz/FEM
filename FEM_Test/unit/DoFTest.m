classdef DoFTest < matlab.unittest.TestCase
    % DoFTest: unit tests of DoF, NdDoF and MoDoF.
    properties
        P2T = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1]; % Vertices of a counter-clockwise triangle (by column).
        P3T = [0.1, 1.2, 0.3, 0.2; -0.1, 0.2, 1.1, 0.3; 0.05, 0.1, 0.2, 1.3]; % Vertices of a positively oriented tetrahedron.
    end
    methods (Test, TestTags = {'Fast'})
        function ndDoFVert(tc)
            D2TElem = MshEnt("D2T").msh;
            % Scalar function at vertex.
            val = eval(NdDoF("D2", D2TElem, 0, [], [0; 0]), Fcn("D2", "x + y"));
            tc.verifyTrue(symEq(val, str2sym("[x1 + y1; x2 + y2; x3 + y3]")));
            val = eval(NdDoF("D2", D2TElem, 0, [], [0; 0], "EntIdx", 2), Fcn("D2", "x + y"));
            tc.verifyTrue(symEq(val, str2sym("x2 + y2")));
            % Derivative at vertex.
            val = eval(NdDoF("D2", D2TElem, 0, [], [1; 0]), Fcn("D2", "x^2 + x*y + y^2"));
            tc.verifyTrue(symEq(val, str2sym("[2*x1 + y1; 2*x2 + y2; 2*x3 + y3]")));
            val = eval(NdDoF("D2", D2TElem, 0, [], [0; 1]), Fcn("D2", "x^2 + x*y + y^2"));
            tc.verifyTrue(symEq(val, str2sym("[x1 + 2*y1; x2 + 2*y2; x3 + 2*y3]")));
            % Vector function at vertex.
            val = eval(NdDoF("D2", D2TElem, 0, [], [0, nan; 0, nan]), Fcn("D2", "[x + y; x - y]"));
            tc.verifyTrue(symEq(val, str2sym("[x1 + y1; x2 + y2; x3 + y3]")));
            val = eval(NdDoF("D2", D2TElem, 0, [], [nan, 0; nan, 0]), Fcn("D2", "[x + y; x - y]"));
            tc.verifyTrue(symEq(val, str2sym("[x1 - y1; x2 - y2; x3 - y3]")));
        end
        function ndDoFEdge(tc)
            D2TElem = MshEnt("D2T").msh;
            % Scalar function on edge.
            val = eval(NdDoF("D2", D2TElem, 1, 1/2, [0; 0]), Fcn("D2", "x"));
            tc.verifyTrue(symEq(val, str2sym("[(x1 + x2)/2; (x2 + x3)/2; (x3 + x1)/2]")));
            val = eval(NdDoF("D2", D2TElem, 1, [0, 1], [0; 0]), Fcn("D2", "x"));
            tc.verifyTrue(symEq(val, str2sym("[x1, x2; x2, x3; x3, x1]")));
            % Normal component of vector function on edge.
            P = tc.P2T;
            val = eval(NdDoF("D2", D2TElem, 1, 1/2, [0, 0; 0, 0], "coef", MshEnt("D2L").UNV), Fcn("D2", "[x; y]"));
            ref = zeros(3, 1);
            for iEdge = 1:3
                EgNd1 = P(:, D2TElem.edge.node(1, iEdge));
                EgNd2 = P(:, D2TElem.edge.node(2, iEdge));
                tan = EgNd2 - EgNd1;
                ref(iEdge) = dot([tan(2); -tan(1)] / norm(tan), (EgNd1 + EgNd2) / 2);
            end
            tc.verifyEqual(numSym(val, "D2T", P), ref, "AbsTol", 1e-12);
            % Trace function on edge.
            val = eval(NdDoF("D2R1", D2TElem, 1, [0, 1], 0), Fcn("D2", "x"));
            tc.verifyTrue(symEq(val, str2sym("[x1, x2; x2, x3; x3, x1]")));
            val = eval(NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1], 0), Fcn("D2R1", "s"));
            tc.verifyTrue(symEq(val, sym([0, 1])));
            % Numerical evaluation agrees with symbolic evaluation.
            val = eval(NdDoF("D2R1", MshEnt("D2LR").msh, 1, [0, 1/3, 1], 0), Fcn("D2R1", "s^2 + 1"), "valType", "num");
            tc.verifyEqual(val, [1, 10/9, 2], "AbsTol", 1e-14);
            msh = mshD2TS([0, 1, 0, 1], 2);
            ndDoF = NdDoF("D2R1", msh, 1, [0, 1/3, 1], 0);
            fcn = Fcn("D2", "x^2 + x*y");
            tc.verifyEqual(ndDoF.eval(fcn, "valType", "num"), double(ndDoF.eval(fcn)), "AbsTol", 1e-14);
        end
        function ndDoFElem(tc)
            D2TElem = MshEnt("D2T").msh;
            val = eval(NdDoF("D2", D2TElem, 2, [1/3; 1/3], [0; 0]), Fcn("D2", "x"));
            tc.verifyTrue(symEq(val, str2sym("(x1 + x2 + x3)/3")));
            val = eval(NdDoF("D2", D2TElem, 2, [0, 0; 1, 0; 0, 1]', [0; 0]), Fcn("D2", "x"));
            tc.verifyTrue(symEq(val, str2sym("[x1, x2, x3]")));
            val = eval(NdDoF("D2", D2TElem, 2, [1/2, 0; 1/2, 1/2; 0, 1/2]', [0; 0]), Fcn("D2", "x"));
            tc.verifyTrue(symEq(val, str2sym("[(x1 + x2)/2, (x2 + x3)/2, (x3 + x1)/2]")));
        end
        function moDoFEdge(tc)
            D2TElem = MshEnt("D2T").msh;
            P = tc.P2T;
            len = zeros(3, 1); dx = zeros(3, 1); tanDif = zeros(3, 1);
            for iEdge = 1:3
                tan = P(:, D2TElem.edge.node(2, iEdge)) - P(:, D2TElem.edge.node(1, iEdge));
                len(iEdge) = norm(tan);
                dx(iEdge) = tan(1);
                tanDif(iEdge) = tan(2) - tan(1);
            end
            gInt = GInt("D2L", 4);
            % Scalar function on edge.
            val = eval(MoDoF("D2", D2TElem, 1, Fcn("D2R1", "1"), [0; 0], "GInt", gInt), Fcn("D2", "1"));
            tc.verifyEqual(numSym(val, "D2T", P), len, "AbsTol", 1e-12);
            val = eval(MoDoF("D2", D2TElem, 1, Fcn("D2R1", "s - 1/2"), [0; 0], "GInt", gInt), Fcn("D2", "1"));
            tc.verifyTrue(symEq(val, sym([0; 0; 0])));
            % int_e x (s - 1/2) ds = |e| (x2 - x1) / 12.
            val = eval(MoDoF("D2", D2TElem, 1, Fcn("D2R1", "s - 1/2"), [0; 0], "GInt", gInt), Fcn("D2", "x"));
            tc.verifyEqual(numSym(val, "D2T", P), len .* dx / 12, "AbsTol", 1e-12);
            % Normal flux: int_e [1; 1] . n ds = t2 - t1.
            val = eval(MoDoF("D2", D2TElem, 1, Fcn("D2R1", "1"), [0, 0; 0, 0], "coef", MshEnt("D2L").UNV, "GInt", gInt), ...
                Fcn("D2", "[1; 1]"));
            tc.verifyEqual(numSym(val, "D2T", P), tanDif, "AbsTol", 1e-12);
            % Trace function.
            val = eval(MoDoF("D2R1", D2TElem, 1, Fcn("D2R1", "1"), 0, "GInt", gInt), Fcn("D2", "1"));
            tc.verifyTrue(symEq(val, sym([1; 1; 1])));
            val = eval(MoDoF("D2R1", D2TElem, 1, Fcn("D2R1", "s - 1/2"), 0, "GInt", gInt), Fcn("D2", "1"));
            tc.verifyTrue(symEq(val, sym([0; 0; 0])));
            val = eval(MoDoF("D2R1", MshEnt("D2LR").msh, 1, Fcn("D2R1", "1"), 0, "GInt", gInt), Fcn("D2R1", "1"));
            tc.verifyTrue(symEq(val, sym(1)));
            val = eval(MoDoF("D2R1", MshEnt("D2LR").msh, 1, Fcn("D2R1", "s - 1/2"), 0, "GInt", gInt), Fcn("D2R1", "s"));
            tc.verifyTrue(symEq(val, sym(1/12)));
        end
        function moDoFElem(tc)
            D2TElem = MshEnt("D2T").msh;
            val = eval(MoDoF("D2", D2TElem, 2, Fcn("D2R", "1"), [0; 0], "GInt", GInt("D2T", 2)), Fcn("D2", "1"));
            tc.verifyTrue(symEq(val, str2sym("((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))/2")));
        end
        function ndDoFVert3D(tc)
            D3TElem = MshEnt("D3T").msh;
            val = eval(NdDoF("D3", D3TElem, 0, [], [0; 0; 0]), Fcn("D3", "x + y + z"));
            tc.verifyTrue(symEq(val, str2sym("[x1 + y1 + z1; x2 + y2 + z2; x3 + y3 + z3; x4 + y4 + z4]")));
            % Derivative at vertex.
            val = eval(NdDoF("D3", D3TElem, 0, [], [0; 0; 1], "EntIdx", [2, 4]), Fcn("D3", "x*z + z^2"));
            tc.verifyTrue(symEq(val, str2sym("[x2 + 2*z2; x4 + 2*z4]")));
            % Vector function at vertex.
            val = eval(NdDoF("D3", D3TElem, 0, [], [nan, nan, 0; nan, nan, 0; nan, nan, 0]), Fcn("D3", "[x; y; x - z]"));
            tc.verifyTrue(symEq(val, str2sym("[x1 - z1; x2 - z2; x3 - z3; x4 - z4]")));
        end
        function ndDoFEdge3D(tc)
            D3TElem = MshEnt("D3T").msh;
            % Local edges: [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4].
            val = eval(NdDoF("D3", D3TElem, 1, 1/2, [0; 0; 0]), Fcn("D3", "z"));
            tc.verifyTrue(symEq(val, str2sym("[z1 + z2; z1 + z3; z1 + z4; z2 + z3; z2 + z4; z3 + z4]/2")));
            val = eval(NdDoF("D3", D3TElem, 1, [0, 1], [0; 0; 0]), Fcn("D3", "y"));
            tc.verifyTrue(symEq(val, str2sym("[y1, y2; y1, y3; y1, y4; y2, y3; y2, y4; y3, y4]")));
            % Tangential component of vector function on edge.
            P = tc.P3T;
            val = eval(NdDoF("D3", D3TElem, 1, 1/2, zeros(3), "coef", MshEnt("D3L").UTV), Fcn("D3", "[y; z; x]"));
            ref = zeros(6, 1);
            for iEdge = 1:6
                EgNd1 = P(:, D3TElem.edge.node(1, iEdge));
                EgNd2 = P(:, D3TElem.edge.node(2, iEdge));
                mid = (EgNd1 + EgNd2) / 2;
                ref(iEdge) = dot((EgNd2 - EgNd1) / norm(EgNd2 - EgNd1), mid([2, 3, 1]));
            end
            tc.verifyEqual(numSym(val, "D3T", P), ref, "AbsTol", 1e-12);
        end
        function ndDoFFace3D(tc)
            D3TElem = MshEnt("D3T").msh;
            % Local face i is opposite to vertex i: [2, 3, 4], [1, 4, 3], [1, 2, 4], [1, 3, 2].
            val = eval(NdDoF("D3", D3TElem, 2, [0, 1, 0; 0, 0, 1], [0; 0; 0]), Fcn("D3", "x"));
            tc.verifyTrue(symEq(val, str2sym("[x2, x3, x4; x1, x4, x3; x1, x2, x4; x1, x3, x2]")));
            val = eval(NdDoF("D3", D3TElem, 2, [1/3; 1/3], [0; 0; 0]), Fcn("D3", "x"));
            tc.verifyTrue(symEq(val, str2sym("[x2 + x3 + x4; x1 + x3 + x4; x1 + x2 + x4; x1 + x2 + x3]/3")));
            % Normal component of vector function on face.
            P = tc.P3T;
            val = eval(NdDoF("D3", D3TElem, 2, [1/3; 1/3], zeros(3), "coef", MshEnt("D3F").UNV), Fcn("D3", "[x; y; z]"));
            ref = zeros(4, 1);
            for iFace = 1:4
                FcNd = P(:, D3TElem.face.node(:, iFace));
                nor = cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1));
                ref(iFace) = dot(nor / norm(nor), mean(FcNd, 2));
                % Outward normal: opposite vertex is on the negative side.
                tc.verifyLessThan(dot(nor, P(:, iFace) - FcNd(:, 1)), 0);
            end
            tc.verifyEqual(numSym(val, "D3T", P), ref, "AbsTol", 1e-12);
        end
        function ndDoFElem3D(tc)
            D3TElem = MshEnt("D3T").msh;
            val = eval(NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], [0; 0; 0]), Fcn("D3", "x"));
            tc.verifyTrue(symEq(val, str2sym("(x1 + x2 + x3 + x4)/4")));
            val = eval(NdDoF("D3", D3TElem, 3, [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1]', [0; 0; 0]), Fcn("D3", "y"));
            tc.verifyTrue(symEq(val, str2sym("[y1, y2, y3, y4]")));
            val = eval(NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], [1; 0; 0]), Fcn("D3", "x^2 + y*z"));
            tc.verifyTrue(symEq(val, str2sym("(x1 + x2 + x3 + x4)/2")));
        end
        function ndDoFNum3D(tc)
            % Numerical evaluation on mesh agrees with symbolic evaluation.
            msh = mshD3TS([0, 1, 0, 1, 0, 1], 1);
            fcn = Fcn("D3", "x^2 + y*z - z");
            crds = {[], 1/2, [1/3; 1/3], [1/4; 1/4; 1/4]};
            for EntDim = 0:3
                ndDoF = NdDoF("D3", msh, EntDim, crds{EntDim + 1}, [0; 0; 0]);
                tc.verifyEqual(ndDoF.nEnt, msh.nEnt(EntDim));
                valNum = ndDoF.eval(fcn, "valType", "num");
                valSym = ndDoF.eval(fcn);
                tc.verifyEqual(valNum, double(valSym), "AbsTol", 1e-14);
                % Reference: function value at barycenter of mesh entity.
                fun = fcn.getFun;
                ref = zeros(ndDoF.nEnt, 1);
                for iEnt = 1:ndDoF.nEnt
                    if EntDim == 0
                        ref(iEnt) = fun(msh.node.coord(:, iEnt));
                    else
                        ref(iEnt) = fun(mean(msh.node.coord(:, msh.ent(EntDim).node(:, iEnt)), 2));
                    end
                end
                tc.verifyEqual(valNum, ref, "AbsTol", 1e-14);
            end
        end
        function moDoF3D(tc)
            D3TElem = MshEnt("D3T").msh;
            P = tc.P3T;
            vol = det(P(:, 2:4) - P(:, 1)) / 6;
            len = zeros(6, 1); area = zeros(4, 1); FcMean = zeros(4, 1);
            for iEdge = 1:6
                len(iEdge) = norm(P(:, D3TElem.edge.node(2, iEdge)) - P(:, D3TElem.edge.node(1, iEdge)));
            end
            for iFace = 1:4
                FcNd = P(:, D3TElem.face.node(:, iFace));
                area(iFace) = norm(cross(FcNd(:, 2) - FcNd(:, 1), FcNd(:, 3) - FcNd(:, 1))) / 2;
                FcMean(iFace) = mean(FcNd(1, :));
            end
            % Element.
            val = eval(MoDoF("D3", D3TElem, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 2)), Fcn("D3", "1"));
            tc.verifyEqual(numSym(val, "D3T", P), vol, "AbsTol", 1e-12);
            val = eval(MoDoF("D3", D3TElem, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 2)), Fcn("D3", "x"));
            tc.verifyEqual(numSym(val, "D3T", P), vol * mean(P(1, :)), "AbsTol", 1e-12);
            % int_K lambda_2 = vol / 4.
            val = eval(MoDoF("D3", D3TElem, 3, Fcn("D3R", "l"), [0; 0; 0], "GInt", GInt("D3T", 2)), Fcn("D3", "1"));
            tc.verifyEqual(numSym(val, "D3T", P), vol / 4, "AbsTol", 1e-12);
            % Face.
            val = eval(MoDoF("D3", D3TElem, 2, Fcn("D3R2", "1"), [0; 0; 0], "GInt", GInt("D3F", 2)), Fcn("D3", "1"));
            tc.verifyEqual(numSym(val, "D3T", P), area, "AbsTol", 1e-12);
            val = eval(MoDoF("D3", D3TElem, 2, Fcn("D3R2", "1"), [0; 0; 0], "GInt", GInt("D3F", 2)), Fcn("D3", "x"));
            tc.verifyEqual(numSym(val, "D3T", P), area .* FcMean, "AbsTol", 1e-12);
            % int_F s = area / 3.
            val = eval(MoDoF("D3", D3TElem, 2, Fcn("D3R2", "s"), [0; 0; 0], "GInt", GInt("D3F", 2)), Fcn("D3", "1"));
            tc.verifyEqual(numSym(val, "D3T", P), area / 3, "AbsTol", 1e-12);
            % Edge.
            val = eval(MoDoF("D3", D3TElem, 1, Fcn("D3R1", "1"), [0; 0; 0], "GInt", GInt("D3L", 2)), Fcn("D3", "1"));
            tc.verifyEqual(numSym(val, "D3T", P), len, "AbsTol", 1e-12);
            val = eval(MoDoF("D3", D3TElem, 1, Fcn("D3R1", "s - 1/2"), [0; 0; 0], "GInt", GInt("D3L", 2)), Fcn("D3", "1"));
            tc.verifyTrue(symEq(val, sym(zeros(6, 1))));
        end
        function moDoFNum3D(tc)
            % Numerical evaluation on Kuhn mesh of unit cube (6 tetrahedra, 18 faces, 19 edges).
            msh = mshD3TS([0, 1, 0, 1, 0, 1], 1);
            val = eval(MoDoF("D3", msh, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 1)), Fcn("D3", "1"), "valType", "num");
            tc.verifyEqual(val, ones(6, 1) / 6, "AbsTol", 1e-14);
            val = eval(MoDoF("D3", msh, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 2)), Fcn("D3", "x^2"), "valType", "num");
            tc.verifyEqual(sum(val), 1/3, "AbsTol", 1e-14);
            moDoF = MoDoF("D3", msh, 2, Fcn("D3R2", "1"), [0; 0; 0], "GInt", GInt("D3F", 1), "EntIdx", msh.bdEnt(2, 1:6));
            val = eval(moDoF, Fcn("D3", "1"), "valType", "num");
            tc.verifyEqual(sum(val), 6, "AbsTol", 1e-14);
            tc.verifyEqual(val, double(eval(moDoF, Fcn("D3", "1"))), "AbsTol", 1e-14);
            val = eval(MoDoF("D3", msh, 1, Fcn("D3R1", "1"), [0; 0; 0], "GInt", GInt("D3L", 1)), Fcn("D3", "1"), "valType", "num");
            tc.verifyEqual(sum(val), 12 + 6 * sqrt(2) + sqrt(3), "AbsTol", 1e-13);
        end
        function prop3D(tc)
            D3TElem = MshEnt("D3T").msh;
            % Default sharing: DoFs on vertices, edges and faces are shared, DoFs on elements are not.
            tc.verifyTrue(NdDoF("D3", D3TElem, 0, [], [0; 0; 0]).share);
            tc.verifyTrue(NdDoF("D3", D3TElem, 1, 1/2, [0; 0; 0]).share);
            tc.verifyTrue(NdDoF("D3", D3TElem, 2, [1/3; 1/3], [0; 0; 0]).share);
            tc.verifyFalse(NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], [0; 0; 0]).share);
            tc.verifyTrue(MoDoF("D3", D3TElem, 2, Fcn("D3R2", "1"), [0; 0; 0], "GInt", GInt("D3F", 1)).share);
            tc.verifyFalse(MoDoF("D3", D3TElem, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 1)).share);
            tc.verifyEqual(MoDoF("D3", D3TElem, 3, Fcn("D3R", "1"), [0; 0; 0], "GInt", GInt("D3T", 1)).nDoF, 1);
            % Invalid properties.
            tc.verifyError(@() NdDoF("D3", MshEnt("D2T").msh, 0, [], [0; 0; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3", D3TElem, 0, [], [0; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3", D3TElem, 2, 1/3, [0; 0; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3", D3TElem, 3, [1/2; 1/2; 1/2], [0; 0; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3", D3TElem, 1, 1/3, [0; 0; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3", D3TElem, 3, [1/4; 1/4; 1/4], [0; 0; 0], "share", true), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3", D3TElem, 2, [1/3; 1/3], [0; 0; 0], "coef", MshEnt("D3L").UTV), ...
                "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D2", D3TElem, 0, [], [0; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() MoDoF("D3", D3TElem, 2, Fcn("D3R2", "1"), [0; 0; 0], "GInt", GInt("D3T", 1)), ...
                "MATLAB:assertion:failed");
            tc.verifyError(@() MoDoF("D3", D3TElem, 2, Fcn("D3R1", "1"), [0; 0; 0], "GInt", GInt("D3F", 1)), ...
                "MATLAB:assertion:failed");
            tc.verifyError(@() MoDoF("D3", D3TElem, 1, Fcn("D3R1", "s"), [0; 0; 0], "GInt", GInt("D3L", 1)), ...
                "MATLAB:assertion:failed");
            % Function on wrong domain.
            tc.verifyError(@() eval(NdDoF("D3", D3TElem, 0, [], [0; 0; 0]), Fcn("D2", "x")), "MATLAB:assertion:failed");
        end
        function ndDoFTrace3D(tc)
            % Trace DoF on reference face.
            val = eval(NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], [0; 0]), Fcn("D3R2", "s + 2*t + 1"));
            tc.verifyTrue(symEq(val, sym([1, 2, 3])));
            % Trace DoF on faces of tetrahedron: function sampled at face vertices (global face order).
            D3TElem = MshEnt("D3T").msh;
            ndDoF = NdDoF("D3R2", D3TElem, 2, [0, 1, 0; 0, 0, 1], [0; 0]);
            tc.verifyFalse(ndDoF.share);
            val = eval(ndDoF, Fcn("D3", "x"));
            tc.verifyTrue(symEq(val, str2sym("[x2, x3, x4; x1, x4, x3; x1, x2, x4; x1, x3, x2]")));
            % Numerical evaluation on mesh.
            msh = mshD3TS([0, 1, 0, 1, 0, 1], 1);
            ndDoF = NdDoF("D3R2", msh, 2, [1/3; 1/3], [0; 0]);
            val = ndDoF.eval(Fcn("D3", "x + y*z"), "valType", "num");
            FcCen = (msh.node.coord(:, msh.face.node(1, :)) + msh.node.coord(:, msh.face.node(2, :)) + ...
                msh.node.coord(:, msh.face.node(3, :))) / 3;
            tc.verifyEqual(val, (FcCen(1, :) + FcCen(2, :) .* FcCen(3, :))', "AbsTol", 1e-14);
            % Numerical evaluation agrees with symbolic evaluation (several sample nodes, global face order).
            ndDoF = NdDoF("D3R2", msh, 2, [0, 1, 0, 1/3; 0, 0, 1, 1/3], [0; 0]);
            fcn = Fcn("D3", "x^2 + 2*y*z - z");
            tc.verifyEqual(ndDoF.eval(fcn, "valType", "num"), double(ndDoF.eval(fcn)), "AbsTol", 1e-14);
            val = eval(NdDoF("D3R2", MshEnt("D3FR").msh, 2, [0, 1, 0; 0, 0, 1], [0; 0]), Fcn("D3R2", "s + 2*t + 1"), "valType", "num");
            tc.verifyEqual(val, [1, 2, 3], "AbsTol", 1e-14);
            % Invalid: derivative, sharing, edge.
            tc.verifyError(@() NdDoF("D3R2", D3TElem, 2, [1/3; 1/3], [1; 0]), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3R2", D3TElem, 2, [1/3; 1/3], [0; 0], "share", true), "MATLAB:assertion:failed");
            tc.verifyError(@() NdDoF("D3R2", D3TElem, 1, 1/2, [0; 0]), "MATLAB:assertion:failed");
        end
        function index(tc)
            D2TElem = MshEnt("D2T").msh;
            dof = [NdDoF("D2", D2TElem, 0, [], [0; 0]), ...
                NdDoF("D2", D2TElem, 1, [1/3, 2/3], [0; 0]), ...
                NdDoF("D2", D2TElem, 2, [1/2, 1/4; 1/4, 1/2; 1/4, 1/4]', [0; 0])];
            tc.verifyEqual(dof.cumDoF, 12);
            tc.verifyEqual(dof.cumDoF(1), 3);
            tc.verifyEqual(dof.cumDoF(2), 9);
            tc.verifyEqual(dof.cumDoF(3), 12);
            tc.verifyEqual(sub2ind(dof, 1), 1:3);
            tc.verifyEqual(sub2ind(dof, 2), 4:9);
            tc.verifyEqual(sub2ind(dof, 3), 10:12);
            tc.verifyEqual(sub2ind(dof, 1, 1:3, 1), 1:3);
            tc.verifyEqual(sub2ind(dof, 2, 1:3, 1), 4:6);
            tc.verifyEqual(sub2ind(dof, 2, 1:3, 2), 7:9);
            tc.verifyEqual(sub2ind(dof, 3, 1, 1:3), 10:12);
            [iDoF, iEnt, iSamp] = ind2sub(dof, 2);
            tc.verifyEqual([iDoF, iEnt, iSamp], [1, 2, 1]);
            [iDoF, iEnt, iSamp] = ind2sub(dof, 5);
            tc.verifyEqual([iDoF, iEnt, iSamp], [2, 2, 1]);
            [iDoF, iEnt, iSamp] = ind2sub(dof, 7);
            tc.verifyEqual([iDoF, iEnt, iSamp], [2, 1, 2]);
            [iDoF, iEnt, iSamp] = ind2sub(dof, 11);
            tc.verifyEqual([iDoF, iEnt, iSamp], [3, 1, 2]);
        end
    end
end
