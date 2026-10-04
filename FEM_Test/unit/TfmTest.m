classdef TfmTest < matlab.unittest.TestCase
    % TfmTest: unit tests of Tfm.
    methods (Test, TestTags = {'Fast'})
        function D2T(tc)
            tfm = Tfm("D2T");
            tc.verifyTrue(symEq(tfm.orgVar, str2sym("[x; y]")));
            tc.verifyTrue(symEq(tfm.refVar, str2sym("[l; m]")));
            % Reference vertices are mapped to original vertices.
            orgTfm = tfm.orgTfm;
            tc.verifyEqual(orgTfm.domn, "D2TR");
            tc.verifyTrue(symEq(orgTfm.eval([0; 0]), str2sym("[x1; y1]")));
            tc.verifyTrue(symEq(orgTfm.eval([1; 0]), str2sym("[x2; y2]")));
            tc.verifyTrue(symEq(orgTfm.eval([0; 1]), str2sym("[x3; y3]")));
            % Original vertices are mapped to reference vertices.
            refTfm = tfm.refTfm;
            tc.verifyEqual(refTfm.domn, "D2T");
            tc.verifyTrue(symEq(refTfm.eval(str2sym("[x1; y1]")), sym([0; 0])));
            tc.verifyTrue(symEq(refTfm.eval(str2sym("[x2; y2]")), sym([1; 0])));
            tc.verifyTrue(symEq(refTfm.eval(str2sym("[x3; y3]")), sym([0; 1])));
            % `toRef` is inverse of `toOrg`.
            tc.verifyTrue(symEq(subs(tfm.toRef, tfm.orgVar, tfm.toOrg), tfm.refVar));
            % Jacobian determinant.
            tc.verifyEqual(tfm.JDet.domn, "D2TR");
            tc.verifyTrue(symEq(tfm.JDet.fun, str2sym("(x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)")));
        end
        function D2L(tc)
            tfm = Tfm("D2L");
            tc.verifyTrue(symEq(tfm.orgVar, str2sym("[x; y]")));
            tc.verifyTrue(symEq(tfm.refVar, str2sym("s")));
            orgTfm = tfm.orgTfm;
            tc.verifyEqual(orgTfm.domn, "D2LR");
            tc.verifyTrue(symEq(orgTfm.eval(0), str2sym("[x1; y1]")));
            tc.verifyTrue(symEq(orgTfm.eval(1/2), str2sym("[(x1 + x2)/2; (y1 + y2)/2]")));
            tc.verifyTrue(symEq(orgTfm.eval(1), str2sym("[x2; y2]")));
            % Jacobian norm is length of line.
            tc.verifyEqual(tfm.JNorm.domn, "D2LR");
            tc.verifyTrue(symEq(tfm.JNorm.fun, str2sym("sqrt((x2 - x1)^2 + (y2 - y1)^2)")));
            tc.verifyEqual(numSym(tfm.JNorm.fun, "D2L", [0, 3; 0, 4]), 5, "AbsTol", 1e-12);
            tc.verifyEqual(numSym(tfm.JNorm.fun, "D2L", [1, -1; 2, 0]), sqrt(8), "AbsTol", 1e-12);
        end
        function D3T(tc)
            P = [0.1, 1.2, 0.3, 0.2; -0.1, 0.2, 1.1, 0.3; 0.05, 0.1, 0.2, 1.3];
            tfm = Tfm("D3T");
            tc.verifyTrue(symEq(tfm.orgVar, str2sym("[x; y; z]")));
            tc.verifyTrue(symEq(tfm.refVar, str2sym("[l; m; n]")));
            % Reference vertices are mapped to original vertices.
            orgTfm = tfm.orgTfm;
            tc.verifyEqual(orgTfm.domn, "D3TR");
            tc.verifyTrue(symEq(orgTfm.eval([0; 0; 0]), str2sym("[x1; y1; z1]")));
            tc.verifyTrue(symEq(orgTfm.eval([1; 0; 0]), str2sym("[x2; y2; z2]")));
            tc.verifyTrue(symEq(orgTfm.eval([0; 1; 0]), str2sym("[x3; y3; z3]")));
            tc.verifyTrue(symEq(orgTfm.eval([0; 0; 1]), str2sym("[x4; y4; z4]")));
            % Original vertices are mapped to reference vertices.
            refTfm = tfm.refTfm;
            tc.verifyEqual(refTfm.domn, "D3T");
            RfNd = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
            for iNode = 1:4
                tc.verifyEqual(numSym(subs(refTfm.fun, tfm.orgVar, P(:, iNode)), "D3T", P), RfNd(:, iNode), "AbsTol", 1e-12);
            end
            % `toRef` is inverse of `toOrg`.
            refPnt = [0.2; 0.3; 0.1];
            tc.verifyEqual(numSym(subs(tfm.toRef, tfm.orgVar, subs(tfm.toOrg, tfm.refVar, refPnt)), "D3T", P), refPnt, "AbsTol", 1e-12);
            % Jacobian determinant is 6 times volume.
            tc.verifyEqual(tfm.JDet.domn, "D3TR");
            tc.verifyEqual(numSym(tfm.JDet.fun, "D3T", P), det(P(:, 2:4) - P(:, 1)), "AbsTol", 1e-12);
        end
        function D3F(tc)
            P = [0.1, 1.2, 0.3; -0.1, 0.2, 1.1; 0.05, 0.1, 0.9];
            tfm = Tfm("D3F");
            tc.verifyTrue(symEq(tfm.refVar, str2sym("[s; t]")));
            orgTfm = tfm.orgTfm;
            tc.verifyEqual(orgTfm.domn, "D3FR");
            tc.verifyTrue(symEq(orgTfm.eval([0; 0]), str2sym("[x1; y1; z1]")));
            tc.verifyTrue(symEq(orgTfm.eval([1; 0]), str2sym("[x2; y2; z2]")));
            tc.verifyTrue(symEq(orgTfm.eval([0; 1]), str2sym("[x3; y3; z3]")));
            % Jacobian norm is 2 times area.
            tc.verifyEqual(tfm.JNorm.domn, "D3FR");
            tc.verifyEqual(numSym(tfm.JNorm.fun, "D3F", P), norm(cross(P(:, 2) - P(:, 1), P(:, 3) - P(:, 1))), "AbsTol", 1e-12);
            % In plane z = 0, Jacobian norm equals Jacobian determinant of D2T (counter-clockwise triangle).
            P2 = [0.1, 1.3, 0.4; -0.2, 0.3, 1.1];
            tc.verifyEqual(numSym(tfm.JNorm.fun, "D3F", [P2; 0, 0, 0]), numSym(Tfm("D2T").JDet.fun, "D2T", P2), "AbsTol", 1e-12);
        end
        function D3L(tc)
            tfm = Tfm("D3L");
            orgTfm = tfm.orgTfm;
            tc.verifyEqual(orgTfm.domn, "D3LR");
            tc.verifyTrue(symEq(orgTfm.eval(0), str2sym("[x1; y1; z1]")));
            tc.verifyTrue(symEq(orgTfm.eval(1), str2sym("[x2; y2; z2]")));
            tc.verifyEqual(numSym(tfm.JNorm.fun, "D3L", [0, 1; 0, 2; 0, 2]), 3, "AbsTol", 1e-12);
        end
    end
end
