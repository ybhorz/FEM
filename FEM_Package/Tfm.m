classdef Tfm
    % Tfm: transformation between original and reference domain.
    properties
        domn {mustBeMember(domn, ["VOID", "D2T", "D2L"])} = "VOID"; % Domain of transformation.
        % D: dimension.
        % T: triangle.
        % L: line.
        orgVar (:, 1) sym; % Original variable.
        refVar (:, 1) sym; % Reference variable.
        toOrg (:, 1) sym; % Transformation to original domain.
        toRef (:, 1) sym; % Transformation to reference domain.
    end
    properties (Dependent)
        JDet Fcn; % Determinant of Jacobian of transformation.
        JNorm Fcn; % Norm of Jacobian of transformation.
        orgTfm Fcn; % Transformation to original domain.
        refTfm Fcn; % Transformation to reference domain.
    end
    methods
        % Constructor.
        function tfm = Tfm(domn)
            arguments
                domn {mustBeMember(domn, ["D2T", "D2L"])};
            end
            tfm.domn = domn;
            tfm.orgVar = MshEnt(domn).var;
            tfm.refVar = MshEnt(dual(domn)).var;
            parm = MshEnt(domn).parm;
            switch domn
                case "D2T"
                    assert(isequal(MshEnt(dual(domn)).node.coord, sym([0, 1, 0; 0, 0, 1])));
                    B = parm(:, 2:3) - parm(:, 1);
                    b = parm(:, 1);
                    tfm.toOrg = B * tfm.refVar + b;
                    tfm.toRef = B \ (tfm.orgVar - b);
                case "D2L"
                    assert(isequal(MshEnt(dual(domn)).node.coord, sym([0, 1])));
                    tfm.toOrg = (parm(:, 2) - parm(:, 1)) * tfm.refVar + parm(:, 1);
            end
        end
        % Get functions.
        function JDet = get.JDet(tfm)
            JDet = Fcn(dual(tfm.domn), det(jacobian(tfm.toOrg, tfm.refVar)));
        end
        function JNorm = get.JNorm(tfm)
            JNorm = Fcn(dual(tfm.domn), norm(jacobian(tfm.toOrg, tfm.refVar)));
        end
        function refTfm = get.refTfm(tfm)
            refTfm = Fcn(tfm.domn, tfm.toRef);
        end
        function orgTfm = get.orgTfm(tfm)
            orgTfm = Fcn(dual(tfm.domn), tfm.toOrg);
        end
    end
end
% Local functions.
function dlDomn = dual(domn)
    % dual: dual domain.
    dlDomn = MshEnt(domn).dual.type;
end
