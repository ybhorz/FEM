classdef GInt
    % GInt: Gauss integration.
    properties
        domn {mustBeMember(domn, ["VOID", "D2T", "D2L"])} = "VOID"; % Domain of integration.
        ord {mustBeMember(ord, [0, 1, 2, 3, 4, 5, 6, 7])}; % Order of integration.
        % - algebraic precision: the rule integrates polynomials up to this degree exactly.
        pnt (:, :); % Coordinates of Gauss points in reference domain (by column).
        wgt (1, :); % Weights of Gauss points in reference domain.
    end
    properties (Access = private)
        orgTfm; % Transformation to original domain.
        Jac; % Jacobian of transformation.
    end
    methods
        % Constructor.
        function GInt = GInt(domn, ord)
            arguments
                domn {mustBeMember(domn, ["D2T", "D2L"])};
                ord {mustBeMember(ord, [0, 1, 2, 3, 4, 5, 6, 7])};
            end
            GInt.domn = domn;
            GInt.ord = ord;
            switch domn
                case "D2T"
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
                    GInt.Jac = Tfm(GInt.domn).JDet.getFun;
                case "D2L"
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
            switch gInt.domn
                case "D2T"
                    if ~isempty(parm)
                        assert(isequal(size(parm), [2, 3]));
                        val = sum(fun(gInt.orgTfm(gInt.pnt, parm)) .* gInt.wgt) .* gInt.Jac([0; 0], parm);
                    else
                        val = sum(fun(gInt.pnt) .* gInt.wgt);
                    end
                case "D2L"
                    if ~isempty(parm)
                        assert(isequal(size(parm), [2, 2]));
                        val = sum(fun(gInt.orgTfm(gInt.pnt, parm)) .* gInt.wgt) .* gInt.Jac(0, parm);
                    else
                        val = sum(fun(gInt.pnt) .* gInt.wgt);
                    end
                    
            end
        end
    end
end
