classdef Fcn
    % Fcn: function.
    properties
        domn {mustBeMember(domn, ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR", ...
            "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"])} = "VOID"; % Domain.
        % D: dimension.
        % T: triangle.
        % L: line.
        % R: reference.
        fun (:, :) sym; % Function.
        coef (:, 1) sym; % Coefficient.
    end
    properties (Dependent)
        var (:, 1) sym; % Variable.
        parm (:, :) sym; % Parameter.
        nVar; % Number of variables.
        sParm (1, 2); % Size of parameter.
        nFun; % Number of functions.
        nCoef; % Number of coefficients.
    end
    methods
        %% Constructor.
        function fcn = Fcn(domn, fun, coef)
            arguments
                domn {mustBeMember(domn, ["VOID", "D2", "D2R", "D2R1", "D2T", "D2TR", "D2L", "D2LR", ...
                    "D3", "D3R", "D3R2", "D3R1", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"])};
                fun (:, :) {mustBeA(fun, ["sym", "string"])};
                coef (:, 1) {mustBeA(coef, ["sym", "string"])} = sym([]);
            end
            fcn.domn = domn;
            if isstring(fun)
                fun = str2sym(fun);
            end
            fcn.fun = fun;
            if isstring(coef)
                coef = str2sym(coef);
            end
            fcn.coef = coef;
        end
        %% Get functions.
        function var = get.var(fcn)
            var = MshEnt(fcn.domn).var;
        end
        function parm = get.parm(fcn)
            parm = MshEnt(fcn.domn).parm;
        end
        function nVar = get.nVar(fcn)
            nVar = MshEnt(fcn.domn).nVar;
        end
        function sParm = get.sParm(fcn)
            sParm = MshEnt(fcn.domn).sParm;
        end
        function nFun = get.nFun(fcn)
            nFun = length(fcn.fun(:));
        end
        function nCoef = get.nCoef(fcn)
            nCoef = length(fcn.coef);
        end
        function domn = getDomn(fcns)
            % Fcn.getDomn: get domain of functions.
            domn = fcns(1).domn;
            for iFcn = 2:length(fcns)
                if ismember(MshEnt(domn), MshEnt(fcns(iFcn).domn))
                    domn = fcns(iFcn).domn;
                else
                    assert(ismember(MshEnt(fcns(iFcn).domn), MshEnt(domn)));
                end
            end
        end
        function funH = getFun(fcn, options)
            % Fcn.getFun: get function handle of `fun`.
            % Default input arguments are `var`, `parm`, `coef` (if non-empty).
            % Multi-point evaluation (vectorized input) is supported for scalar-valued function.
            arguments
                fcn Fcn;
                options.parm (1, :) cell = {}; % Manually set parameter.
                options.coef (1, :) cell = {}; % Manually set coefficient.
                options.vec logical = false; % Batch evaluation over points and mesh entities.
                % If true, each input argument `arg` (symbolic array with n entries) is given as an n x nPnt x nEnt or
                % n x 1 x nEnt numeric array, whose first dimension lists the entries of `arg` in column-major order;
                % e.g. `var` as d x nPnt x nEnt points, a 3 x 4 parameter as 12 x 1 x nEnt.
                % Output is a 1 x nPnt x nEnt array (n x 1 x nEnt arguments are replicated over points).
                % Only supported for scalar-valued function.
                % `fcn` may be an array of functions with the same arguments: output is then nFcn x nPnt x nEnt, and
                % only one function handle is generated (`matlabFunction` costs about 0.1 s per call).
            end
            if options.vec && ~isscalar(fcn)
                for iFcn = 2:numel(fcn)
                    assert(isequal(fcn(iFcn).domn, fcn(1).domn) && isequal(fcn(iFcn).coef, fcn(1).coef));
                end
                funs = arrayfun(@(f) f.fun, fcn, "UniformOutput", false);
                fcn = fcn(1);
            else
                assert(isscalar(fcn));
                funs = {fcn.fun};
            end
            argLst = {};
            if ~isempty(fcn.var)
                argLst{end + 1} = fcn.var;
            end
            if ~isempty(options.parm)
                for i = 1:length(options.parm)
                    assert(isa(options.parm{i}, "sym") || isa(options.parm{i}, "string"));
                    if isstring(options.parm{i})
                        options.parm{i} = str2sym(options.parm{i});
                    end
                end
                argLst = [argLst, options.parm];
            else
                if ~isempty(fcn.parm)
                    argLst{end + 1} = fcn.parm;
                end
            end
            if ~isempty(options.coef)
                for i = 1:length(options.coef)
                    assert(isa(options.coef{i}, "sym") || isa(options.coef{i}, "string"));
                    if isstring(options.coef{i})
                        options.coef{i} = str2sym(options.coef{i});
                    end
                end
                argLst = [argLst, options.coef];
            else
                if ~isempty(fcn.coef)
                    argLst{end + 1} = fcn.coef;
                end
            end
            if options.vec
                % Each argument is passed as a column vector: the generated code reads its k-th entry as row `in(k, :)`,
                % so a batch of values (by column) is evaluated at once.
                assert(all(cellfun(@isscalar, funs)));
                nSym = cellfun(@numel, argLst);
                for i = 1:length(argLst)
                    if ~iscolumn(argLst{i})
                        argLst{i} = reshape(argLst{i}, [], 1);
                    end
                end
                % One output per function.
                vecH = matlabFunction(funs{:}, "Vars", argLst);
                funH = @(varargin) vecEval(vecH, nSym, numel(funs), varargin);
                return;
            end
            funH = matlabFunction(fcn.fun, "Vars", argLst);
            if ~isempty(fcn.var) && all(~ismember(fcn.var, symvar(fcn.fun)))
                % If `fun` not depend on `var`, broadcast the constant output across all evaluation points.
                argStr = "";
                for i = 1:length(argLst)
                    if i > 1
                        argStr = argStr + ", ";
                    end
                    argStr = argStr + "var" + i;
                end
                funStr = "@(" + argStr + ") funH(" + argStr + ") .* ones(1, size(var1, 2))";
                funH = eval(funStr);
            end
        end
        %% Overload intrinsic operators.
        function fcn = plus(fcn1, fcn2)
            % Fcn.plus: overload + operator.
            % Coefficients are concatenated in plus operation.
            arguments (Input)
                fcn1 Fcn;
                fcn2 Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn = mrgFcn(fcn1, fcn2, "domn");
            fcn.fun = fcn1.fun + fcn2.fun;
            fcn.coef = [fcn1.coef; fcn2.coef];
        end
        function fcn = minus(fcn1, fcn2)
            % Fcn.minus: overload - operator.
            % Coefficients are concatenated in minus operation.
            arguments (Input)
                fcn1 Fcn;
                fcn2 Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn = mrgFcn(fcn1, fcn2, "domn");
            fcn.fun = fcn1.fun - fcn2.fun;
            fcn.coef = [fcn1.coef; fcn2.coef];
        end
        function fcn = uminus(fcn)
            % Fcn.uminus: overload - operator.
            fcn.fun = -fcn.fun;
        end
        function fcn = mtimes(fcn1, fcn2)
            % Fcn.mtimes: overload * operator.
            arguments (Input)
                fcn1 Fcn;
                fcn2 {mustBeA(fcn2, ["Fcn", "double"])};
            end
            arguments (Output)
                fcn Fcn;
            end
            if isnumeric(fcn2)
                fcn2 = Fcn.cst(fcn2);
            end
            fcn = mrgFcn(fcn1, fcn2, "domn", "coef");
            fcn.fun = fcn1.fun * fcn2.fun;
        end
        function fcn = mldivide(fcn1, fcn2)
            % Fcn.mldivide: overload \ operator.
            arguments (Input)
                fcn1 Fcn;
                fcn2 {mustBeA(fcn2, ["Fcn", "double"])};
            end
            arguments (Output)
                fcn Fcn;
            end
            if isnumeric(fcn2)
                fcn2 = Fcn.cst(fcn2);
            end
            fcn = mrgFcn(fcn1, fcn2, "domn", "coef");
            fcn.fun = fcn1.fun \ fcn2.fun;
        end
        function fcn = times(fcn1, fcn2)
            % Fcn.times: overload .* operator.
            arguments (Input)
                fcn1 Fcn;
                fcn2 Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn = mrgFcn(fcn1, fcn2, "domn", "coef");
            fcn.fun = fcn1.fun .* fcn2.fun;
        end
        function fcn = ldivide(fcn1, fcn2)
            % Fcn.ldivide: overload .\ operator.
            arguments (Input)
                fcn1 Fcn;
                fcn2 Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn = mrgFcn(fcn1, fcn2, "domn", "coef");
            fcn.fun = fcn1.fun .\ fcn2.fun;
        end
        function fcn = power(fcn, r)
            % Fcn.power: overload .^ operator.
            arguments (Input)
                fcn Fcn;
                r {mustBeNumeric(r)};
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn.fun = (fcn.fun) .^ r;
        end
        function fcn = transpose(fcn)
            % Fcn.transpose: overload .' operator.
            arguments (Input)
                fcn Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn.fun = fcn.fun.';
        end
        %% Overload intrinsic functions.
        function fcn = sum(fcn, dim)
            % Fcn.sum: overload `sum`.
            arguments (Input)
                fcn Fcn;
                dim = []; % Dimension to sum over.
                % If `dim` is not specified, sum over all dimensions.
            end
            arguments (Output)
                fcn Fcn;
            end
            if isempty(dim)
                fcn.fun = sum(fcn.fun, "all");
            else
                fcn.fun = sum(fcn.fun, dim);
            end
        end
        function fcn = abs(fcn)
            % Fcn.abs: overload `abs`.
            arguments (Input)
                fcn Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn.fun = abs(fcn.fun);
        end
        function fcn = dot(fcn1, fcn2)
            % Fcn.dot: overload `dot`.
            arguments (Input)
                fcn1 Fcn;
                fcn2 Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn = mrgFcn(fcn1, fcn2, "domn", "coef");
            fcn.fun = sum(fcn1.fun .* fcn2.fun);
        end
        function fcn = diff(fcn, var, ord)
            % Fcn.diff: overload `diff`.
            arguments (Input)
                fcn Fcn;
                var {mustBeA(var, ["sym", "string"])};
                ord = 1;
            end
            arguments (Output)
                fcn Fcn;
            end
            if isstring(var)
                var = str2sym(var);
            end
            fcn.fun = diff(fcn.fun, var, ord);
        end
        function fcn = simplify(fcn)
            % Fcn.simplify: overload `simplify`.
            arguments (Input)
                fcn Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn.fun = simplify(fcn.fun);
        end
        function fcn = subs(fcn, old, new)
            % Fcn.subs: overload `subs`.
            arguments (Input)
                fcn Fcn;
                old {mustBeA(old, ["sym", "string"])};
                new {mustBeA(new, ["double", "sym", "string"])};
            end
            arguments (Output)
                fcn Fcn;
            end
            if isnumeric(old)
                old = sym(old);
            end
            if isstring(old)
                old = str2sym(old);
            end
            if isstring(new)
                new = str2sym(new);
            end
            fcn.fun = simplify(subs(fcn.fun, old, new));
        end
        %% Mathematical operations.
        function funVal = eval(fcn, val)
            % Fcn.eval: evaluate function at given value.
            arguments
                fcn Fcn;
                val (:, 1) {mustBeA(val, ["double", "sym", "string"])};
            end
            if isstring(val)
                val = str2sym(val);
            end
            assert(length(val) == fcn.nVar);
            funVal = simplify(subs(fcn.fun, fcn.var, val));
        end
        function fcn = comb(fcns, coef)
            % Fcn.comb: linear combination of functions w.r.t. given coefficients.
            arguments (Input)
                fcns (1, :) Fcn;
                coef (:, 1) {mustBeA(coef, ["double", "sym", "string"])};
            end
            arguments (Output)
                fcn Fcn;
            end
            if isstring(coef)
                coef = str2sym(coef);
            end
            fcn = Fcn.cst(0);
            for iFcn = 1:length(fcns)
                fcn = mrgFcn(fcn, fcns(iFcn), "domn");
            end
            for iFcn = 1:length(fcns)
                fcn.fun = fcn.fun + fcns(iFcn).fun * coef(iFcn);
                fcn.coef = [fcn.coef; fcns(iFcn).coef];
            end
            if ~isnumeric(coef)
                assert(isempty(fcn.coef));
                fcn.coef = coef;
            end
        end
        function fcns = tfm(fcns, domn)
            % Fcn.tfm: transform function to original/reference domain.

            % Valid input arguments:
            % fcn.domn | tfm.domn
            %----------+-----------
            % D2       | D2TR D2LR
            % D2R      | D2T
            % D2T      | D2TR
            % D2L      | D2LR
            % D2TR     | D2T
            % D3       | D3TR D3FR D3LR
            % D3R      | D3T
            % D3T      | D3TR
            % D3F      | D3FR
            % D3L      | D3LR
            % D3TR     | D3T

            arguments (Input)
                fcns Fcn;
                domn {mustBeMember(domn, ["D2T", "D2TR", "D2LR", "D3T", "D3TR", "D3FR", "D3LR"])}; % Target domain.
            end
            arguments (Output)
                fcns Fcn;
            end
            % Transformed functions are plain `Fcn` objects, even if input functions are of a subclass (e.g. `FEF`).
            tfmFcns = repmat(Fcn.cst(0), size(fcns));
            for iFcn = 1:length(fcns)
                fcn = fcns(iFcn);
                switch domn
                    case "D2T"
                        assert(ismember(fcn.domn, ["VOID", "D2R", "D2TR"]));
                        fcn = fcn.ensDomn("D2R");
                        fcn = fcn.compose(Tfm("D2T").refTfm);
                    case "D2TR"
                        assert(ismember(fcn.domn, ["VOID", "D2", "D2T"]));
                        fcn = fcn.ensDomn("D2");
                        fcn = fcn.compose(Tfm("D2T").orgTfm);
                    case "D2LR"
                        assert(ismember(fcn.domn, ["VOID", "D2", "D2L"]));
                        fcn = fcn.ensDomn("D2");
                        fcn = fcn.compose(Tfm("D2L").orgTfm);
                    case "D3T"
                        assert(ismember(fcn.domn, ["VOID", "D3R", "D3TR"]));
                        fcn = fcn.ensDomn("D3R");
                        fcn = fcn.compose(Tfm("D3T").refTfm);
                    case "D3TR"
                        assert(ismember(fcn.domn, ["VOID", "D3", "D3T"]));
                        fcn = fcn.ensDomn("D3");
                        fcn = fcn.compose(Tfm("D3T").orgTfm);
                    case "D3FR"
                        assert(ismember(fcn.domn, ["VOID", "D3", "D3F"]));
                        fcn = fcn.ensDomn("D3");
                        fcn = fcn.compose(Tfm("D3F").orgTfm);
                    case "D3LR"
                        assert(ismember(fcn.domn, ["VOID", "D3", "D3L"]));
                        fcn = fcn.ensDomn("D3");
                        fcn = fcn.compose(Tfm("D3L").orgTfm);
                end
                tfmFcns(iFcn) = fcn;
            end
            fcns = tfmFcns;
        end
        function fcns = dif(fcns, ord)
            % Fcn.dif: differentiate functions w.r.t. given order.
            arguments (Input)
                fcns (1, :) Fcn;
                ord (:, :, :); % ord: order of derivative.
                % Basic case: `ord` is 2D
                % - ord(i, j) = k: take k-th derivative of j-th function component w.r.t. i-th variable;
                % - ord(i, j) is not non-negative number: set j-th function component to 0.
                % Stacking Case: `ord` is 3D
                % - The 3rd dimension represents a collection of differential operations. Each slice ord(:, :, m) follows the 2D logic above.
                % - If input function is scalar/vector, the output is vector/matrix.
            end
            arguments (Output)
                fcns (1, :) Fcn;
            end
            for iFcn = 1:length(fcns)
                fcn = fcns(iFcn);
                assert(size(ord, 1) == fcn.nVar && size(ord, 2) == fcn.nFun);
                if ndims(ord) == 2
                    for iFun = 1:fcn.nFun
                        if all(ord(:, iFun) >= 0)
                            for iVar = 1:fcn.nVar
                                fcn.fun(iFun) = diff(fcn.fun(iFun), fcn.var(iVar), ord(iVar, iFun));
                            end
                        else
                            fcn.fun(iFun) = sym(0);
                            continue;
                        end
                    end
                elseif ndims(ord) == 3
                    if isscalar(fcn.fun)
                        nOutFun = size(ord, 3);
                        outFun = repmat(fcn.fun, nOutFun, 1);
                        for iOutFun = 1:nOutFun
                            if all(ord(:, 1, iOutFun) >= 0)
                                for iVar = 1:fcn.nVar
                                    outFun(iOutFun) = diff(outFun(iOutFun), fcn.var(iVar), ord(iVar, 1, iOutFun));
                                end
                            else
                                outFun(iOutFun) = sym(0);
                                continue;
                            end
                        end
                        fcn.fun = outFun;
                    elseif iscolumn(fcn.fun)
                        nOutFun = size(ord, 3);
                        outFun = repmat(fcn.fun, 1, nOutFun);
                        for iOutFun = 1:nOutFun
                            for iFun = 1:fcn.nFun
                                if all(ord(:, iFun, iOutFun) >= 0)
                                    for iVar = 1:fcn.nVar
                                        outFun(iFun, iOutFun) = diff(outFun(iFun, iOutFun), fcn.var(iVar), ord(iVar, iFun, iOutFun));
                                    end
                                else
                                    outFun(iFun, iOutFun) = sym(0);
                                    continue;
                                end
                            end
                        end
                        fcn.fun = outFun;
                    else
                        error('Invalid size of ord.');
                    end
                else
                    error('Invalid size of ord.');
                end
                fcns(iFcn) = fcn;
            end
        end
        function intVal = int(fcn, domn, idx)
            % Fcn.int: integrate function over given domain.

            % Valid input arguments and corresponding output:
            % fcn.domn | int.domn | intVal.domn
            % ---------+----------+-----------
            % D2       | D2T      | D2T
            % D2       | D2L      | D2L
            % D2R      | D2TR     | D2R
            % D2R      | D2LR     | D2R
            % D2R1     | D2LR     | D2R1
            % D2T      | D2T      | D2T
            % D2T      | D2L      | D2T
            % D2L      | D2L      | D2L
            % D2TR     | D2TR     | D2TR
            % D2TR     | D2LR     | D2TR
            % D2LR     | D2LR     | D2LR
            % D3       | D3T      | D3T
            % D3       | D3F      | D3F
            % D3       | D3L      | D3L
            % D3R      | D3TR     | D3R
            % D3R      | D3FR     | D3R
            % D3R      | D3LR     | D3R
            % D3R2     | D3FR     | D3R2
            % D3R1     | D3LR     | D3R1
            % D3T      | D3T      | D3T
            % D3T      | D3F      | D3T
            % D3T      | D3L      | D3T
            % D3F      | D3F      | D3F
            % D3L      | D3L      | D3L
            % D3TR     | D3TR     | D3TR
            % D3TR     | D3FR     | D3TR
            % D3TR     | D3LR     | D3TR
            % D3FR     | D3FR     | D3FR
            % D3LR     | D3LR     | D3LR

            % Input arguments:
            arguments (Input)
                fcn Fcn;
                domn {mustBeMember(domn, ["D2T", "D2TR", "D2L", "D2LR", "D3T", "D3TR", "D3F", "D3FR", "D3L", "D3LR"])}; % Integrated domain.
                idx = []; % When integrating a function over a sub-entity of its original domain, specify the index of the sub-entity.
                % For example, when integrating a function defined on a triangle over its 2-nd edge, set `idx` to 2;
                % when integrating a function defined on a tetrahedron over its 3-rd face (edge), set `idx` to 3.
            end
            arguments (Output)
                intVal Fcn; % Value of integral.
                % `intVal` is a `Fcn` object where `fun` is independent of `var`.
            end
            switch domn
                case "D2T"
                    assert(ismember(fcn.domn, ["VOID", "D2", "D2T"]));
                    fcn = fcn.ensDomn("D2");
                    intFcn = fcn.tfm("D2TR") .* Tfm("D2T").JDet;
                    assert(isequal(MshEnt("D2TR").node.coord, sym([0, 1, 0; 0, 0, 1])));
                    intVal = Fcn("D2T", int(int(intFcn.fun, intFcn.var(2), 0, 1 - intFcn.var(1)), intFcn.var(1), 0, 1), fcn.coef);
                case "D2TR"
                    assert(ismember(fcn.domn, ["VOID", "D2R", "D2TR"]));
                    fcn = fcn.ensDomn("D2R");
                    assert(isequal(MshEnt("D2TR").node.coord, sym([0, 1, 0; 0, 0, 1])));
                    intVal = Fcn(fcn.domn, int(int(fcn.fun, fcn.var(2), 0, 1 - fcn.var(1)), fcn.var(1), 0, 1), fcn.coef);
                case "D2L"
                    assert(ismember(fcn.domn, ["VOID", "D2", "D2L", "D2T"]));
                    switch fcn.domn
                        case {"VOID", "D2", "D2L"}
                            fcn = fcn.ensDomn("D2");
                            intFcn = fcn.tfm("D2LR") .* Tfm("D2L").JNorm;
                            assert(isequal(MshEnt("D2LR").node.coord, sym([0, 1])));
                            intVal = Fcn("D2L", int(intFcn.fun, intFcn.var, 0, 1), fcn.coef);
                        case "D2T"
                            assert(~isempty(idx));
                            mshEnt = MshEnt("D2T"); EgParm = mshEnt.node.coord(:, mshEnt.edge.node(:, idx));
                            orgTfm = Tfm("D2L").orgTfm.subParm(EgParm);
                            JNorm = Tfm("D2L").JNorm.subParm(EgParm);
                            intFcn = fcn.compose(orgTfm) .* JNorm;
                            assert(isequal(MshEnt("D2LR").node.coord, sym([0, 1])));
                            intVal = Fcn("D2T", int(intFcn.fun, intFcn.var, 0, 1), fcn.coef);
                    end
                case "D2LR"
                    assert(ismember(fcn.domn, ["VOID", "D2R1", "D2LR", "D2R", "D2TR"]));
                    switch fcn.domn
                        case {"VOID", "D2R1", "D2LR"}
                            fcn = fcn.ensDomn("D2R1");
                            assert(isequal(MshEnt("D2LR").node.coord, sym([0, 1])));
                            intVal = Fcn(fcn.domn, int(fcn.fun, fcn.var, 0, 1), fcn.coef);
                        case {"D2R", "D2TR"}
                            assert(~isempty(idx));
                            mshEnt = MshEnt("D2TR"); EgParm = mshEnt.node.coord(:, mshEnt.edge.node(:, idx));
                            orgTfm = Tfm("D2L").orgTfm.subParm(EgParm);
                            JNorm = Tfm("D2L").JNorm.subParm(EgParm);
                            intFcn = fcn.compose(orgTfm) .* JNorm;
                            assert(isequal(MshEnt("D2LR").node.coord, sym([0, 1])));
                            intVal = Fcn(fcn.domn, int(intFcn.fun, intFcn.var, 0, 1), fcn.coef);
                    end
                case "D3T"
                    assert(ismember(fcn.domn, ["VOID", "D3", "D3T"]));
                    fcn = fcn.ensDomn("D3");
                    intFcn = fcn.tfm("D3TR") .* Tfm("D3T").JDet;
                    assert(isequal(MshEnt("D3TR").node.coord, sym([0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1])));
                    intVal = Fcn("D3T", int(int(int(intFcn.fun, intFcn.var(3), 0, 1 - intFcn.var(1) - intFcn.var(2)), ...
                        intFcn.var(2), 0, 1 - intFcn.var(1)), intFcn.var(1), 0, 1), fcn.coef);
                case "D3TR"
                    assert(ismember(fcn.domn, ["VOID", "D3R", "D3TR"]));
                    fcn = fcn.ensDomn("D3R");
                    assert(isequal(MshEnt("D3TR").node.coord, sym([0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1])));
                    intVal = Fcn(fcn.domn, int(int(int(fcn.fun, fcn.var(3), 0, 1 - fcn.var(1) - fcn.var(2)), ...
                        fcn.var(2), 0, 1 - fcn.var(1)), fcn.var(1), 0, 1), fcn.coef);
                case "D3F"
                    assert(ismember(fcn.domn, ["VOID", "D3", "D3F", "D3T"]));
                    switch fcn.domn
                        case {"VOID", "D3", "D3F"}
                            fcn = fcn.ensDomn("D3");
                            intFcn = fcn.tfm("D3FR") .* Tfm("D3F").JNorm;
                            assert(isequal(MshEnt("D3FR").node.coord, sym([0, 1, 0; 0, 0, 1])));
                            intVal = Fcn("D3F", int(int(intFcn.fun, intFcn.var(2), 0, 1 - intFcn.var(1)), intFcn.var(1), 0, 1), fcn.coef);
                        case "D3T"
                            assert(~isempty(idx));
                            mshEnt = MshEnt("D3T"); FcParm = mshEnt.node.coord(:, mshEnt.face.node(:, idx));
                            orgTfm = Tfm("D3F").orgTfm.subParm(FcParm);
                            JNorm = Tfm("D3F").JNorm.subParm(FcParm);
                            intFcn = fcn.compose(orgTfm) .* JNorm;
                            assert(isequal(MshEnt("D3FR").node.coord, sym([0, 1, 0; 0, 0, 1])));
                            intVal = Fcn("D3T", int(int(intFcn.fun, intFcn.var(2), 0, 1 - intFcn.var(1)), intFcn.var(1), 0, 1), fcn.coef);
                    end
                case "D3FR"
                    assert(ismember(fcn.domn, ["VOID", "D3R2", "D3FR", "D3R", "D3TR"]));
                    switch fcn.domn
                        case {"VOID", "D3R2", "D3FR"}
                            fcn = fcn.ensDomn("D3R2");
                            assert(isequal(MshEnt("D3FR").node.coord, sym([0, 1, 0; 0, 0, 1])));
                            intVal = Fcn(fcn.domn, int(int(fcn.fun, fcn.var(2), 0, 1 - fcn.var(1)), fcn.var(1), 0, 1), fcn.coef);
                        case {"D3R", "D3TR"}
                            assert(~isempty(idx));
                            mshEnt = MshEnt("D3TR"); FcParm = mshEnt.node.coord(:, mshEnt.face.node(:, idx));
                            orgTfm = Tfm("D3F").orgTfm.subParm(FcParm);
                            JNorm = Tfm("D3F").JNorm.subParm(FcParm);
                            intFcn = fcn.compose(orgTfm) .* JNorm;
                            assert(isequal(MshEnt("D3FR").node.coord, sym([0, 1, 0; 0, 0, 1])));
                            intVal = Fcn(fcn.domn, int(int(intFcn.fun, intFcn.var(2), 0, 1 - intFcn.var(1)), intFcn.var(1), 0, 1), fcn.coef);
                    end
                case "D3L"
                    assert(ismember(fcn.domn, ["VOID", "D3", "D3L", "D3T"]));
                    switch fcn.domn
                        case {"VOID", "D3", "D3L"}
                            fcn = fcn.ensDomn("D3");
                            intFcn = fcn.tfm("D3LR") .* Tfm("D3L").JNorm;
                            assert(isequal(MshEnt("D3LR").node.coord, sym([0, 1])));
                            intVal = Fcn("D3L", int(intFcn.fun, intFcn.var, 0, 1), fcn.coef);
                        case "D3T"
                            assert(~isempty(idx));
                            mshEnt = MshEnt("D3T"); EgParm = mshEnt.node.coord(:, mshEnt.edge.node(:, idx));
                            orgTfm = Tfm("D3L").orgTfm.subParm(EgParm);
                            JNorm = Tfm("D3L").JNorm.subParm(EgParm);
                            intFcn = fcn.compose(orgTfm) .* JNorm;
                            assert(isequal(MshEnt("D3LR").node.coord, sym([0, 1])));
                            intVal = Fcn("D3T", int(intFcn.fun, intFcn.var, 0, 1), fcn.coef);
                    end
                case "D3LR"
                    assert(ismember(fcn.domn, ["VOID", "D3R1", "D3LR", "D3R", "D3TR"]));
                    switch fcn.domn
                        case {"VOID", "D3R1", "D3LR"}
                            fcn = fcn.ensDomn("D3R1");
                            assert(isequal(MshEnt("D3LR").node.coord, sym([0, 1])));
                            intVal = Fcn(fcn.domn, int(fcn.fun, fcn.var, 0, 1), fcn.coef);
                        case {"D3R", "D3TR"}
                            assert(~isempty(idx));
                            mshEnt = MshEnt("D3TR"); EgParm = mshEnt.node.coord(:, mshEnt.edge.node(:, idx));
                            orgTfm = Tfm("D3L").orgTfm.subParm(EgParm);
                            JNorm = Tfm("D3L").JNorm.subParm(EgParm);
                            intFcn = fcn.compose(orgTfm) .* JNorm;
                            assert(isequal(MshEnt("D3LR").node.coord, sym([0, 1])));
                            intVal = Fcn(fcn.domn, int(intFcn.fun, intFcn.var, 0, 1), fcn.coef);
                    end
            end
        end
        %% Parameter and coefficient management.
        function fcns = subParm(fcns, parm)
            % Fcn.subParm: substitute parameter.
            % Fcn.clrParm is automatically called after substitution.
            arguments (Input)
                fcns Fcn;
                parm (:, :) {mustBeA(parm, ["sym", "string"])};
            end
            arguments (Output)
                fcns Fcn;
            end
            if isstring(parm)
                parm = str2sym(parm);
            end
            for iFcn = 1:length(fcns)
                fcn = fcns(iFcn);
                if ~isempty(fcn.parm)
                    assert(isequal(size(parm), fcn.sParm));
                    fcn.fun = subs(fcn.fun, fcn.parm, parm);
                    fcn = fcn.clrParm;
                end
                fcns(iFcn) = fcn;
            end
        end
        function fcn = subCoef(fcn, coef)
            % Fcn.subCoef: substitute coefficient.
            % Fcn.clrCoef is automatically called after substitution.
            arguments (Input)
                fcn Fcn;
                coef (:, 1) {mustBeA(coef, ["sym", "string"])};
            end
            arguments (Output)
                fcn Fcn;
            end
            if isstring(coef)
                coef = str2sym(coef);
            end
            assert(length(coef) == fcn.nCoef);
            fcn.fun = subs(fcn.fun, fcn.coef, coef);
            fcn = fcn.clrCoef;
        end
        function fcn = clrParm(fcn)
            % Fcn.clrParm: clear parameter.
            arguments (Input)
                fcn Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn.domn = MshEnt.getInfo(fcn.domn).free;
        end
        function fcn = clrCoef(fcn)
            % Fcn.clrCoef: clear coefficient.
            arguments (Input)
                fcn Fcn;
            end
            arguments (Output)
                fcn Fcn;
            end
            fcn.coef = sym([]);
        end
    end
    %% Static functions.
    methods (Static)
        function fcn = cst(val)
            % Fcn.cst: constant function.
            fcn = Fcn("VOID", sym(val));
        end
    end
    %% Private functions.
    methods (Access = private)
        function fcn = ensDomn(fcn, domn)
            % Fcn.ensDomn: ensure domain compatibility.
            if ismember(MshEnt(fcn.domn), MshEnt(domn))
                fcn = Fcn(domn, fcn.fun);
            else
                assert(ismember(MshEnt(domn), MshEnt(fcn.domn)));
            end
        end
        function fcn = compose(fcn1, fcn2)
            % Fcn.compose: compose functions.
            % fcn = fcn1(fcn2).
            fcn = Fcn(fcn2.domn, subs(fcn1.fun, fcn1.var, fcn2.fun), fcn1.coef);
        end
    end
end
%% Local functions.
function val = vecEval(vecH, nSym, nFun, args)
    % vecEval: evaluate function handle generated with column-vector arguments on batched arguments
    % (see `Fcn.getFun` with "vec"): each argument n x nPnt x nEnt or n x 1 x nEnt, output nFun x nPnt x nEnt.
    assert(length(args) == length(nSym));
    nPnt = max(cellfun(@(arg) size(arg, 2), args));
    nEnt = max(cellfun(@(arg) size(arg, 3), args));
    inp = cell(1, length(args));
    for iArg = 1:length(args)
        arg = args{iArg};
        assert(size(arg, 1) == nSym(iArg));
        if size(arg, 2) < nPnt
            arg = repmat(arg, 1, nPnt / size(arg, 2), 1);
        end
        inp{iArg} = reshape(arg, nSym(iArg), nPnt * nEnt);
    end
    out = cell(1, nFun);
    [out{:}] = vecH(inp{:});
    val = zeros(nFun, nPnt * nEnt);
    for iFun = 1:nFun
        val(iFun, :) = out{iFun};
    end
    val = reshape(val, nFun, nPnt, nEnt);
end
function fcn = mrgFcn(fcn1, fcn2, varargin)
    % mrgFcn: merge functions.
    % varargin: names of properties to be merged.
    fcn = Fcn.cst(0);
    for iProp = 1:length(varargin)
        prop = varargin{iProp};
        fcn.(prop) = mrgProp(fcn1.(prop), fcn2.(prop), prop);
    end
end
function prop = mrgProp(prop1, prop2, name)
    % mrgProp: merge properties.
    switch name
        case "domn"
            if ismember(MshEnt(prop1), MshEnt(prop2))
                prop = prop2;
            elseif ismember(MshEnt(prop2), MshEnt(prop1))
                prop = prop1;
            else
                error('Fcn.%s are not compatible.', name);
            end
        otherwise
            % Either prop1 or prop2 is empty, or prop1 is equal to prop2.
            if isempty(prop1) && isempty(prop2)
                prop = [];
            elseif isempty(prop1)
                prop = prop2;
            elseif isempty(prop2)
                prop = prop1;
            else
                if isequal(prop1, prop2)
                    prop = prop1;
                else
                    error('Fcn.%s are not identical.', name);
                end
            end
    end
end
