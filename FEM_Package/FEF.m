classdef FEF < Fcn
    % FEF: Finite Element function.
    properties
        msh Msh; % Mesh.
        elem {mustBeMember(elem, ["VOID", "D2T", "D2LR", "D3T", "D3FR"])} = "VOID"; % Element type.
        ElParm (:, :, :); % Parameter of function on each element.
        % ElParm(:, :, i): parameter for i-th element.
        ElCoef (:, :); % Coefficient of function on each element.
        % ElCoef(:, i): coefficient for i-th element.
        map {mustBeMember(map, ["none", "affine", "piolaDiv", "piolaCurl"])} = "none"; % Mapping of base functions (see `FE`).
        RefBase (1, :) Fcn; % Base function on reference element before mapping (empty if map is "none").
        RefKey (1, 1) string = ""; % Key of reference base functions (see `FE`).
    end
    properties (Dependent)
        sElParm (1, 2); % Size of parameter on one element.
        nElCoef; % Number of coefficients on one element.
    end
    methods
        % Constructor.
        function FEF = FEF(fES, DoFVal)
            arguments
                fES FES; % Finite Element space.
                DoFVal (:, 1); % DoF values of function.
            end
            assert(ismember(fES.msh.type, ["D2T", "D3T"]));
            assert(ismember(fES.elem, ["D2T", "D2LR", "D3T", "D3FR"]));
            assert(fES.nGlDoF == length(DoFVal));
            fcn = comb(fES.LcBase, sym("coef", [fES.nLcDoF, 1]));
            FEF = FEF@Fcn(fcn.domn, fcn.fun, fcn.coef);
            FEF.msh = fES.msh;
            FEF.elem = fES.elem;
            FEF.ElParm = fES.ElParm;
            % Zero entries of `Lc2Gl` (facets not selected in a trace space): coefficient 0.
            L2G = abs(fES.Lc2Gl);
            FEF.ElCoef = zeros(size(L2G));
            FEF.ElCoef(L2G > 0) = DoFVal(L2G(L2G > 0));
            FEF.ElCoef = FEF.ElCoef .* sign(fES.Lc2Gl);
            FEF.map = fES.map;
            FEF.RefBase = fES.RefBase;
            FEF.RefKey = fES.RefKey;
        end
        % Get functions.
        function sElParm = get.sElParm(FEF)
            sElParm = size(FEF.ElParm);
            sElParm = sElParm(1:2);
        end
        function nElCoef = get.nElCoef(FEF)
            nElCoef = size(FEF.ElCoef, 1);
        end
        function msh = getMsh(FEFs)
            msh = FEFs(1).msh;
            for iFEF = 2:length(FEFs)
                assert(isequal(FEFs(iFEF).msh, msh));
            end
        end
        % Public functions.
        function fcn = subt(fcn1, fcn2)
            % FEF.subt: subtract `fcn2` from `fcn1`.
            arguments (Input)
                fcn1 FEF;
                fcn2 FEF;
            end
            arguments (Output)
                fcn FEF;
            end
            assert(isEqualProp(fcn1, fcn2, "domn", "fun", "coef", "elem", "ElParm", "map"));
            fcn = fcn1;
            fcn.ElCoef = fcn1.ElCoef - fcn2.ElCoef;
        end
        function val = evalPnt(fEF, X, ElIdx, ord)
            % FEF.evalPnt: values of FE function (or of its derivative) at given points of given elements, e.g. to
            % evaluate a function of one mesh at points of another mesh. Points are not located here: `ElIdx` gives
            % the element containing each point. One function handle is generated per call: pass all points at once.
            arguments (Input)
                fEF FEF;
                X (:, :); % Points (by column, physical coordinates).
                ElIdx (1, :); % Index of the element containing each point.
                ord (:, :, :) = zeros(fEF.nVar, fEF.nFun); % Order of derivative (see `Fcn.dif`); default: value.
            end
            arguments (Output)
                val (:, :); % val(:, j): value at j-th point (components in column-major order).
            end
            % Element spaces only (not trace spaces).
            assert(isequal(fEF.elem, fEF.msh.ElDomn));
            nPnt = size(X, 2);
            assert(size(X, 1) == fEF.msh.dim && length(ElIdx) == nPnt);
            assert(all(ElIdx >= 1 & ElIdx <= fEF.msh.nElem));
            div = fEF.dif(ord);
            comps = repmat(Fcn.cst(0), 1, numel(div.fun));
            for iComp = 1:numel(div.fun)
                comps(iComp) = Fcn(div.domn, div.fun(iComp), div.coef);
            end
            funH = comps.getFun("vec", true);
            val = funH(reshape(X, [], 1, nPnt), reshape(fEF.ElParm(:, :, ElIdx), [], 1, nPnt), ...
                reshape(fEF.ElCoef(:, ElIdx), [], 1, nPnt));
            val = reshape(val, [], nPnt);
        end
    end
    % Static functions.
    methods (Static)
        function varargout = multi(FESs, DoFVal)
            % FEF.multi: construct multiple `FEF` objects.
            arguments (Input)
                FESs (1, :) FES; % Finite Element spaces.
                DoFVal (:, 1); % DoF values of function.
            end
            FEFs = cell(1, length(FESs));
            for iFES = 1:length(FESs)
                FEFs{iFES} = FEF(FESs(iFES), DoFVal(FESs.cumDoF(iFES - 1) + 1:FESs.cumDoF(iFES)));
            end
            if nargout <= 1
                varargout = {[FEFs{:}]};
            else
                varargout = FEFs(1:min(nargout, length(FEFs)));
            end
        end
    end
end
% Local functions.
function flag = isEqualProp(fcn1, fcn2, varargin)
    % isEqualProp: compare FE functions by the named properties.
    flag = true;
    for iProp = 1:length(varargin)
        prop = varargin{iProp};
        if ~isequal(fcn1.(prop), fcn2.(prop))
            flag = false;
            return;
        end
    end
end
