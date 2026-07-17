classdef FEF < Fcn
    % FEF: Finite Element function.
    properties
        msh Msh; % Mesh.
        elem {mustBeMember(elem, ["VOID", "D2T", "D2LR"])} = "VOID"; % Element type.
        ElParm (:, :, :); % Parameter of function on each element.
        % ElParm(:, :, i): parameter for i-th element.
        ElCoef (:, :); % Coefficient of function on each element.
        % ElCoef(:, i): coefficient for i-th element.
    end
    properties (Dependent)
        sElParm (1, 2); % Size of parameter on element.
        nElCoef; % Number of coefficients on element.
    end
    methods
        % Constructor.
        function FEF = FEF(fES, DoFVal)
            arguments
                fES FES; % Finite Element space.
                DoFVal (:, 1); % DoF values of function.
            end
            assert(ismember(fES.msh.type, "D2T"));
            assert(ismember(fES.elem, ["D2T", "D2LR"]));
            assert(fES.nGlDoF == length(DoFVal));
            fcn = comb(fES.LcBase, sym("coef", [fES.nLcDoF, 1]));
            FEF = FEF@Fcn(fcn.domn, fcn.fun, fcn.coef);
            FEF.msh = fES.msh;
            FEF.elem = fES.elem;
            FEF.ElParm = fES.ElParm;
            FEF.ElCoef = reshape(DoFVal(abs(fES.Lc2Gl)), size(abs(fES.Lc2Gl))) .* sign(fES.Lc2Gl);
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
            % FEF.subt: subtract.
            arguments (Input)
                fcn1 FEF;
                fcn2 FEF;
            end
            arguments (Output)
                fcn FEF;
            end
            assert(isequal(fcn1, fcn2, "domn", "fun", "coef", "elem", "ElParm"));
            fcn = fcn1;
            fcn.ElCoef = fcn1.ElCoef - fcn2.ElCoef;
        end
    end
    % Static functions.
    methods (Static)
        function varargout = multi(FESs, DoFVal)
            % FEF.multi: construct multiple FEF.
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
function flag = isequal(fcn1, fcn2, varargin)
    % isequal: compare FE functions.
    % varargin: names of properties to be compared.
    flag = true;
    for iProp = 1:length(varargin)
        prop = varargin{iProp};
        if ~isequal(fcn1.(prop), fcn2.(prop))
            flag = false;
            return;
        end
    end
end
