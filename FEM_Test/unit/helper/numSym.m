function val = numSym(expr, domn, parm)
    % numSym: evaluate symbolic expression numerically by substituting parameter of domain.
    arguments (Input)
        expr; % Symbolic expression.
        domn (1, 1) string; % Domain whose parameter is substituted, e.g. "D2T".
        parm (:, :) double; % Value of parameter.
    end
    arguments (Output)
        val double;
    end
    symParm = MshEnt(domn).parm;
    assert(isequal(size(symParm), size(parm)));
    val = double(subs(sym(expr), symParm, parm));
end
