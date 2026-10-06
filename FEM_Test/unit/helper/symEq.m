function flag = symEq(act, ref)
    % symEq: check whether two symbolic expressions are identical.
    arguments (Input)
        act; % Actual expression (sym or double).
        ref; % Reference expression (sym or double).
    end
    arguments (Output)
        flag (1, 1) logical;
    end
    act = sym(act);
    ref = sym(ref);
    if ~isequal(size(act), size(ref))
        flag = false;
        return;
    end
    flag = all(isAlways(simplify(act - ref) == 0, "Unknown", "false"), "all");
end
