function results = runTests(tag, options)
    % runTests: run unit tests of FEM package.

    % Tag        | Content
    % ---------- | --------------------------------------------------
    % Fast       | unit tests and small-mesh regression (default)
    % Slow       | examples and convergence tests
    % Regression | comparison with 2D baseline (small mesh and examples)
    % Style      | coding style checks
    % All        | all tests
    %
    % Examples:
    %   runTests;                                   % all Fast tests
    %   runTests("Fast", "class", ["FETest", "FESTest"]); % Fast tests of the given test classes only
    %   runTests("Slow", "class", "Poisson3DTest");

    arguments (Input)
        tag (1, 1) string {mustBeMember(tag, ["Fast", "Slow", "Regression", "Style", "All"])} = "Fast"; % Tag of tests.
        options.class (1, :) string = string.empty; % Names of test classes to run (all if empty).
    end
    arguments (Output)
        results matlab.unittest.TestResult; % Test results.
    end
    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner;
    import matlab.unittest.Verbosity;
    import matlab.unittest.selectors.HasTag;
    root = fileparts(mfilename("fullpath"));
    addpath(fullfile(root, "FEM_Package"), fullfile(root, "FEM_Test", "unit"), fullfile(root, "FEM_Test", "unit", "helper"), ...
        fullfile(root, "FEM_Test", "baseline"));
    figVis = get(groot, "DefaultFigureVisible");
    set(groot, "DefaultFigureVisible", "off");
    cleanup = onCleanup(@() set(groot, "DefaultFigureVisible", figVis));
    % Cached finite elements and function handles are rebuilt in each run (the package may have changed).
    clear stdFE mapVal
    suite = TestSuite.fromFolder(fullfile(root, "FEM_Test", "unit"));
    if tag ~= "All"
        suite = suite.selectIf(HasTag(tag));
    end
    if ~isempty(options.class)
        suite = suite(ismember(extractBefore(string({suite.Name}), "/"), options.class));
    end
    runner = TestRunner.withTextOutput("OutputDetail", Verbosity.Concise);
    results = runner.run(suite);
    % Summary.
    fprintf("\n==== runTests(""%s"") ====\n", tag);
    if ~isempty(options.class)
        fprintf("Classes: %s\n", strjoin(options.class, ", "));
    end
    fprintf("%s: %s\n", "MATLAB", version);
    fprintf("%d passed, %d failed, %d incomplete, %.1f s\n", ...
        nnz([results.Passed]), nnz([results.Failed]), nnz([results.Incomplete]), sum([results.Duration]));
    for iRes = 1:length(results)
        if results(iRes).Failed
            fprintf("FAILED: %s\n", results(iRes).Name);
        elseif results(iRes).Incomplete
            fprintf("INCOMPLETE: %s\n", results(iRes).Name);
        end
    end
    % Duration of each test class.
    if ~isempty(results)
        [clsName, ~, clsIdx] = unique(extractBefore(string({results.Name}), "/"));
        clsDur = accumarray(clsIdx(:), [results.Duration]');
        for iCls = 1:length(clsName)
            fprintf("  %-20s %8.1f s\n", clsName(iCls), clsDur(iCls));
        end
    end
end
