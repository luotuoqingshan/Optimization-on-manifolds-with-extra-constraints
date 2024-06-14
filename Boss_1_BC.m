function Boss_1_BC(varargin) 
    format long g;
    p = inputParser;
    addOptional(p, 'graph', 'G1', @ischar);
    addOptional(p, 'seed', 0, @isnumeric);
    % 2: RALM 3: Q_LSE 4: LSE
    addOptional(p, 'solver', 'RALM', @ischar);
    parse(p, varargin{:});

    seed = p.Results.seed;
    minbisec_data = p.Results.graph;
    solver = p.Results.solver;
    %--------------------------Balanced Cut------------------------------------
    %specifier.matlabversion = 0; %0 if older than 2015 1 otherwise

    rank = 2;     %Graph Bisection
    rng(seed, 'twister');

    fprintf("Solving Minimum Bisection SDP for %s\n", minbisec_data);

    %% modify the path before running
    data = load(['data/MinimumBisection/', minbisec_data, '.mat']);
    A = data.A;
    n = size(A, 1);
    L = spdiags(A*ones(n,1),0,n,n) - A;

    %________Experiment_____
    options.maxOuterIter = 100000000;
    options.maxtime = 3600;
    options.minstepsize = 1e-10;


    solver_str = convertCharsToStrings(solver);
    if solver_str == "RALM"
        specifier = 2;
    elseif solver_str == "Q_LQH"
        specifier = 3;
    else % "Q_LSE"
        specifier = 4;
    end
    result = clientconstraint_oblique_balancedcut(L, rank, options, specifier);
    disp(result);
    if ~exist(['output/MinimumBisection/', minbisec_data, '/', solver],'dir') 
        mkdir(['output/MinimumBisection/', minbisec_data, '/', solver]); 
    end
    save(['output/MinimumBisection/', minbisec_data,...
         '/', solver, '/', solver, '-seed-', num2str(seed), '.mat'],'result','-v7.3');
end