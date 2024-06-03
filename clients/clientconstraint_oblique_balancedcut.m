function data = clientconstraint_oblique_balancedcut(L, rankY, methodoptions, specifier)

% data stores the output
% data(1) stores the entry of primal violation with maximum magnitude 
% data(2) stores objective
% data(3) stores time
% data(4) stores the value of the cut
% data(5) stores relative primal violation
data = NaN(5, 1);

[N, ~] = size(L);
manifold = obliquefactory(rankY, N);
problem.M = manifold;
problem.cost = @(u) costFun(u);
problem.egrad = @(u) gradFun(u);
x0 = problem.M.rand();

%     DEBUG only
%     checkgradient(problem);

%-------------------------Set-up Constraints-----------------------
colones = ones(N, 1);

eq_constraints_cost = cell(1,1);
eq_constraints_cost{1} = @(U) norm(U*colones,2)^2;

eq_constraints_grad = cell(1,1);
eq_constraints_grad{1} = @(U) meanzero_eq_constraint(U);

problem.eq_constraint_cost = eq_constraints_cost;
problem.eq_constraint_grad = eq_constraints_grad;

%     Debug Only
%     checkconstraints(problem)

condet = constraintsdetail(problem);

%     ------------------------- Solving ---------------------------
    options = methodoptions;
    
    if specifier == 1
        %MINI-SUM-MAX
        fprintf('Starting Mini-sum-max \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaMinimax(problem, x0, options);
        time = toc(timetic);

        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1) = maxviolation;
        % original minimum bisection SDP
        data(2) = 0.25 * cost; 
        data(3) = time;
        data(4) = random_rounding(L, xfinal);
        data(5) = rel_primal_vio(L, xfinal);
    end
    
    if specifier == 2
        %ALM
        fprintf('Starting ALM \n');
        timetic = tic();
        [xfinal, info] = almbddmultiplier(problem, x0, options);
        time = toc(timetic);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1) = maxviolation;
        % original minimum bisection SDP
        data(2) = 0.25*cost;
        data(3) = time;
        data(4) = random_rounding(L, xfinal);
        data(5) = rel_primal_vio(L, xfinal);
    end

    if specifier == 3
        %LQH
        fprintf('Starting LQH \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaSmoothinglqh(problem, x0, options);
        time = toc(timetic);

        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1) = maxviolation;
        data(2) = 0.25*cost;
        data(3) = time;
        data(4) = random_rounding(L, xfinal);
        data(5) = rel_primal_vio(L, xfinal);
    end
    
    
    if specifier == 4
        %LSE
        fprintf('Starting LSE \n');
        timetic = tic();
        [xfinal, info] = exactpenaltyViaSmoothinglse(problem, x0, options);
        time = toc(timetic);
        
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        maxviolation = max(maxviolation, manifoldViolation(xfinal));
        data(1) = maxviolation;
        data(2) = 0.25*cost;
        data(3) = time;
        data(4) = random_rounding(L, xfinal);
        data(5) = rel_primal_vio(L, xfinal);
    end
    
    if specifier == 5
        %FMINCON
        fprintf('Starting fmincon \n');
        maxiter = 1000000;
        if specifier.matlabversion == 0
            % Use this if you are at 2015a or older.
            options = optimoptions('fmincon', 'MaxIter', maxiter, 'MaxFunEvals', maxiter,...
                'GradObj', 'on', 'GradConstr', 'on', 'OutputFcn', @outfun,...
                'TolX', methodoptions.minstepsize);
        else
            % Use this otherwise
            options = optimoptions('fmincon', 'MaxIterations', maxiter, 'MaxFunctionEvaluations', maxiter,...
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'OutputFcn', @outfun,...
                'StepTolerance', methodoptions.minstepsize);
        end
        timetic = tic();
        [xfinal, fval, exitflag, output] = fmincon(@(v) costFunfmincon(v), x0(:), [], [], [], [], [], [], @nonlcon, options);
        time = toc(timetic);
        
        xfinal = reshape(xfinal, [rankY, N]);
        [maxviolation, meanviolation, cost] = evaluation(problem, xfinal, condet);
        data(1) = output.constrviolation;
        data(2) = 0.25*cost;
        data(3) = time;
        data(4) = random_rounding(L, xfinal);
        data(5) = rel_primal_vio(L, xfinal);
    end
     
     %------------------------sub functions-----------
    
    function stop = outfun(x, optimValues, state)
        stop = false;
        if toc(timetic) > methodoptions.maxtime
            stop = true;
        end
    end 
     
    function [f, g] = costFunfmincon(v)
        Y = reshape(v, [rankY, N]);
        f = trace((Y) * L * (Y.') );
        if nargout > 1
            g = Y*L + Y*L.';
            g = g(:);
        end
    end 

    function [c, ceq, gradc, gradceq] = nonlcon(v)
        Y = reshape(v, [rankY, N]);
        ceq = zeros(N+1,1);
        for rowCeq = 1: N
            ceq(rowCeq,1) = Y(:,rowCeq).'*Y(:,rowCeq) - 1;
        end
        ceq(N+1,1) = colones.' *(Y.')* Y * colones;
        c = [];
        if nargout > 2
            gradc = [];
            gradceq = zeros(rankY*N, N+1);
            grad = 2* repmat(Y*colones, 1, N);
            gradceq(:, N+1) = grad(:);
            for rowCeq = 1: N
                grad = zeros(rankY, N);
                grad(:, rowCeq) = 2*Y(:, rowCeq);
                gradceq(:, rowCeq) = grad(:);
            end
        end
    end

    function val = meanzero_eq_constraint(U)
        val = 2* repmat(U*colones, 1, N);
    end

    function val = costFun(u)
        val = trace((u) * L * (u.') );
    end

    function val = gradFun(u)
        val = u*L + u*L.';
    end

    function manvio = manifoldViolation(x)
        %Oblique Factory:
        manvio = max(abs(diag(x.'*x)-colones));
    end

    function cutvalue = random_rounding(L, x)
        cutvalue = Inf;
        n = size(L, 1);
        for repeat = 1:100
            y = x'*randn(size(x, 1), 1);
            [~, ord] = sort(y);
            partition_vec = zeros(n, 1);
            for i = 1:n
                if i <= n/2
                    partition_vec(ord(i)) = 1;
                else
                    partition_vec(ord(i)) = -1;
                end
            end
            rankvalue = (0.25*partition_vec'*(L*partition_vec));
            cutvalue = min(cutvalue, rankvalue);
        end
    end

    function rel_pvio = rel_primal_vio(L, x)
        n = size(L, 1);
        allones = ones(n, 1);
        pvio = ones(n+1, 1);
        pvio(1:n) = diag(x.'*x) - allones;
        pvio(n+1) = norm(x*allones, 'fro')^2;
        rel_pvio = norm(pvio, 'fro') / norm(allones, 'fro'); 
    end
end

