% Build in ultimatum

function [xCur, xCurCost, info, options] = bfgs_Smooth_release_version(problem, xCur, options)
% Riemannian BFGS solver for smooth objective function.
%
% function [x, cost, info, options] = bfgsSmooth(problem)
% function [x, cost, info, options] = bfgsSmooth(problem, x0)
% function [x, cost, info, options] = bfgsSmooth(problem, x0, options)
% function [x, cost, info, options] = bfgsSmooth(problem, [], options)
%
%
% This is Riemannian BFGS solver (quasi-Newton method), which aims to
% minimize the cost function in problem structure problem.cost. It needs 
% gradient of the cost function.
%
%
% For a description of the algorithm and theorems offering convergence
% guarantees, see the references below.
%
% The initial iterate is xCur if it is provided. Otherwise, a random point on
% the manifold is picked. To specify options whilst not specifying an
% initial iterate, give xCur as [] (the empty matrix).
%
% The two outputs 'xCur' and 'xCurcost' are the last reached point on the manifold
% and its cost. Notice that x is not necessarily the best reached point,
% because this solver is not forced to be a descent method. In particular,
% very close to convergence, it is sometimes preferable to accept very
% slight increases in the cost value (on the order of the machine epsilon)
% in the process of reaching fine convergence. In practice, this is not a
% limiting factor, as normally one does not need fine enough convergence
% that this becomes an issue.
% 
% The output 'info' is a struct-array which contains information about the
% iterations:
%   iter (integer)
%       The (outer) iteration number, or number of steps considered
%       (whether accepted or rejected). The initial guess is 0.
%	cost (double)
%       The corresponding cost value.
%	gradnorm (double)
%       The (Riemannian) norm of the gradient.
%	time (double)
%       The total elapsed time in seconds to reach the corresponding cost.
%	stepsize (double)
%       The size of the step from the previous to the new iterate.
%   accepted (boolean)
%       1 if the current step is accepted in the cautious update. 0 otherwise
%   And possibly additional information logged by options.statsfun.
% For example, type [info.gradnorm] to obtain a vector of the successive
% gradient norms reached at each (outer) iteration.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%
%   tolgradnorm (1e-6)
%       The algorithm terminates if the norm of the gradient drops below
%       this. For well-scaled problems, a rule of thumb is that you can
%       expect to reduce the gradient norm by 8 orders of magnitude
%       (sqrt(eps)) compared to the gradient norm at a "typical" point (a
%       rough initial iterate for example). Further decrease is sometimes
%       possible, but inexact floating point arithmetic will eventually
%       limit the final accuracy. If tolgradnorm is set too low, the
%       algorithm may end up iterating forever (or at least until another
%       stopping criterion triggers).
%   maxiter (1000)
%       The algorithm terminates if maxiter (outer) iterations were executed.
%   maxtime (Inf)
%       The algorithm terminates if maxtime seconds elapsed.
%	miniter (3)
%       Minimum number of outer iterations (used only if useRand is true).
%   minstepsize (1e-10)
%     The minimum norm of the tangent vector that points from the current
%     point to the next point. If the norm is less than minstepsize, the 
%     program will terminate.
%   memory(30)
%     The number of previous iterations the program remembers in LBFGS. This is used 
%     to approximate the Hessian at the current point. Because of difficulty
%     of maintaining a representation of hessian in terms of coordinates, and
%     thus a recursive computation for the direction pointing to the next
%     point is done by considering approximating Hessian as an operator that takes
%     a vector and outputs a vector in the tangent space. Theoretically, a
%     vector recurse back memory size number of times and thus memory size 
%     is linear with the time taken to compute directions towards the next
%     point.
%     It can take any value >= 0, or Inf (which will then take value maxiter).
%   strict_inc_func(@(x) x)
%     The Cautious step needs a real function that has value 0 at x = 0, and 
%     strictly increasing. See details in Wen Huang's paper
%     "A Riemannian BFGS Method without Differentiated Retraction for 
%     Nonconvex Optimization Problems
%   statsfun (none)
%       Function handle to a function that will be called after each
%       iteration to provide the opportunity to log additional statistics.
%       They will be returned in the info struct. See the generic Manopt
%       documentation about solvers for further information. statsfun is
%       called with the point x that was reached last, after the
%       accept/reject decision. See comment below.
%   stopfun (none)
%       Function handle to a function that will be called at each iteration
%       to provide the opportunity to specify additional stopping criteria.
%       See the generic Manopt documentation about solvers for further
%       information.
%   verbosity (2)
%       Integer number used to tune the amount of output the algorithm
%       generates during execution (mostly as text in the command window).
%       The higher, the more output. 0 means silent. 3 and above includes a
%       display of the options structure at the beginning of the execution.
%   debug (false)
%       Set to true to allow the algorithm to perform additional
%       computations for debugging purposes. If a debugging test fails, you
%       will be informed of it, usually via the command window. Be aware
%       that these additional computations appear in the algorithm timings
%       too, and may interfere with operations such as counting the number
%       of cost evaluations, etc. (the debug calls get storedb too).
%   storedepth (20)
%       Maximum number of different points x of the manifold for which a
%       store structure will be kept in memory in the storedb. If the
%       caching features of Manopt are not used, this is irrelevant. If
%       memory usage is an issue, you may try to lower this number.
%       Profiling may then help to investigate if a performance hit was
%       incured as a result.
%
% Notice that statsfun is called with the point x that was reached last,
% after the accept/reject decision. Hence: if the step was accepted, we get
% that new x, with a store which only saw the call for the cost and for the
% gradient. If the step was rejected, we get the same x as previously, with
% the store structure containing everything that was computed at that point
% (possibly including previous rejects at that same point). Hence, statsfun
% should not be used in conjunction with the store to count operations for
% example. Instead, you should use storedb's shared memory for such
% purposes (either via storedb.shared, or via store.shared, see
% online documentation). It is however possible to use statsfun with the
% store to compute, for example, other merit functions on the point x
% (other than the actual cost function, that is).
%
%
% Please cite the Manopt paper as well as the research paper:
%     @TECHREPORT{HAG2017,
%     author = "Wen Huang and P.-A. Absil and K. A. Gallivan",
%     title = "A Riemannian BFGS Method without Differentiated Retraction for Nonconvex Optimization Problems",
%     institution = "U.C.Louvain",
%     number = "UCL-INMA-2017.04",
%     year = 2017,
%     }
%


% This file is part of Manopt: www.manopt.org.
% Change log: 
%
%   CL July 15, 2017:
%        Finished the first released version

    % Verify that the problem description is sufficient for the solver.
    if ~canGetCost(problem)
        warning('manopt:getCost', ...
            'No cost provided. The algorithm will likely abort.');
    end
    if ~canGetGradient(problem) && ~canGetApproxGradient(problem)
        % Note: we do not give a warning if an approximate gradient is
        % explicitly given in the problem description, as in that case the user
        % seems to be aware of the issue.
        warning('manopt:getGradient:approx', ...
            ['No gradient provided. The correctness is not garanteed.\n'....
              'The program is not designed for approximating gradient.\n' ...
            'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
        problem.approxgrad = approxgradientFD(problem);
    end
    
    % Local defaults for the program
    localdefaults.minstepsize = 1e-10;
    localdefaults.maxiter = 1000;
    localdefaults.tolgradnorm = 1e-6;
    localdefaults.memory = 30;
    localdefaults.strict_inc_func = @(x) x;
    localdefaults.max_iter_line_search = 25;
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % To make sure memory in range [0, Inf)
    if options.memory < 0
        options.memory = 0;
    elseif options.memory == Inf
        options.memory = options.maxiter;
    end
    
    M = problem.M;
    
    % Create a random starting point if no starting point
    % is provided.
    if ~exist('xCur','var')|| isempty(xCur)
        xCur = M.rand(); 
    end
    
    timetic = tic();
    
    % Create a store database and get a key for the current x
    storedb = StoreDB(options.storedepth);
    key = storedb.getNewKey();
    
    % __________Initialization of variables______________
    % number of current element in memory
    k = 0;  
    %number of total iteration in BFGS
    iter = 0; 
    % saves vector that represents x_{k}'s projection on 
    % x_{k+1}'s tangent space. And transport it to most
    % current point's tangent space after every iteration.
    sHistory = cell(1, options.memory);
    % saves gradient of x_{k} by transporting it to  
    % x_{k+1}'s tangent space. And transport it to most
    % current point's tangent space after every iteration.
    yHistory = cell(1, options.memory);
    % saves inner(sk,yk)
    rhoHistory = cell(1, options.memory);
    % scaling of direction given by getDirection for acceptable step
    alpha = 1; 
    % scaling of initial matrix, BB.
    scaleFactor = 1;
    % norm of the step
    stepsize = 1;
    accepted = 1;
    xCurGradient = getGradient(problem, xCur, storedb, key);
    xCurGradNorm = M.norm(xCur, xCurGradient);
    xCurCost = getCost(problem, xCur);
    lsstats = [];
    ultimatum = 0;
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    if options.verbosity >= 2
    fprintf(' iter\t               cost val\t                 grad. norm\t        alpha \n');
    end
    
    while (1)
%------------------------ROUTINE----------------------------

        % Display iteration information
        if options.verbosity >= 2
        %_______Print Information and stop information________
        fprintf('%5d\t%+.16e\t%.8e\t %.4e\n', iter, xCurCost, xCurGradNorm, alpha);
        end
        
        % Start timing this iteration
        timetic = tic();
        
        % Run standard stopping criterion checks
        [stop, reason] = stoppingcriterion(problem, xCur, options, ...
            info, iter+1);
        
        % If none triggered, run specific stopping criterion check
        if ~stop 
            if stats.stepsize < options.minstepsize
                if ultimatum == 0
                    fprintf(['stepsize is too small, restart the bfgs procedure' ...
                        'with the current point\n']);
                    k = 0;
                    ultimatum = 1;
                else
                    stop = true;
                    reason = sprintf(['Last stepsize smaller than minimum '  ...
                        'allowed; options.minstepsize = %g.'], ...
                        options.minstepsize);
                end
            else
                ultimatum = 0;
            end
        end  
        
        if stop
            if options.verbosity >= 1
                fprintf([reason '\n']);
            end
            break;
        end

        %_______Get Direction___________________________

        p = getDirection(M, xCur, xCurGradient, sHistory,...
            yHistory, rhoHistory, scaleFactor, min(k, options.memory));

        %_______Line Search____________________________
        [alpha, xNext, xNextCost, lsstats] = linesearchArmijo_start_with_alpha_eq_one(problem,...
            xCur, p, xCurCost, M.inner(xCur,xCurGradient,p), options.max_iter_line_search); 
        step = M.lincomb(xCur, alpha, p);
        stepsize = M.norm(xCur, p)*alpha;

        
        %_______Updating the next iteration_______________
        newkey = storedb.getNewKey();
        xNextGradient = getGradient(problem, xNext, storedb, newkey);
        sk = M.transp(xCur, xNext, step);
        yk = M.lincomb(xNext, 1, xNextGradient,...
            -1, M.transp(xCur, xNext, xCurGradient));

        inner_sk_yk = M.inner(xNext, yk, sk);
        % If cautious step is not accepted, then we do no take the
        % current sk, yk into account. Otherwise, we record it 
        % and use it in approximating hessian.
        % sk, yk are maintained in the most recent point's 
        % tangent space by transport.
        if (inner_sk_yk / M.inner(xNext, sk, sk))>= options.strict_inc_func(xCurGradNorm)
            accepted = 1;
            rhok = 1/inner_sk_yk;
            scaleFactor = inner_sk_yk / M.inner(xNext, yk, yk);
            if (k>= options.memory)
                % sk and yk are saved from 1 to the end
                % with the most currently recorded to the 
                % rightmost hand side of the cells that are
                % occupied. When memory is full, do a shift
                % so that the rightmost is earliest and replace
                % it with the most recent sk, yk.
                for  i = 2:options.memory
                    sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                    yHistory{i} = M.transp(xCur, xNext, yHistory{i});
                end
                if options.memory > 1
                sHistory = sHistory([2:end 1]);
                yHistory = yHistory([2:end 1]);
                rhoHistory = rhoHistory([2:end 1]);
                end
                if options.memory > 0
                    sHistory{options.memory} = sk;
                    yHistory{options.memory} = yk;
                    rhoHistory{options.memory} = rhok;
                end
            else
                for  i = 1:k
                    sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                    yHistory{i} = M.transp(xCur, xNext, yHistory{i});
                end
                sHistory{k+1} = sk;
                yHistory{k+1} = yk;
                rhoHistory{k+1} = rhok;
            end
            k = k+1;
        else
            accepted = 0;
            for  i = 1:min(k,options.memory)
                sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                yHistory{i} = M.transp(xCur, xNext, yHistory{i});
            end
        end
        iter = iter + 1;
        xCur = xNext;
        xCurGradient = xNextGradient;
        xCurGradNorm = M.norm(xCur, xNextGradient);
        xCurCost = xNextCost;
        
        % Make sure we don't use too much memory for the store database
        storedb.purge();
        
        key = newkey;
        
        % Log statistics for freshly executed iteration
        stats = savestats();
        info(iter+1) = stats; 
        
    end

    
    info = info(1:iter+1);

    if options.verbosity >= 1
        fprintf('Total time is %f [s] (excludes statsfun)\n', ...
                info(end).time);
    end

    % Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.iter = iter;
        stats.cost = xCurCost;
        stats.gradnorm = xCurGradNorm;
        if iter == 0
            stats.stepsize = NaN;
            stats.accepted = NaN;
            stats.time = toc(timetic);
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
            stats.accepted = accepted;
        end
        stats.linesearch = lsstats;
        stats = applyStatsfun(problem, xCur, storedb, key, options, stats);
    end

end

% BFGS step, look at Wen's paper for details, it is equivalent
% to transporting all the way back to the tangent spaces
% of previous points and to the operation if isotransp is
% enforced, but in practice, there is no observed difference
% in them, if your problem requires isotransp, it may be good
% to replace transp with isotransp. There are built in isotransp
% for sphere and Obliquefactory
function dir = getDirection(M, xCur, xCurGradient, sHistory, yHistory, rhoHistory, scaleFactor, k)
    q = xCurGradient;
    inner_s_q = cell(1, k);
    for i = k : -1: 1
        inner_s_q{i} = rhoHistory{i}*M.inner(xCur, sHistory{i},q);
        q = M.lincomb(xCur, 1, q, -inner_s_q{i}, yHistory{i});
    end
    r = M.lincomb(xCur, scaleFactor, q);
    for i = 1: k
         omega = rhoHistory{i}*M.inner(xCur, yHistory{i},r);
         r = M.lincomb(xCur, 1, r, inner_s_q{i}-omega, sHistory{i});
    end
    dir = M.lincomb(xCur, -1, r);
end


function [alpha, xNext, xNextCost, lsstats] = ...
                  linesearchArmijo_start_with_alpha_eq_one(problem, x, d, f0, df0, max_iter_line_search)

    % Backtracking default parameters. These can be overwritten in the
    % options structure which is passed to the solver.
    contraction_factor = .5;
    suff_decr = 1e-4;
    
    % At first, we have no idea of what the step size should be.
    alpha = 1;

    % Make the chosen step and compute the cost there.
    xNext = problem.M.retr(x, d, alpha);
    xNextCost = getCost(problem, xNext);
    num_cost_eval = 1;
    
    % Backtrack while the Armijo criterion is not satisfied
    while xNextCost > f0 + suff_decr*alpha*df0
        
        % Reduce the step size,
        alpha = contraction_factor * alpha;
        
        % and look closer down the line
        xNext = problem.M.retr(x, d, alpha);
        xNextCost = getCost(problem, xNext);
        num_cost_eval = num_cost_eval + 1;
        
        % Make sure we don't run out of budget
        if num_cost_eval >= max_iter_line_search
            break;
        end
        
    end
    
    lsstats.num_cost_eval = num_cost_eval;
    % If we got here without obtaining a decrease, we reject the step.
    if xNextCost > f0
        alpha = 0;
        xNext = x;
        xNextCost = f0; 
    end
end
