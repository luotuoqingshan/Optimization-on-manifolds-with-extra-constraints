% TO DO List: approximation of gradient?
%  To write a manual
% Output info structs (Add acceptness)
% options adda function handle for  >theta(x) accept the step inequality
% To fix memory == 0, 1, Inf case
% Delete other linesearch versions.
% Initialstepsize is not 1 actually.
% lsstats
% maxsteps in linesearch is supposed to be in options.

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
#If no gradient is provided, an approximation of the gradient is computed,
# but this can be slow for manifolds of high dimension.
%
% For a description of the algorithm and theorems offering convergence
% guarantees, see the references below. Documentation for this solver is
% available online at:
%
# http://www.manopt.org/solver_documentation_trustregions.html
%
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
# The output 'info' is a struct-array which contains information about the
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
%     It can take any value >= 0, or Inf (which is effectively RBFGS).
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
%     @Article{genrtr,
%       Title    = {Trust-region methods on {Riemannian} manifolds},
%       Author   = {Absil, P.-A. and Baker, C. G. and Gallivan, K. A.},
%       Journal  = {Foundations of Computational Mathematics},
%       Year     = {2007},
%       Number   = {3},
%       Pages    = {303--330},
%       Volume   = {7},
%       Doi      = {10.1007/s10208-005-0179-9}
%     }
%
% See also: steepestdescent conjugategradient manopt/examples

% An explicit, general listing of this algorithm, with preconditioning,
% can be found in the following paper:
%     @Article{boumal2015lowrank,
%       Title   = {Low-rank matrix completion via preconditioned optimization on the {G}rassmann manifold},
%       Author  = {Boumal, N. and Absil, P.-A.},
%       Journal = {Linear Algebra and its Applications},
%       Year    = {2015},
%       Pages   = {200--239},
%       Volume  = {475},
%       Doi     = {10.1016/j.laa.2015.02.027},
%     }

% When the Hessian is not specified, it is approximated with
% finite-differences of the gradient. The resulting method is called
% RTR-FD. Some convergence theory for it is available in this paper:
% @incollection{boumal2015rtrfd
% 	author={Boumal, N.},
% 	title={Riemannian trust regions with finite-difference Hessian approximations are globally convergent},
% 	year={2015},
% 	booktitle={Geometric Science of Information}
% }


% This file is part of Manopt: www.manopt.org.
% This code is an adaptation to Manopt of the original GenRTR code:
% RTR - Riemannian Trust-Region
% (c) 2004-2007, P.-A. Absil, C. G. Baker, K. A. Gallivan
% Florida State University
% School of Computational Science
% (http://www.math.fsu.edu/~cbaker/GenRTR/?page=download)
% See accompanying license file.
% The adaptation was executed by Nicolas Boumal.
%
%
% Change log: 
%
%   NB April 3, 2013:
%       




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
            ['No gradient provided. Using an FD approximation instead (slow).\n' ...
            'It may be necessary to increase options.tolgradnorm.\n' ...
            'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
        problem.approxgrad = approxgradientFD(problem);
    end

    
    % Set local defaults here
    %linesearchVersion: 0 - Basic Armijo starting with step size 1.
    %                         1 - Manopt Armijo with guessing step size
    %                         2 - Qi's linesearch
    %                         3 - alpha fix to be 1.
    
    localdefaults.minstepsize = 1e-10;
    localdefaults.maxiter = 1000;
    localdefaults.tolgradnorm = 1e-6;
    localdefaults.memory = 30;
    localdefaults.linesearchVersion = 0;
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    timetic = tic();
    
    
    % Create a store database and get a key for the current x
    storedb = StoreDB(options.storedepth);
    key = storedb.getNewKey();
    
    M = problem.M;

    if ~exist('xCur','var')|| isempty(xCur)
        xCur = xCur;
    else
        xCur = M.rand();
    end

    
    k = 0;
    iter = 0;
    sHistory = cell(1, options.memory);
    yHistory = cell(1, options.memory);
    rhoHistory = cell(1, options.memory);
    alpha = 1;
    scaleFactor = 1;
    stepsize = 1;
    xCurGradient = getGradient(problem, xCur);
    xCurGradNorm = M.norm(xCur, xCurGradient);
    xCurCost = getCost(problem, xCur);
    
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
        if ~stop && stats.stepsize < options.minstepsize
            stop = true;
            reason = sprintf(['Last stepsize smaller than minimum '  ...
                'allowed; options.minstepsize = %g.'], ...
                options.minstepsize);
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
        switch options.linesearchVersion
            case 0
                [alpha, xNext, xNextCost] = linesearchArmijo_start_with_alpha_eq_one(problem,...
                    xCur, p, xCurCost, M.inner(xCur,xCurGradient,p)); %Check if df0 is right
                step = M.lincomb(xCur, alpha, p);
                stepsize = M.norm(xCur, p)*alpha;
            case 1
                [stepsize, xNext, newkey, lsstats] =linesearch(problem, xCur, p, xCurCost, M.inner(xCur,xCurGradient,p));
                alpha = stepsize/M.norm(xCur, p);
                step = M.lincomb(xCur, alpha, p);
                xNextCost = getCost(problem, xNext);
            case 2
                [xNextCost,alpha] = linesearch_qi(problem, M, xCur, p, M.inner(xCur,xCurGradient,p), alpha);
                step = M.lincomb(xCur, alpha, p);
                stepsize = M.norm(xCur, step);
                xNext = M.exp(xCur, step, 1);
            otherwise
                alpha = 1;
                step = M.lincomb(xCur, alpha, p);
                stepsize = M.norm(xCur, step);
                xNext = M.retr(xCur, step, 1);
                xNextCost = getCost(problem, xNext);
        end
        
        %_______Updating the next iteration_______________
        newkey = storedb.getNewKey();
        xNextGradient = getGradient(problem, xNext);
        sk = M.transp(xCur, xNext, step);
        yk = M.lincomb(xNext, 1, xNextGradient,...
            -1, M.transp(xCur, xNext, xCurGradient));

        inner_sk_yk = M.inner(xNext, yk, sk);
        if (inner_sk_yk / M.inner(xNext, sk, sk))>= xCurGradNorm
            rhok = 1/inner_sk_yk;
            scaleFactor = inner_sk_yk / M.inner(xNext, yk, yk);
            if (k>= options.memory)
                for  i = 2:options.memory
                    sHistory{i} = M.transp(xCur, xNext, sHistory{i});
                    yHistory{i} = M.transp(xCur, xNext, yHistory{i});
                end
                sHistory = sHistory([2:end 1]);
                sHistory{options.memory} = sk;
                yHistory = yHistory([2:end 1]);
                yHistory{options.memory} = yk;
                rhoHistory = rhoHistory([2:end 1]);
                rhoHistory{options.memory} = rhok;
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
            stats.time = toc(timetic);
        else
            stats.stepsize = stepsize;
            stats.time = info(iter).time + toc(timetic);
        end
        stats.linesearch = [];
        stats = applyStatsfun(problem, xCur, storedb, key, options, stats);
    end

end

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


function [alpha, xNext, xNextCost] = ...
                  linesearchArmijo_start_with_alpha_eq_one(problem, x, d, f0, df0)

    % Backtracking default parameters. These can be overwritten in the
    % options structure which is passed to the solver.
    contraction_factor = .5;
    suff_decr = 1e-4;
    max_steps = 25;
    
    % At first, we have no idea of what the step size should be.
    alpha = 1;

    % Make the chosen step and compute the cost there.
    xNext = problem.M.retr(x, d, alpha);
    xNextCost = getCost(problem, xNext);
    cost_evaluations = 1;
    
    % Backtrack while the Armijo criterion is not satisfied
    while xNextCost > f0 + suff_decr*alpha*df0
        
        % Reduce the step size,
        alpha = contraction_factor * alpha;
        
        % and look closer down the line
        xNext = problem.M.retr(x, d, alpha);
        xNextCost = getCost(problem, xNext);
        cost_evaluations = cost_evaluations + 1;
        
        % Make sure we don't run out of budget
        if cost_evaluations >= max_steps
            break;
        end
        
    end
    
    % If we got here without obtaining a decrease, we reject the step.
    if xNextCost > f0
        alpha = 0;
        xNext = x;
        xNextCost = f0; 
    end
    
%     fprintf('alpha = %.16e\n', alpha)
end


function [costNext,alpha] = linesearch_qi(problem, M, x, p, df0, alphaprev)

    alpha = alphaprev;
    costAtx = getCost(problem,x);
    while (getCost(problem,M.exp(x,p,2*alpha))-costAtx < alpha*df0)
        alpha = 2*alpha;
    end
    costNext = getCost(problem,M.exp(x,p,alpha));
    diff = costNext - costAtx;
    while (diff>= 0.5*alpha*df0)
        if (diff == 0)
            alpha = 0;
            break;
        end
        alpha = 0.5 * alpha;
        costNext = getCost(problem,M.exp(x,p,alpha));
        diff = costNext - costAtx;
    end
%     fprintf('alpha = %.16e\n',alpha);    
end