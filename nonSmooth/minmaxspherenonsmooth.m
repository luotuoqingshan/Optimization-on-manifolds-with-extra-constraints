function minmaxspherenonsmooth

    d = 10;
    n = 20;
    % Create the problem structure.
    manifold = spherefactory(d);
    problem.M = manifold;
    M = problem.M;
    data = zeros(d,n);
    for pointnum = 1 : n
        data(:,pointnum) = problem.M.rand();  %not applicable to other manifolds
        data(1,:) = abs(data(1,:));
    end
    
    cost = @(X) costFun(data, X);
    grad = @(X) gradFun(data, X);
    
    
    % Define the problem cost function and its Euclidean gradient.
    problem.cost  = cost;
    problem.egrad = grad;
    problem.reallygrad = grad;
    
    %Set options
    options.linesearchVersion = 4;
    options.memory = 400;

    X1 = problem.M.rand();
    X2 = problem.M.rand();
    X3 = problem.M.rand();
    options.assumedoptX = problem.M.rand();
    while (1)
        xCur = problem.M.rand();
        [gradnorms, alphas, stepsizes, costs, distToAssumedOptX, xHistory, X1, time] = bfgsnonsmooth(problem, xCur, options);
        if M.dist(X1,X2)+M.dist(X1,X3)+M.dist(X2,X3) <= 1e-5
            break;
        else
            X3 = X2;
            X2 = X1;
        end
    end
    profile clear;
    profile on;
    options.assumedoptX = X1;
    xCur = problem.M.rand();
    xCur = [-1; 0; 0];
    [gradnorms, alphas, stepsizes, costs, distToAssumedOptX, xHistory, xCur, time] = bfgsnonsmooth(problem, xCur, options);
    
    profile off;
    profile report
    
    disp(xCur)
    figure;
    
    subplot(2,2,1)
    semilogy(gradnorms, '.-');
    xlabel('Iter');
    ylabel('GradNorms');

    titletest = sprintf('Time: %f', time);
    title(titletest);
    
    subplot(2,2,2)
    plot(alphas, '.-');
    xlabel('Iter');
    ylabel('Alphas');

    subplot(2,2,3)
    semilogy(stepsizes, '.-');
    xlabel('Iter');
    ylabel('stepsizes');

    subplot(2,2,4)
    semilogy(distToAssumedOptX, '.-');
    xlabel('Iter');
    ylabel('costs');

    
    figure
    surfprofile(problem, xCur);
    
    if d == 3
        figure;
        % Plot the sphere
        [sphere_x, sphere_y, sphere_z] = sphere(50);
        handle = surf(sphere_x, sphere_y, sphere_z);
        set(handle, 'FaceColor', [152,186,220]/255);
        set(handle, 'FaceAlpha', .5);
        set(handle, 'EdgeColor', [152,186,220]/255);
        set(handle, 'EdgeAlpha', .5);
        daspect([1 1 1]);
        box off;
        axis off;
        hold on;
        % Add the chosen points
        Y = cell2mat(xHistory);
        Y = 1.02*Y;
        [row, col] = size(Y);
        plot3(Y(1,:), Y(2,:), Y(3,:), 'r.', 'MarkerSize', 5);
        plot3(data(1,:), data(2,:), data(3,:), 'g.', 'MarkerSize', 20);
        % And connect the points which are at minimal distance,
%         % within some tolerance.
        for k = 1 : col-1
            i = k; j = k+1;
            plot3(Y(1, [i j]), Y(2, [i j]), Y(3, [i j]), 'k-');
        end
        hold off;
    end
    
    
    
        function val = costFun(data, x)
            Inner = - x.'*data;
            val = max(Inner(:));
        end

        function val = gradFun(data, x)
            Inner = - x.'*data;
            [maxval,pos] = max(Inner(:));
            val = - data(:, pos);
        end
end