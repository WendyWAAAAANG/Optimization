function [x, k, y_opt] = newton_backtracking(f, x0, tol)
    max_iter = 500;
    x = x0;
    % counter -- time of iteration.
    k = 0;
    
     % find gradient and hessian matrix.
    syms x1 x2
    vars = [x1; x2];
    % calculate gradient and hessian function of f.
    grad = matlabFunction(gradient(f(vars), vars));
    hess = matlabFunction(hessian(f(vars), vars));
    
    % iteration.
    while norm(grad(x(1), x(2))) > tol && k < max_iter
        alpha = 1;
        % for this iteration, calculate gradient and hessian value.
        g = grad(x(1), x(2));
        h = hess(x(1), x(2));
        % calculate update term.
        %s = -h \ g;
        s = -inv(h) * g;
        x_temp = x + s*alpha;
        if f(x_temp) > f(x)
            alpha = alpha / 2;
        end
        x = x + s*alpha;
        k = k + 1;
    end
    y_opt = f(x);
end