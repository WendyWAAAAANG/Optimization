function [x, k, y_opt] = BFGS_newton(f, x0, B0, tol)
    max_iter = 500;
    x = x0;
    % counter -- time of iteration.
    k = 0;
    H = B0;
    
    % find gradient and hessian matrix.
    syms x1 x2
    vars = [x1; x2];
    % calculate gradient and hessian function of f.
    grad = matlabFunction(gradient(f(vars), vars));
    
    % iteration.
    while norm(grad(x(1), x(2))) > tol && k < max_iter
        % alpha initialization.
        alpha = 1;
        % for this iteration, calculate gradient and hessian value.
        g = grad(x(1), x(2));
        
        % calculate the update term.
        p = -H \ g;
        
        % using line search to find suitable learning rate.
        s_temp = alpha * p;
        x_temp = x + s_temp;
        if f(x_temp) > f(x)
            alpha = alpha / 2;
        end
        
        s = alpha * p;
        x_next = x + s;
        y = grad(x_next(1), x_next(2)) - g;
        
        % calculate the update term.
        H = H - ((H*s)*(H*s)')/(s'*H*s) + (y*y') / (s'*y);
        
        % find next iteration point.
        x = x + s;
        k = k + 1;
    end
    y_opt = f(x);
end