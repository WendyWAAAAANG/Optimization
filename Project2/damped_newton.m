function [x, k, y_opt] = damped_newton(f, x0, tol)
    max_iter = 500;
    x = x0;
    % counter -- time of iteration.
    k = 0;
    % initialize beta to 0.
    beta = 0;
    
    % find gradient and hessian matrix.
    syms x1 x2
    vars = [x1; x2];
    % calculate gradient and hessian function of f.
    grad = matlabFunction(gradient(f(vars), vars));
    hess = matlabFunction(hessian(f(vars), vars));
    
    % iteration.
    while norm(grad(x(1), x(2))) > tol && k < max_iter
        %d = -inv(hess(x))*grad(x);
        % for this iteration, calculate gradient and hessian value.
        g = grad(x(1), x(2));
        h = hess(x(1), x(2));
        % find beta that let H be a positive definite matrix.
        H = h + beta*eye(length(x));
        while all(eig(H) > 0) == false
            beta = beta + 1;
            H = h + beta*eye(length(x));
        end
        
        % calculate the update term.
        d = -H \ g;
        
        % initialize learning rate to 1.
        alpha = 1;
        
        % using line search to find suitable learning rate.
        x_temp = x + d*alpha;
        if f(x_temp) > f(x)
            alpha = alpha / 2;
        end
        x = x + d*alpha;
        k = k + 1;
    end
    y_opt = f(x);
end