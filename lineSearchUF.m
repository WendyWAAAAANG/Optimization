function [ui, k] = lineSearchUF(f, st, tol)
%% function[ui, k] = lineSearchUF(f, st, tol) is the function.
%search method for unimode function.
% Input:
% Output: 

% SEE ALSO: newton
% copyright @ shengxin ZHU

% reference: https://www.bilibili.com/read/cv20451589

IntervalLength = st(2) - st(1);
k = 1;
while IntervalLength > tol
    tau = 0.25  * IntervalLength;
    x1 = st(1) + tau;
    x2 = st(2) - tau;
    if f(x1) < f(x2)
        st(2) = x2;
    else
        if f(x1) > f(x2)
            st(1) = x1;
        else
            st(1) = x1;
            st(2) = x2;
        end
    end
    IntervalLength = st(2) - st(1);
    k = k+1;
end
ui = st;
end