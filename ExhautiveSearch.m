function [ui, k] = ExhautiveSearch(f, st, tol)
L0 = st(2) - st(1);
N = 2*L0 / tol;
x =  linspace(st(1), st(2), N+2);
fe = f(x);
[~, k] = min(fe);
ui = x([k-1, k+1]);
end