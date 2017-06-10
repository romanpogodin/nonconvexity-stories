function ns = norm_schatten(X, p, eps, make_unit_diag, diag_index)
% NORM_SCHATTEN Computes Schatten norm of a matrix
%   ns = NORM_SCHATTEN(X, q, eps) for rank(X) approximaton with smoothing
%   parameter eps and norm parameter p
%   Other parameters are:
%   make_unit_diag -- whether to force diagonal to be unit
%   diag_index -- if make_unit_diag == true, pass indices of the diagonal
%   elements (speeds up computations)

if nargin < 4
    make_unit_diag = false;
end

if nargin < 3
    eps = 0.0;
end

if nargin < 2
    p = 1.0;
end

if make_unit_diag
   X(diag_index) = 1; 
end

ns = sum((svd(X * transpose(X)) + eps) .^ (p / 2)) ^ (1 / p);

end