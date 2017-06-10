function ns = norm_singval(X, q, eps, make_unit_diag, diag_index)
% NORM_SINGVAL Computes singular value rank approximation
%   ns = NORM_SINGVAL(X, q, eps) for rank(X) approximaton with smoothing
%   parameter eps and multiplier (1 + eps ^ q)
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
    q = 1.0;
end

if make_unit_diag
   X(diag_index) = 1; 
end

ns = (1 + eps ^ q) * (size(X, 1) - eps * sum((svd(X * transpose(X)) + eps) .^ -1));

end