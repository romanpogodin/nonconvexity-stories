function ns = norm_singval(X, q, eps)

if nargin < 3
    eps = 0.0;
end

if nargin < 2
    q = 1.0;
end

ns = (1 + eps ^ q) * (size(X, 1) - eps * sum((svd(X * transpose(X)) + eps) .^ -1));

end