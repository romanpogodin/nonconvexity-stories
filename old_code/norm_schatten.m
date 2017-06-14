function ns = norm_schatten(X, p, eps)

if nargin < 3
    eps = 0.0;
end

if nargin < 2
    p = 1.0;
end

ns = sum((svd(X * transpose(X)) + eps) .^ (p / 2)) ^ (1 / p);

end