function ns = norm_schatten(X, p, eps)

if nargin < 3
    eps = 0.1;
end

if nargin < 2
    p = 1.0;
end

ns = sum((svd(X) + eps) .^ p) ^ (1 / p);

end