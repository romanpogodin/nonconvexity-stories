function grad = compute_singval_grad(X, q, eps)
% COMPUTE_SINGVAL_GRAD Computes gradient of singular value rank
% approximation
%   ns = COMPUTE_SINGVAL_GRAD(X, p, eps) to computer gradient with
%   smoothing parameter eps and multiplier (1 + eps ^ q)
    grad = 2 * eps * (1 + eps ^ q) * mpower(X * transpose(X) + ...
        eps * eye(size(X, 1)), -2) * X;
end