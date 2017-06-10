function grad = compute_schatten_grad(X, p, eps)
% COMPUTE_SCHATTEN_GRAD Computes gradient of Schatten p-norm in the p-th
% power
%   ns = COMPUTE_SCHATTEN_GRAD(X, p, eps) to computer gradient with
%   smoothing parameter eps of p-norm in p-th power
    grad = p * X * mpower(transpose(X) * X + ...
        eps * eye(size(X, 1)), (p - 2.0) / 2.0);
end