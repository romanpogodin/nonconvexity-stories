function grad = compute_singval_grad(matrix, q, eps)
    grad = 2 * eps * (1 + eps ^ q) * mpower(matrix * transpose(matrix) + ...
        eps * eye(size(matrix, 1)), -2) * matrix;
end