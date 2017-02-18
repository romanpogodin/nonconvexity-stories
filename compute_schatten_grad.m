function grad = compute_schatten_grad(matrix, p, eps)
    grad = matrix * mpower(transpose(matrix) * matrix + ...
        eps * eye(size(matrix, 1)), (p - 2.0) / 2.0);
end