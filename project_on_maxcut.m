function new_x = project_on_maxcut(curr_x, laplacian_matrix, ...
    cut_optval, sdp_optval, is_cvx_quiet, is_constraint_relaxed)
%PROJECT_ON_MAXCUT Make a Frobenius norm projection
if is_constraint_relaxed
    if is_cvx_quiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end
    variable Y(size(curr_x)) symmetric
        minimize norm(Y - curr_x, 'fro')
        subject to
            Y >= 0
            diag(Y) == 1
            4 * cut_optval <= trace(laplacian_matrix * Y) <= 4 * sdp_optval
    cvx_end
    new_x = Y;
else
    if is_cvx_quiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end
    variable Y(size(curr_x)) symmetric
        minimize norm(Y - curr_x, 'fro')
        subject to
            Y >= 0
            diag(Y) == 1
            trace(laplacian_matrix * Y) == 4 * sdp_optval
    cvx_end
    new_x = Y;
end
end

