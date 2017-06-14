mat_size = [7, 7];
diagonal_index = logical(eye(mat_size));
std_err = 1e-8;

% for i = 1:100000
    w = triu((rand(mat_size) - 0.5) * 1.1, 1);
    w = w + transpose(w);
    w(diagonal_index) = 1;

    [~, is_negative] = cholcov(w);
    
%     if ~is_negative
        counter = 0;
        for j = 1:10000
           noise =  abs(triu(randn(mat_size), 1) * std_err);
           if min(eig(sdp + noise + transpose(noise))) >= 0
              counter = 1 + counter; 
           end
        end
        disp(counter / 10000);
%     end
% end