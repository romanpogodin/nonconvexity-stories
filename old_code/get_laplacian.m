function laplacian_matrix = get_laplacian(type, first_param, second_param)
assert(ischar(type))

if strcmp(type, 'random')
    assert(nargin == 3)
    adj_matrix = random('bino', 1, first_param, [second_param, second_param]);
    adj_matrix = triu(adj_matrix, 1) + transpose(triu(adj_matrix, 1));
    laplacian_matrix = full(laplacian(graph(adj_matrix)));
    
elseif strcmp(type, 'petersen')
    adj_matrix = zeros(10); % Petersen graph
    numbers = [[1, 1, 1, 1, 2, 2, 2, 3, 3, 4]; ...
               [2, 3, 4, 5, 3, 4, 5, 4, 5, 5]];
    for i = 1:10
        for j = i+1:10
            if length(unique([numbers(:, i), numbers(:, j)]), 1) == 4
                adj_matrix(i, j) = 1;
                adj_matrix(j, i) = 1;
            end
        end
    end
    laplacian_matrix = full(laplacian(graph(adj_matrix)));
    
elseif strcmp(type, 'k_nn')
    assert(nargin == 2)
    adj_matrix = cat(2, cat(1, zeros(first_param), ones(first_param)), ...
        cat(1, ones(first_param), zeros(first_param))); %K_{n,n}
    laplacian_matrix = full(laplacian(graph(adj_matrix)));
    
elseif strcmp(type, 'cycle')
    assert(nargin == 2)
    adj_matrix = zeros(first_param);
    for i = 1:first_param-1
        adj_matrix(i, i+1) = 1;
        adj_matrix(i+1, i) = 1;
    end
    
    adj_matrix(1, first_param) = 1;
    adj_matrix(first_param, 1) = 1;
    
    laplacian_matrix = full(laplacian(graph(adj_matrix)));

elseif strcmp(type, 'type5')
    laplacian_matrix = [3, -1, 0, 0, -1, -1;
                -1, 3, -1, -1, 0, 0;
                0, -1, 3, -1, 0, -1;
                0, -1, -1, 3, -1, 0;
                -1, 0, 0, -1, 3, -1;
                -1, 0, -1, 0, -1, 3]; % topology 5
            
elseif strcmp(type, 'type4')     
    laplacian_matrix = [3, -1, 0, 0, -1, -1;
                -1, 3, -1, -1, 0, 0;
                0, -1, 2 , -1, 0, 0;
                0, -1, -1, 3, -1, 0;
                -1, 0, 0, -1, 3, -1;
                -1, 0, 0, 0, -1, 2]; % topology 4
            
elseif strcmp(type, 'type3')
    laplacian_matrix = [3, -1, -1, -1; 
                        -1, 3, -1, -1;
                        -1, -1, 3, -1;
                        -1, -1, -1, 3]; % topology 3
    
elseif strcmp(type, 'type2')
    laplacian_matrix = [2, -1, 0, -1; 
                        -1, 3, -1, -1;
                        0, -1, 2, -1; 
                        -1, -1, -1, 3]; % topology 2
    
elseif strcmp(type, 'type1')
    laplacian_matrix = [2, -1, -1, 0;
                        -1, 3, -1, -1;
                        -1, -1, 3, -1;
                        0, -1, -1, 2]; % topology 1
    
elseif strcmp(type, 'triangle')
    laplacian_matrix = [2, -1, -1;
                        -1, 2, -1;
                        -1, -1 ,2]; % triangle
end
end