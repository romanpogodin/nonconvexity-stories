function objective_matrix = read_bigmac_problem(problem_size, problem_number)
%READ_BIGMAC_PROBLEM Reads a Bigmac 0-1 problem
%   objective_matrix = READ_BIGMAC_PROBLEM(problem_size, problem_number) 
%   reads a 0-1 min problem, invertes it to make a max problem    
%   http://biqmac.uni-klu.ac.at/library/biq/beasley/bqp#1-#2.sparse, where #1 =
%   problems_size, #2 = problem_number

url = strcat('http://biqmac.uni-klu.ac.at/library/biq/beasley/bqp', ...
    int2str(problem_size), '-', int2str(problem_number), '.sparse');
options = weboptions('ContentType', 'text'); % does not work with table type
data = webread(url, options);

data = strsplit(data); 
data = cellfun(@str2num, data(1:length(data) - 1));

data = data(3:end);
data = transpose(reshape(data, 3, length(data) / 3));

objective_matrix = sparse(data(1:end, 1), ...
    data(1:end, 2), ...
    data(1:end, 3), ...
    problem_size, problem_size);

objective_matrix = objective_matrix + objective_matrix.' - ...
    diag(diag(objective_matrix));

objective_matrix = -objective_matrix; % to have a max problem
end
