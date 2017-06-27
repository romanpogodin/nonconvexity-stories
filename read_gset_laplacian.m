function graph_laplacian = read_gset_laplacian(graph_number)
%READ_GSET_LAPLACIAN Reads Gset graph into a Laplacian matrix
%   graph_laplacian = READ_GSET_LAPLACIAN(graph_number) reads a graph from 
%   https://web.stanford.edu/~yyye/yyye/Gset/G#, where # = graph_number,
%   into a Laplacian matrix

url = strcat('https://web.stanford.edu/%7Eyyye/yyye/Gset/G', ...
    int2str(graph_number));
options = weboptions('ContentType','table');
data = webread(url, options);
graph_laplacian = full(laplacian(graph(table2array(data(2:end, 1)), ...
    table2array(data(2:end, 2)))));
end

