function sim = intra_group_similarity(A)

%{
This function takes an M x N matrix A, where M (rows) are the elements of
the flattened adjacency matrix for a certain spike, and N (columns) are the
number of adjacency matrices, and it outputs sim, a measure of the average
correlation between each pair of N flattened adjacency matrices.
%}

% The pairwise Pearson correlation coefficient between each pair of columns
% in A (produces an NxN matrix)
rho = corr(A);

% Take the upper triangle of the rho matrix (it is a symmetric matrix, and the
% diagonal has no useful information)
upper = triu(rho,1);

% Take the average of the non zero elements (this assumes that none of the
% actual correlations are exactly 0).
elements = upper(upper~=0);
sim = mean(elements);

end