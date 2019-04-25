function sync = synchronizability_sp(A)

%% Calculate node strength D
% Get the degree vector of the adjacency matrix
dvector = sum(A,1);

% convert it into a diagonal matrix D
D = diag(dvector);

%% Calculate Laplacian L
L = D - A;

%% Compute the eigenspectrum of L(t)
e = eig(L);

%% Calculate synchronizabilty s

% Sort smallest to largest eigenvalue
e = sort(e);

% Compute synchronizability: ratio of second smallest eigenvalue to largest
% eigenvalue
sync = abs(e(2)/e(end));


end