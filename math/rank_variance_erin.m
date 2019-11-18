function std_rank = rank_variance_erin(X)

%{
This function takes a matrix X, in which each column of X has a ranking of
elements (N elements represented by rows), and returns a std std_rank,
representing the std in the ranks
%}


% calculate the standard deviation across each row
std_rows = std(X,0,2);

% Take the average standard deviation
std_rank = mean(std_rows);

end