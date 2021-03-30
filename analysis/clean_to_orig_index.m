function orig_index = clean_to_orig_index(p,s,sp_or_not,bad)

%{
This function converts cleaned indices to original indices
%}

bad_indices = bad.pt(p).sp_or_not(sp_or_not).bad;
not_bad_indices = find(bad_indices==0);
orig_index = not_bad_indices(s);


end