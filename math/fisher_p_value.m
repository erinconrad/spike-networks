function comb_p = fisher_p_value(all_p)

X_2 = -2 * sum(log(all_p));
comb_p = 1-chi2cdf(X_2,2*length(all_p));

end