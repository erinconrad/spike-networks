function count_sig_pts(metrics,met)
alpha = 0.05;
nfreq = 3;
for f = 1:nfreq
  ps = metrics.time.freq(f).(met).auc.individual_ranksum_p;  
  fname = metrics.time.freq(f).name;
  alpha_adj = alpha/length(ps);
  n_sig = sum(ps<alpha);
  fprintf('\nFor %s %s, %d of %d had significant pre-IED rise.\n',...
      fname,met,n_sig,length(ps));
    
end


end