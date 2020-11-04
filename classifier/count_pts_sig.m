function count_pts_sig(metrics,met,n,f,t,alpha)
    
    curr_met = metrics(n).time(t).freq(f).(met);
    all_p = [];
    for i = 1:length(curr_met.pt)
        all_p = [all_p; metrics(n).time(t).freq(f).(met).pt(i).test.p];
    end

    n_sig = sum(all_p < alpha);
    n_pts = length(all_p);
    fprintf('\nOf %d patients, %d had significant pre-spike rise in %s.\n',n_pts,n_sig,met);
end