function rel_change = calc_rel_change(data,before_rise,rel_times)


nspikes = size(data,1);
rel_change = nan(nspikes,1);

for s = 1:nspikes
    % get pre-spike times
    br = before_rise(s,:);
    br_data = data(s,br == 1);
    
    last = br_data(end);
    baseline = mean(br_data(rel_times));
    
    rel_change_current = (last-baseline)/abs(baseline);
    rel_change(s) = rel_change_current;
end

end