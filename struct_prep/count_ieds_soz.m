function count_ieds_soz(is_spike_soz)

ns = zeros(length(is_spike_soz),2);

for p = 1:length(is_spike_soz)
    n_soz = sum(is_spike_soz(p).is_soz);
    n_not_soz = sum(~is_spike_soz(p).is_soz);
    ns(p,:) = [n_soz,n_not_soz];
    
end

total_ns = sum(ns,1);
all_total = sum(total_ns);
total_soz = total_ns(1);
[~,pval,~,stats] = ttest(ns(:,1),ns(:,2));

fprintf('\n%d of %d (%1.1f%%) IEDs were in the SOZ (p = %1.3f).\n',...
    total_soz,all_total,total_soz/all_total*100,pval);




end