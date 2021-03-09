function [all,pt_all] = compare_two_reviewers(pre_spike)

all = [];
first= [];

for p = 1:length(pre_spike)
    both = pre_spike(p).windows.before_rise_both_times;
    all = [all;both];
    first = [first;min(both,[],2)];
    pt_all(p).both = both;
end

if 0
figure
histogram(all(:,1)-all(:,2))
xlabel('Difference in time of earliest rise (seconds)')
end


fprintf('\nReviewer 1 marked the spike rise on average %1.1f ms before reviewer 2.\n',...
    -mean(all(:,1)-all(:,2))*1e3)

fprintf('\nThe average earliest pre-spike rise was M = %1.1f ms (SD = %1.1f ms).\n',...
    mean(first)*1e3,std(first)*1e3);
end