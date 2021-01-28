function all = compare_two_reviewers(pre_spike)

all = [];

for p = 1:length(pre_spike)
    both = pre_spike(p).windows.before_rise_both_times;
    all = [all;both];

end

if 0
figure
histogram(all(:,1)-all(:,2))
xlabel('Difference in time of earliest rise (seconds)')
end


fprintf('\nReviewer 1 marked the spike rise on average %1.1f ms before reviewer 2.\n',...
    -mean(all(:,1)-all(:,2))*1e3)
end