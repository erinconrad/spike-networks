function spike_info(times)

np = 0;
nspikes = [];
for i = 1:length(times)
    
    if isempty(times(i).spike_times) == 1
        continue
    else
        np = np + 1;
        nspikes = [nspikes;sum(~isnan(times(i).spike_times))];
    end
    
    
end

fprintf('There are %d patients\n',np);

end