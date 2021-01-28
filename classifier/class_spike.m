function tbl=class_spike(stats,met)

nfreq = length(stats.time.freq);
np = length(stats.time.freq(1).(met).pt);


%% Aggregate the normalized metrics for all spikes, all patients, all freq, spike and not a spike
for f = 1:nfreq
    freq(f).name = stats.time.freq(f).name;
    freq(f).all_sp = [];
    freq(f).all_no = [];
    for p = 1:np
        met_sp = stats.time.freq(f).(met).pt(p).spike.auc;
        met_no = stats.time.freq(f).(met).pt(p).not.auc;
        
        %% Get the last non nan time
        all_nan_columns = sum(isnan(met_sp),1) == size(met_sp,1);
        last_non_nan = find(all_nan_columns);
        last_non_nan(last_non_nan == 1) = [];
        if isempty(last_non_nan)
            last_non_nan = size(met_sp,2);
        else
            last_non_nan = last_non_nan(1)-1;
        end
        met_sp = met_sp(:,last_non_nan);
        met_no = met_no(:,last_non_nan);
        
        freq(f).all_sp = [freq(f).all_sp;met_sp];
        freq(f).all_no = [freq(f).all_no;met_no];
        
    end
    
    %% Shuffle the spikes and non spikes
    freq(f).all_sp = freq(f).all_sp(randperm(length(freq(f).all_sp)));
    freq(f).all_no = freq(f).all_no(randperm(length(freq(f).all_no)));
    
    %% Get identities
    identities = cell(length(freq(f).all_sp)+length(freq(f).all_no),1);
    identities(1:length(freq(f).all_sp)) = {'spike'};
    identities(length(freq(f).all_sp)+1:end) = {'not'};
    
    %% Put it into a table
    if f == 1
        tbl = table(identities);
    end
    
    tbl = addvars(tbl,[freq(f).all_sp;freq(f).all_no],...
        'NewVariableNames',freq(f).name);
    
    
end






end