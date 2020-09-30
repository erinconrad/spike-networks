function [sig_dev,stats,metrics] = reduce_to_ns(sig_dev,stats,metrics,paired)

network_count = length(stats);
time_count = length(stats(1).time);


% First confirm that times align
for n = 1:network_count
    for t = 1:time_count
        if t > length(sig_dev), continue; end
        
        perm_times = stats(n).time(t).times;
        sd_times = sig_dev(t).is_spike(1).times;
        
        if ~isequal(perm_times,sd_times), error('Non-aligning times'); end
    end
end

% Now reduce power
for t = 1:time_count
        
    if t > length(sig_dev), continue; end
    
    % get the indices of significant power change
    sig = sig_dev(t).tests.(paired).sig; % only care about spikes, not not-spikes

    % Reduce time and power change to non-significant indices
    for s = 1:2 % reduce both spike and not spike times based on spike power
        if s>length(sig_dev(t).is_spike), continue; end
        sig_dev(t).is_spike(s).times = sig_dev(t).is_spike(s).times(~sig);
        sig_dev(t).is_spike(s).t_stat_all = sig_dev(t).is_spike(s).t_stat_all(:,~sig);
    end
    
end

% Now reduce perm
for n = 1:network_count
    nfreq = length(stats(n).time(1).freq);
    for t = 1:time_count
        
        if t > length(sig_dev), continue; end
        
        % get the indices of significant power change
        sig = sig_dev(t).tests.(paired).sig;
        
        % reduce time to non-significant indices
        stats(n).time(t).times = stats(n).time(t).times(~sig);
        
        for f = 1:nfreq
            % reduce perm to non-significant indices
            stats(n).time(t).freq(f).F_all = stats(n).time(t).freq(f).F_all(:,~sig,:);
            
        end
    end
end

% Now reduce ns
for n = 1:network_count
    nfreq = length(metrics(n).time(1).freq);
    for t = 1:time_count
        
        if t > length(metrics(n).time), continue; end
        
        % get the indices of significant power change
        sig = sig_dev(t).tests.(paired).sig;
        
        if isempty(metrics(n).time(t).times), continue; end
        
        % reduce time to non-significant indices
        metrics(n).time(t).times = metrics(n).time(t).times(~sig);
        
        for f = 1:nfreq
            % reduce metrics to non-significant indices
            metrics(n).time(t).freq(f).ge.data = metrics(n).time(t).freq(f).ge.data(:,~sig,:);
            metrics(n).time(t).freq(f).ns_avg.data = metrics(n).time(t).freq(f).ns_avg.data(:,~sig,:);
            metrics(n).time(t).freq(f).ns_big.data = metrics(n).time(t).freq(f).ns_big.data(:,~sig,:);
            
        end
    end
end

end