function stats = assess_net_change(stats,alpha)

network_count = length(stats);
time_count = length(stats(1).time);

%% Get slopes for network change over time
for n = 1:network_count
    nfreq = length(stats(n).time(1).freq);
    for t = 1:time_count
        
        for f = 1:nfreq
            
            if isempty(stats(n).time(t).freq), continue; end
            
            % get pseudo F statistics
            F_curr = stats(n).time(t).freq(f).F_all;
            
            % Skip if incomplete
            if ndims(F_curr) < 3, continue; end
            
            % z score to normalize within pt
            z_curr = (F_curr-mean(F_curr,2))./std(F_curr,0,2);
            
            
            % initialize slopes
            slopes = zeros(size(z_curr,1),2); % n_pt x 2 (spike and not spike)
            
            % Loop over patients and get slopes
            np = size(z_curr,1);
            for i = 1:np
                for s = 1:2 % spike and not_spike
                    
                    % Get z scores over times
                    z = squeeze(z_curr(i,:,s))';
                    
                    % Linear regression to get best fit line
                    x = [ones(length(z),1), (1:length(z))'];
                    b = x\z;
                    slopes(i,s) = b(2); % slope of line
                
                end
            end
            stats(n).time(t).freq(f).z_curr = z_curr;
            stats(n).time(t).freq(f).tests.slopes = slopes;
        end
    end
end

%% Paired t-test comparing slopes in spike and not-spike
total_n_freq = 1+length(stats(1).time(1).freq);
for n = 1:network_count
    nfreq = length(stats(n).time(1).freq);
    for t = 1:time_count
        if isempty(stats(n).time(t).freq), continue; end
        for f = 1:nfreq
            if ~isfield(stats(n).time(t).freq(f),'tests'), continue; end
            
            slopes = stats(n).time(t).freq(f).tests.slopes;
            [~,p,~,stats1] = ttest(slopes(:,1),slopes(:,2)); % expect positive t stats because higher slopes in spike
            stats(n).time(t).freq(f).tests.paired.p = p;
            stats(n).time(t).freq(f).tests.paired.stats = stats1;
            
            adj_alpha = alpha/total_n_freq;
            stats(n).time(t).freq(f).tests.paired.h = p < adj_alpha;
            
            % Unpaired test
            for s = 1:2
                [~,p,~,stats1] = ttest(slopes(:,s)); % expect positive t stats because positive slopes
                stats(n).time(t).freq(f).tests.unpaired.p(s) = p;
                stats(n).time(t).freq(f).tests.unpaired.stats(s).stats = stats1;

                adj_alpha = alpha/total_n_freq;
                stats(n).time(t).freq(f).tests.unpaired.h(s) = p < adj_alpha;
            end
                            
        end
    end
end

end