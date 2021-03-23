function sig_dev = assess_power_change(sig_dev,alpha)

time_count = length(sig_dev);

%% Get slopes for network change over time

for t = 1:time_count
    sig_dev(t).tests.slopes = zeros(size(sig_dev(t).is_spike(1).t_stat_all,1),2);
    for s = 1:2 % spike and not spike
        
        if s>length(sig_dev(t).is_spike), continue; end
        
        % get t-stats
        tstats = sig_dev(t).is_spike(s).t_stat_all;

        % z score to normalize within pt
        z_curr = (tstats-mean(tstats,2))./std(tstats,0,2);
        
        sig_dev(t).is_spike(s).z_curr = z_curr;

        % initialize slopes
        slopes = zeros(size(z_curr,1),1); % n_pt x 1 

        % Loop over patients and get slopes
        np = size(z_curr,1);
        for i = 1:np

            % Get z scores over times
            z = squeeze(z_curr(i,:))';

            % Linear regression to get best fit line
            x = [ones(length(z),1), (1:length(z))'];
            b = x\z;
            slopes(i) = b(2); % slope of line
            
        end
        if length(slopes) ~= size(sig_dev(t).tests.slopes,1)
            fprintf('Warning, inconsistent number of patients.\n');
            continue;
        end
        sig_dev(t).tests.slopes(:,s) = slopes;
    end
end

%% Paired t-test comparing slopes in spike and not-spike
for t = 1:time_count
    slopes = sig_dev(t).tests.slopes;
    [~,p,~,stats1] = ttest(slopes(:,1),slopes(:,2)); 
    sig_dev(t).tests.paired.p = p;
    sig_dev(t).tests.paired.stats = stats1;
    adj_alpha = alpha;
    sig_dev(t).tests.paired.h = p < adj_alpha;
end

%% Unpaired t-test looking at slopes in spike
for t = 1:time_count
    slopes = sig_dev(t).tests.slopes;
    for s = 1:2
        [~,p,~,stats1] = ttest(slopes(:,s)); 
        sig_dev(t).tests.unpaired.p(s) = p;
        sig_dev(t).tests.unpaired.stats(s).stats = stats1;
        adj_alpha = alpha;
        sig_dev(t).tests.unpaired.h(s) = p < adj_alpha;
    end
end


end