function spike_classifier(stats,sig_dev)

which_time_windows = 1;
which_networks = 1; %1 = coherence, 2 = simple
which_freqs = 1:8; 

all_feats.s = [];
all_feats.ns = [];

for n = (which_networks)

    for t = (which_time_windows)
        
        temp_times = stats(n).time(t).index_windows(:,1)/stats(n).time(t).fs-3;
        times = realign_times(temp_times,nan);
        % Fix rounding error
        times = round(times*1e2)/(1e2);

        % Find the corresponding sig_dev time window index
        for tt = 1:length(sig_dev)
            if strcmp(sig_dev(tt).name,stats(n).time(t).name) == 1
                sig_dev_idx = tt;
                break
            end
        end

        % The times in the signal deviation structure should line up with
        % the times in the network structure
        sig_power_change_bin = sig_dev(sig_dev_idx).sig;
        sig_dev_times = sig_dev(sig_dev_idx).times;
        if ~isequal(sig_dev_times,times)
            error('Non-aligning times');
        end

        
        
        for f = (which_freqs)
            F_curr = stats(n).time(t).freq(f).F_all;
             % Just take times without significant power change
            F_curr = F_curr(:,~sig_power_change_bin,:);
            
            % z score to normalize within pt
            z_curr = (F_curr-nanmean(F_curr,2))./nanstd(F_curr,0,2); %nan because first time is nans
            slopes = zeros(size(z_curr,1),2);
            
            % loop over patients to get slopes
            for i = 1:size(z_curr,1)
                for s = 1:2
                    % Get slope for each patient
                    y = squeeze(z_curr(i,:,s))';
                    x = [ones(length(y),1), (1:length(y))'];
                    % do regression to find best fit line through the F stats
                    % for that patient
                    b = x\y; 
                    slopes(i,s) = b(2);
                end
                
            end
            
            all_feats.s = [all_feats.s,slopes(:,1)];
            all_feats.ns = [all_feats.ns,slopes(:,2)];
            
        end

    end
end

%% Logistic regression
% binary for s vs ns
s_or_no = [ones(length(all_feats.s),1);zeros(length(all_feats.ns),1)];
s_or_no = logical(s_or_no);
mdl = glmfit([all_feats.s;all_feats.ns],s_or_no,'binomial','link','logit');

end