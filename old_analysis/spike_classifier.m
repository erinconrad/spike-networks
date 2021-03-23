function spike = spike_classifier(stats)

which_time_windows = 1;
which_networks = 1; %1 = coherence, 2 = simple
which_freqs = 1:8; 

all_feats.s = [];
all_feats.ns = [];

for n = (which_networks)

    for t = (which_time_windows)

        for f = (which_freqs)
            slopes = stats(n).time(t).freq(f).tests.slopes;
            
            
            all_feats.s = [all_feats.s,slopes(:,1)];
            all_feats.ns = [all_feats.ns,slopes(:,2)];
            
        end

    end
end

%% Logistic regression
% binary for s vs ns
s_or_no = [ones(length(all_feats.s),1);zeros(length(all_feats.ns),1)];
s_or_no = logical(s_or_no);

spike = array2table([all_feats.s;all_feats.ns]);
spike.Group = s_or_no;
%[b,dev,stats] = glmfit([all_feats.s;all_feats.ns],s_or_no,'binomial','link','logit');

end