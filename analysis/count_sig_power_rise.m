function rise_pts = count_sig_power_rise(main,spike,pre,clinical)

alpha = 0.05;
np = length(spike);

% Initialize p-value array
pval_array = nan(np,2);
agg_array = nan(np,2);
all_change = cell(np,1);


% Loop over patients 
for p = 1:np
    
    %% Get max spike power electrode
    sp_power = spike(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);

    %% Get fake max spike power electrode for non-spike periods
    n_not_spikes = size(main.pt(p).sp_or_not(2).data,1);
    % Randomly sample n_not_spikes from max chs
    fake_max_ch = datasample(max_ch,n_not_spikes);

    %% Get fake pre-times for non-spike periods
    for s = 1:n_not_spikes
        pre(p).fake_before_rise(s,:) = mode(pre(p).before_rise,1);
    end
    
    for sp = 1:2
        data = main.pt(p).sp_or_not(sp).data;
        main_ch_data = nan(size(data,1),size(data,2));
        first = nan(size(data,1),1);
        last = nan(size(data,1),1);

        % Loop over spikes
        for s = 1:size(data,1)

            if sp == 1
                % Get the power in the biggest spike power electrode
                main_ch_data(s,:) = data(s,:,max_ch(s));
                current_pre = pre(p).before_rise(s,:);
            else
                % Get the power in the biggest fake power electrode
                main_ch_data(s,:) = data(s,:,fake_max_ch(s));
                current_pre = pre(p).fake_before_rise(s,:);
            end

            % remove non-pre times
            last_pre = find(current_pre ==0);
            last_pre = last_pre(1) - 1; % the last pre-spike window
            main_ch_data(s,current_pre==0) = nan;
            

            first(s) = main_ch_data(s,1);
            last(s) = main_ch_data(s,last_pre);

        end
        
        if sp == 1
            all_change{p} = last-first;
        end
        
        % Paired t-test comparing first and last time
        [~,pval,~,stats] = ttest(last,first);
        pval_array(p,sp) = pval;
        
        agg_array(p,sp) = stats.tstat;
        
        % investigation
        %{
        if sp == 1 && p == 13, error('look'); end
        for t = 1:size(main_ch_data,2)
            plot(t,main_ch_data(:,t),'o');
            hold on
        end
        %}
    end
    
end

%% Plot
% Turn this into a violin plot
for i = 1:length(all_change)
    plot(i+0.05*rand(length(all_change{i}),1),all_change{i},'ko')
    hold on
    plot(xlim,[0 0],'k--')
end

fprintf('\n%d of %d patients had a significant pre-IED rise.\n',...
    sum(pval_array(:,1)<alpha),size(pval_array,1));

fprintf('\n%d of %d patients had a significant rise in IED-free periods.\n',...
    sum(pval_array(:,2)<alpha),size(pval_array,1));

% Do aggregate test
[~,pval_spike] = ttest(agg_array(:,1));
[~,pval_not] = ttest(agg_array(:,2));

fprintf('\nThe aggregate p-vals are %1.3f for IED and %1.3f for IED-free.\n',...
    pval_spike,pval_not);

rise_pts = find(pval_array(:,1)<alpha);

%% Do some clinical correlations
% Build clinical table
names = cell(np,1);
ageSurgery = cell(np,1);
ageOnset = cell(np,1);
sex = cell(np,1);
localization = cell(np,1);
pathology = cell(np,1);
lesional = cell(np,1);
for p = 1:np
    names{p} = clinical(p).name;
    ageSurgery{p} = (clinical(p).clinical.ageSurgery);
    ageOnset{p} = (clinical(p).clinical.ageOnset);
    sex{p} = clinical(p).clinical.sex;
    localization{p} = clinical(p).clinical.seizureOnset;
    pathology{p} = clinical(p).clinical.pathology;
    lesional{p} = clinical(p).clinical.lesionstatus;
end

table(names,pval_array(:,1),ageSurgery,ageOnset,sex,localization,pathology,lesional)
% Doesn't look like it matters


end