function [rise_pts,all_change] = count_sig_power_rise(main,spike,pre,clinical,out_folder)

alpha = 0.05;
np = length(spike);

% Initialize p-value array
pval_array = nan(np,2);
agg_array = nan(np,2);
all_change = cell(np,1);
first_last = cell(np,1);
first_last_no = cell(np,1);
mean_first_last = nan(np,2);


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
            first_last{p} = [first,last];
            mean_first_last(p,:) = [mean(first) mean(last)];
        elseif sp == 2
            first_last_no{p} = [first,last];
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




%% Stats
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

%% Plot figure showing aggregate across patients
if 0
    cols = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980];
    figure
    for i = 1:2
        plot(i+0.05*rand(size(mean_first_last,1),1),mean_first_last(:,i),'o',...
            'markersize',15)
        hold on
    end
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Baseline','pre-IED'})
    ylabel('Mean power (\muV^2)');
    sigstar([1 2],pval_spike)
    set(gca,'fontsize',20)
    %print(gcf,[out_folder,'all_pt_change'],'-dpng');
end

%% Plot showing each patient individually
if 1
    for sp = [1:2]
        if sp == 1
            thing = first_last;
        else
            thing = first_last_no;
        end
    figure
    set(gcf,'position',[440 462 988 336])
    cols = [0, 0.4470, 0.7410;...
        0.8500, 0.3250, 0.0980];
    for i = 1:length(first_last)
        fp = errorbar(i-0.1,mean(thing{i}(:,1)),std(thing{i}(:,1)),'o',...
            'color',cols(1,:),'markersize',15,'linewidth',2);
        hold on
        lp = errorbar(i+0.1,mean(thing{i}(:,2)),std(thing{i}(:,2)),'o',...
            'color',cols(2,:),'markersize',15,'linewidth',2);
    end
    if sp == 1
        legend([fp lp],{'Baseline','pre-IED'},'fontsize',20,'location','southwest')
    else
        legend([fp lp],{'Baseline','pre-IED'},'fontsize',20,'location','northwest')
    end
    xticklabels([])
    xlabel('Patient')
    ylabel('Power (\muV^2)')
    set(gca,'fontsize',20)
    yl = ylim;
    ylim([yl(1) yl(2)*1.1]);
    yl = ylim;
    for i = 1:length(first_last)
        text_out = get_asterisks(pval_array(i,sp),1);
        text(i,yl(1) + 0.9*(yl(2)-yl(1)),text_out,'fontsize',35,...
            'horizontalalignment','center')
    end
    if sp == 1
        print(gcf,[out_folder,'all_pt_change'],'-dpng');
    else
        print(gcf,[out_folder,'all_pt_change_no_ied'],'-dpng');
    end
    end
end


%% Plot showing heterogeneity in spikes across rise patients
if 0
figure
set(gcf,'position',[440 342 576 456])
for i = 1:length(rise_pts)
    p = rise_pts(i);
    plot(i+0.05*rand(length(all_change{p}),1),all_change{p},'ko')
    hold on
    
    
end
xticklabels([])
xlim([0 length(rise_pts)+1])
plot(xlim,[0 0],'k--')
xlabel('Patients with pre-IED rise')
ylabel('Pre-IED - baseline power difference (\muV^2)')
set(gca,'fontsize',20)
print(gcf,[out_folder,'sp_heterogeneity'],'-dpng');
end



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

pretty_p = arrayfun(@(x) sprintf('%1.3f',x),pval_array(:,1),'UniformOutput',false);
pretty_p(strcmp(pretty_p,'0.000')) = {'<0.001'};
t = table(names,ageSurgery,ageOnset,sex,localization,pathology,lesional,pretty_p);
writetable(t,[out_folder,'clinical_table.csv']);
% Doesn't look like it matters


end