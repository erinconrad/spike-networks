function detailed_sp_analyses(rise_pts,spike_power,pre_spike,abs_power,ers,ns,pearson,out_folder)


alpha = 0.05;
nfreq = 3;

%% Locations
locations = spike_network_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
pt_folder = [data_folder,'spike_structures/'];

%% Load pt structure
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;


%% Get pre rise times
pre_rise_times = nan(length(rise_pts),1);
for i = 1:length(rise_pts)
    p = rise_pts(i);
    both = pre_spike(p).before_rise_both_times;
    pre_rise_times(i) = mean(min(both,[],2));
end

%% Get locs
all_locs = cell(length(rise_pts),1);
for i = 1:length(rise_pts)
    p = rise_pts(i);
    locs = [];
    name = abs_power.freq(1).pt(p).name;
    for j = 1:length(pt)
        if strcmp(pt(j).name,name)
            locs = pt(j).new_elecs.locs;
        end
    end
    all_locs{i} = locs;
end

%% When does the change occur?
rel_change_thresh = -inf;
all_earliest = nan(length(rise_pts),1);
times_ex = pre_spike(1).times(pre_spike(1).times<=0);
all_dat_for_plot = nan(length(rise_pts),60,length(times_ex));
% Loop over rise pts
for i = 1:length(rise_pts)
    
    p = rise_pts(i);

    %% Get max spike power electrode
    sp_power = spike_power(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);

    %% Get before rise data
    before_rise = pre_spike(p).before_rise;
    
    % Get last time before any rise
    any_before_rise = sum(any(before_rise,1));

    %% Find those spikes with a large relative power change
    % Get spike abs power data
    data = abs_power.freq(1).pt(p).sp_or_not(1).data;
    
    % get the data on the channel with the highest spike power change
    data_main_ch = get_data_from_spec_ch(data,max_ch);
    
    % calcaulte the relative power change compared to baseline
    rel_change = calc_rel_change(data_main_ch,before_rise,1:5);
    
    % decide if it's above a threshold
    change_sp_idx = find(rel_change > rel_change_thresh);
    

    % get the data occurring before the rise (setting all other times to
    % nans)
    data_change_sp_idx = data_main_ch(change_sp_idx,:);
    before_rise_change_sp_idx = before_rise(change_sp_idx,:);
    all_data_before_rise = data_change_sp_idx;
    all_data_before_rise(before_rise_change_sp_idx==0) = nan; % set times after the rise to nans
    
    %nsp = size(abs_power.freq(1).pt(p).sp_or_not(1).data,1);
    all_dat_for_plot(i,1:size(data_main_ch,1),:) = all_data_before_rise(:,1:length(times_ex));
    
    time_pvals = nan(size(all_data_before_rise,2),1);
    % Loop over times to see which have significant difference from first
    % time
    for t = 2:size(all_data_before_rise,2)
        
        % Do paired ttest
        [~,pval] = ttest(all_data_before_rise(:,1),all_data_before_rise(:,t));
        time_pvals(t) = pval;
    end

    % See which have p < 0.05 for it and all subsequent ones
    earliest = nan;
    less_alpha = time_pvals < alpha;
    for t = 2:any_before_rise
        num_left = any_before_rise - t + 1; % how many left including this one
        num_sub_alpha = sum(less_alpha(t:length(less_alpha)));
        if num_left == num_sub_alpha
            earliest = t;
            break
        end
    end
    
    if isnan(earliest)
        earliest_time = nan;
    else
        earliest_time = abs_power.times(earliest);
    end
    
    all_earliest(i) = earliest_time;
end

%% Display earliest time of rise
if 0
    all_earliest
    figure
    set(gcf,'position',[100 100 700 1000]);
    [ha,~] = tight_subplot(length(rise_pts),1,[0.045 0.01],[0.07 0.04],[0.085 0.01]);
    for i = 1:length(rise_pts)  
        axes(ha(i))
        errorbar(times_ex,nanmean(squeeze(all_dat_for_plot(i,:,:)),1),...
            nanstd(squeeze(all_dat_for_plot(i,:,:)),[],1),'k');  
        hold on
        
        % get sig times
        sig_times = find(times_ex >= all_earliest(i) & times_ex < -0.1);
        yl = ylim;
        ylim([yl(1) (yl(2)-yl(1))*1.1])
        for t = 1:length(sig_times)
            text(times_ex(sig_times(t)),yl(1)+(yl(2)-yl(1))*0.95,'*','fontsize',30,...
                'horizontalalignment','center')
        end
        
        plot([pre_rise_times(i) pre_rise_times(i)],ylim,'r--','linewidth',2)
        xlim([-2 0])
        if i ~= length(rise_pts)
            xticklabels([])
        else
            xlabel('Time (s)');
        end
        if i == round(length(rise_pts)/2)
            ylabel(('Mean power (\muV^2)'));
        end
        set(gca,'fontsize',20)
    end
    print(gcf,[out_folder,'timing'],'-dpng')
end
        
%% What frequencies are involved in the change
change_by_f = cell(length(rise_pts),1);
first_last_by_f = cell(length(rise_pts),1);
for i = 1:length(rise_pts)
    p = rise_pts(i);
    change_by_f{i} = nan(size(abs_power.freq(1).pt(p).sp_or_not(1).data,1),nfreq);
    first_last_by_f{i} = nan(size(abs_power.freq(1).pt(p).sp_or_not(1).data,1),nfreq,2);
end
for f = 1:nfreq
    for i = 1:length(rise_pts)
        
        p = rise_pts(i);

        %% Get max spike power electrode
        sp_power = spike_power(p).spike_powers;
        [~,max_ch] = max(sp_power,[],2);

        %% Get before rise data
        before_rise = pre_spike(p).before_rise;

        %% Find those spikes with a large relative power change
        % Get spike abs power data
        data = abs_power.freq(1).pt(p).sp_or_not(1).data;

        % get the data on the channel with the highest spike power change
        data_main_ch = get_data_from_spec_ch(data,max_ch);

        % calcaulte the relative power change compared to baseline
        rel_change = calc_rel_change(data_main_ch,before_rise,1:5);

        % decide if it's above a threshold
        change_sp_idx = find(rel_change > rel_change_thresh);

        % get the data occurring before the rise (setting all other times to
        % nans)
        data = ers.freq(f).pt(p).sp_or_not(1).data;
        data_main_ch = get_data_from_spec_ch(data,max_ch);
        data_change_sp_idx = data_main_ch(change_sp_idx,:);
        before_rise_change_sp_idx = before_rise(change_sp_idx,:);
        all_data_before_rise = data_change_sp_idx;
        all_data_before_rise(before_rise_change_sp_idx==0) = nan; % set times after the rise to nans
    
        % loop through each spike and get first and last
        for s = 1:size(all_data_before_rise,1)
            curr = all_data_before_rise(s,:);
            first = curr(1);
            last_idx = find(isnan(curr));
            last_idx = last_idx(1) - 1;
            last = curr(last_idx);
            change_by_f{i}(s,f) = last-first;
            first_last_by_f{i}(s,f,:) = [first,last];
        end
        
        
        
        
    end

           
end



%% individ patient - freq
if 0
    cols = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980];
    figure
    set(gcf,'position',[100 100 800 1000]);
    [ha,~] = tight_subplot(length(rise_pts),3,[0.04 0.04],[0.04 0.04],[0.085 0.01]);
    for i = 1:length(rise_pts)
        curr_dat = first_last_by_f{i};
        for f = 1:size(curr_dat,2)
            axes(ha((i-1)*3+f))
            errorbar(1,mean(curr_dat(:,f,1)),std(curr_dat(:,f,1)),'o',...
                'color',cols(1,:),'linewidth',2)
            hold on
            errorbar(2,mean(curr_dat(:,f,2)),std(curr_dat(:,f,2)),'o',...
                'color',cols(2,:),'linewidth',2)
            xlim([0.5 2.5])
            
            % t-test
            [~,pval] = ttest(curr_dat(:,f,1),curr_dat(:,f,2));
            pp = pretty_p(pval,3); % bonferroni divide by 3
            yl = ylim;
            ylim([yl(1) (yl(2)-yl(1))*1.1]);
            yl = ylim;
            plot([1 2],[yl(1)+(yl(2)-yl(1))*0.8 yl(1)+(yl(2)-yl(1))*0.8],...
                'k');
            text([1.5],[yl(1)+(yl(2)-yl(1))*0.9],...
                pp,'horizontalalignment','center','fontsize',20);
            
            if i == 1
                if f == 1
                    title('Sub-gamma','fontsize',15,'fontweight','normal');
                elseif f ==2
                    title('Low gamma','fontsize',15,'fontweight','normal');
                elseif f==3
                    title('High gamma','fontsize',15,'fontweight','normal');
                end
            end
            
            if i == length(rise_pts)
                xticks([1 2])
                xticklabels({'Baseline','Pre-IED'});
            else
                xticklabels([])
            end
            
            if i == round(length(rise_pts)/2) && f == 1
                ylabel('Mean power (\muV^2)');
            end
            set(gca,'fontsize',20)
        end
    
    end
    print(gcf,[out_folder,'frequency'],'-dpng');
end

% compare with anova - within patient
if 0
mean_changes = nan(length(rise_pts),nfreq);
for i = 1:length(rise_pts)
    changes = change_by_f{i};
    pval = anova1(changes,[],'off');
    
    if 0
        plot(1+0.05*rand(size(change_by_f{i},1),1),change_by_f{i}(:,1),'o')
        hold on
        plot(2+0.05*rand(size(change_by_f{i},1),1),change_by_f{i}(:,2),'o')
        plot(3+0.05*rand(size(change_by_f{i},1),1),change_by_f{i}(:,3),'o')
        title(sprintf('sub-gamma: %1.1e, low-gamma: %1.1e, high-gamma: %1.1e\np = %1.3f',...
            mean(changes(:,1)),mean(changes(:,2)),mean(changes(:,3)),pval))
        pause
        hold off
    end
    
    for f = 1:nfreq
        mean_changes(i,f) = mean(changes(:,f));
    end
    
end


% Compare frequency changes across patients
pval = anova1(mean_changes,[],'off');
if 0
    for f = 1:nfreq
        plot(f+0.05*rand(size(mean_changes,1),1),mean_changes(:,f),'o')
        hold on
        title(sprintf('sub-gamma: %1.1e, low-gamma: %1.1e, high-gamma: %1.1e\np = %1.3f',...
            mean(mean_changes(:,1)),mean(mean_changes(:,2)),mean(mean_changes(:,3)),pval))
    end
end
end

%% What electrodes are involved
all_rho_and_p = nan(length(rise_pts),2);
all_powers = cell(length(rise_pts),2);
% Loop over rise pts
for i = 1:length(rise_pts)
    p = rise_pts(i);

    
    % Get power in spike periods across channels
    sp_power = spike_power(p).spike_powers;
    
    % Get pre-spike abs power data
    data = abs_power.freq(1).pt(p).sp_or_not(1).data;
    
    % before rise indicators
    before_rise = pre_spike(p).before_rise;

    % Loop over spikes
    nsp = size(data,1);
    nch = size(data,3);
    sp_rhos = nan(nsp,1);
    sp_change = nan(nsp,nch);
    sp_rel_change = nan(nsp,nch);
    for s = 1:nsp
        curr_sp_power = squeeze(sp_power(s,:));
        curr_sp_pre_data = squeeze(data(s,:,:));
        curr_before_rise = squeeze(before_rise(s,:));
        
        last_pre_rise = find(curr_before_rise == 0);
        last_pre_rise = last_pre_rise(1) - 1; % move one back to get last pre rise
        
        curr_sp_first = curr_sp_pre_data(1,:);
        curr_sp_last = curr_sp_pre_data(last_pre_rise,:);
        change = (curr_sp_last - curr_sp_first);
        rel_change = (curr_sp_last - curr_sp_first)./abs(curr_sp_first);
        
        curr_rho = corr(curr_sp_power',change');
        sp_rhos(s) = curr_rho;
        sp_change(s,:) = change;
        sp_rel_change(s,:) = rel_change;
    end
    
    % Median change
    median_sp_change = median(sp_change,1);
    mean_sp_rel_change = mean(sp_rel_change,1);
    median_sp_power = median(sp_power,1);
    
    [rho,pval] = corr(median_sp_power',median_sp_change');
    all_rho_and_p(i,:) = [rho,pval];
    all_powers{i,1} = median_sp_power';
    all_powers{i,2} = median_sp_change';
    
end

if 0
figure
set(gcf,'position',[100 100 1000 1000]);
[ha,~] = tight_subplot(length(rise_pts),2,[0.03 0.02],[0.01 0.05],[0.03 0.04]);
for i = 1:length(rise_pts)
    for j = 1:2
        axes(ha((i-1)*2+j))
        scatter3(all_locs{i}(:,1),all_locs{i}(:,2),all_locs{i}(:,3),...
            100,all_powers{i,j},'filled');
        hold on
        scatter3(all_locs{i}(:,1),all_locs{i}(:,2),all_locs{i}(:,3),...
            100,'k');
        if i == 1
            view([-236.4833   -1.0483]);
            if j == 1
                title('Median IED power (\muV^2)','fontweight','normal');
            elseif j ==2
                title('Median pre-IED power change (\muV^2)',...
                    'fontweight','normal');
            end
        elseif i == 2
            view([111.4028  -10.3355]);
        elseif i == 3
            view([-158.7083  -16.3805]);
        elseif i == 4
            view([-49.0094  -31.9073]);
        elseif i == 5
            view([-33.3085   -1.2509]);
        elseif i == 6
            view([-172.7596  -12.7132]);
        end
        
        if j == 1
            pnum = all_rho_and_p(i,2);
            if pnum < 0.001
                ptext = '<0.001';
            else
                ptext = sprintf('%1.3f',pnum);
            end
            text(-0.05,0.85,sprintf('r = %1.2f\np = %s',...
                all_rho_and_p(i,1),ptext),...
                'Units','normalized','fontsize',20)
        end
        colorbar
        set(gca,'fontsize',20)
        axis('off')
    end
end
print(gcf,[out_folder,'electrode_change'],'-dpng');
end

%% Is there a corresponding network change?
if 1
% Loop over rise pts
all_data = cell(length(rise_pts),1);
for i = 1:length(rise_pts)
    p = rise_pts(i);
    
    % Get power in spike periods across channels
    sp_power = spike_power(p).spike_powers;
    [~,max_ch] = max(sp_power,[],2);
    
    % before rise indicators
    before_rise = pre_spike(p).before_rise;
    nsp = size(sp_power,1);
    
    first_last = nan(nsp,nfreq,2);
    
    for f = 1:1%nfreq
    
        % Get node strength
        data = pearson.freq(f).pt(p).sp_or_not(1).data;

        % get the data on the channel with the highest spike power change
        data_main_ch = get_data_from_spec_ch(data,max_ch);
    
        nsp = size(data_main_ch,1);
        
        
        % Loop over spikes
        for s = 1:nsp
            data_sp = squeeze(data_main_ch(s,:));
            %data_sp = squeeze(data(s,:,:));
            before_rise_sp = squeeze(before_rise(s,:));
            last_time = find(before_rise_sp==0);
            last_time = last_time(1)-1; % 1 before is the last before rise
            
            
            %first = data_sp(1,:);
            %last = data_sp(last_time,:);
            %first_last(s,f,:) = [mean(first) mean(last)];
            
            first = data_sp(1);
            last = data_sp(last_time);
            first_last(s,f,:) = [first last];
        end
        
    end
    
    all_data{i} = first_last;
    
    % Plot and significance testing on single pt level
    if 0
        figure
        set(gcf,'position',[440 490 1001 308])
        for f = 1:nfreq
            [~,pval] = ttest(squeeze(first_last(:,f,2)),squeeze(first_last(:,f,1)));
            
            plot(f+0.05*rand(size(first_last,1),1),...
                first_last(:,f,2)-first_last(:,f,1),'o') 
            hold on
            text(f,max(first_last(:,f,2)-first_last(:,f,1)),...
                sprintf('%1.3f',pval),'fontsize',20)
        end
        plot(xlim,[0 0],'k')
        pause
        close(gcf)
    end
    
end

figure
for i = 1:length(rise_pts)
    dat = squeeze(all_data{i}(:,1,1:2));
    plot(i+0.1*rand(size(dat,1),1),dat(:,2)-dat(:,1),'ko','markersize',10)  
    hold on
end
xlim([0 length(rise_pts)+1])
plot(xlim,[0 0],'k--')
xticklabels([])
xlabel('Patient')
ylabel({'Pre-IED - baseline','node strength difference'})

if 1
yl = ylim;
ylim([yl(1) yl(1) + (yl(2)-yl(1))*1.1])
yl = ylim;
for i = 1:length(rise_pts)
    dat = squeeze(all_data{i}(:,1,1:2));
    [~,pval] = ttest(dat(:,1),dat(:,2));
    if pval >= 0.05
        text(i,yl(1) + 0.9*(yl(2)-yl(1)),'n.s.','fontsize',20,...
            'horizontalalignment','center');
    end
end
end

set(gca,'fontsize',20)
print(gcf,[out_folder,'pearson'],'-dpng')

end




end


