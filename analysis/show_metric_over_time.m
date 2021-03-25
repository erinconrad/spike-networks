function show_metric_over_time(main,spike,clean_pre,f,p)

%% Parameters
only_pre = 1;

%% Get max spike power electrode
sp_power = spike(p).spike_powers;
[~,max_ch] = max(sp_power,[],2);

%% Get fake max spike power electrode for non-spike periods
n_not_spikes = size(main.freq(f).pt(p).sp_or_not(2).data,1);
% Randomly sample n_not_spikes from max chs
fake_max_ch = datasample(max_ch,n_not_spikes);

%% Get fake pre-times for non-spike periods
for s = 1:n_not_spikes
    clean_pre(p).fake_before_rise(s,:) = mode(clean_pre(p).before_rise,1);
end

figure
set(gcf,'position',[173 334 1268 464])
%% Get info for plotting spikes
for sp = 1:2
    data = main.freq(f).pt(p).sp_or_not(sp).data;
    main_ch_data = nan(size(data,1),size(data,2));
    first = nan(size(data,1),1);
    last = nan(size(data,1),1);
    
    % Loop over spikes
    for s = 1:size(data,1)

        if sp == 1
            % Get the power in the biggest spike power electrode
            main_ch_data(s,:) = data(s,:,max_ch(s));
            current_pre = clean_pre(p).before_rise(s,:);
        else
            % Get the power in the biggest fake power electrode
            main_ch_data(s,:) = data(s,:,fake_max_ch(s));
            current_pre = clean_pre(p).fake_before_rise(s,:);
        end

        % remove non-pre times
        last_pre = find(current_pre ==0);
        last_pre = last_pre(1) - 1; % the last pre-spike window
        if only_pre
            main_ch_data(s,current_pre==0) = nan;
        end

        first(s) = main_ch_data(s,1);
        last(s) = main_ch_data(s,last_pre);

    end


    
    subplot(2,2,(sp-1)*2+1)
    % Histogram of first vs last time
    
    plot(1+0.05*rand(length(first),1),first,'o')
    hold on
    plot(2+0.05*rand(length(first),1),last,'o')
    xticks([1 2])
    xticklabels({'First time','Last time'})
    [~,pval] = ttest(first,last);
    title(sprintf('p = %1.3f',pval))
    xlim([0 3])
    yticklabels([])
    if sp == 1
        ylabel('IED period')
    else
        ylabel('IED-free period')
    end
    set(gca,'fontsize',20)
    


    subplot(2,2,(sp-1)*2+2)
    % Plot
    if 1
    for t = 1:size(main_ch_data,2)
        curr_time = main_ch_data(:,t);
        errorbar(main.times(t),nanmean(curr_time),nanstd(curr_time),'o');
        hold on
    end
    end
    xlabel('Time relative to IED peak (s)')
    yticklabels([])
    set(gca,'fontsize',20)
    
end





end