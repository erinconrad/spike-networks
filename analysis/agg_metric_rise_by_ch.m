function agg_metric_rise_by_ch(main,spike,clean_pre,f,alpha)

sp = 1; % spike only

np = length(spike);

all_pre = cell(np,1);
all_sp = cell(np,1);
rhos = nan(np,1);
zs = nan(np,1);

for p = 1:np

    % Get electrode spike power
    sp_power = spike(p).spike_powers;
    mean_sp_power = mean(sp_power,1);

    % Get pre-spike data
    data = main.freq(f).pt(p).sp_or_not(sp).data;

    % initialize
    first = nan(size(data,1),size(data,3));
    last = nan(size(data,1),size(data,3));

    % Loop over spikes
    for s = 1:size(data,1)
        if sp == 1
            current_pre = clean_pre(p).before_rise(s,:);
        else
            current_pre = clean_pre(p).fake_before_rise(s,:);
        end
        last_pre = find(current_pre ==0);
        last_pre = last_pre(1) - 1; % the last pre-spike window

        first(s,:) = squeeze(data(s,1,:));
        last(s,:) = squeeze(data(s,last_pre,:));
    end


    % plot last-first for all electrodes, see if larger change in electrodes
    % with greater spike power
    change = last-first;
    mean_change = mean(change,1);
    
    normalized_mean_change = (mean_change-mean(mean_change))/std(mean_change);
    normalized_mean_sp_power = (mean_sp_power-mean(mean_sp_power))/std(mean_sp_power);
    
    all_pre{p} = normalized_mean_change;
    all_sp{p} = normalized_mean_sp_power;
    rhos(p) = corr(mean_sp_power',mean_change');
    zs(p) = atanh(corr(mean_sp_power',mean_change'));
end


mean_z = mean(zs);
mean_r = tanh(mean_z);

figure
set(gcf,'position',[78 500 1330 305])
for p = 1:np
    plot(p,rhos(p),'o')
    %plot(all_sp{p},all_pre{p},'o')
    hold on
end

end