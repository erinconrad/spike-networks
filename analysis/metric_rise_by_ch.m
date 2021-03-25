function metric_rise_by_ch(main,spike,clean_pre,p,f,alpha)

thresh_change = 10;
nchs = size(main.freq(f).pt(p).sp_or_not(1).data,3);

% Get electrode spike power
sp_power = spike(p).spike_powers;
mean_sp_power = mean(sp_power,1);

% For each spike, what electrode has the max spike power
[~,max_ch] = max(sp_power,[],2);

% Get fake pre-times for non-spike periods
n_not_spikes = size(main.freq(f).pt(p).sp_or_not(2).data,1);
for s = 1:n_not_spikes
    clean_pre(p).fake_before_rise(s,:) = mode(clean_pre(p).before_rise,1);
end

%
sp = 1;
data = main.freq(f).pt(p).sp_or_not(sp).data;
first = nan(size(data,1),size(data,3));
last = nan(size(data,1),size(data,3));
max_ch_change = nan(size(data,1),1);

% Loop over spikes
for s = 1:size(data,1)
    current_pre = clean_pre(p).before_rise(s,:);
    
    last_pre = find(current_pre ==0);
    last_pre = last_pre(1) - 1; % the last pre-spike window

    first(s,:) = squeeze(data(s,1,:));
    last(s,:) = squeeze(data(s,last_pre,:));
    
    % Change in electrode with max power
    max_ch_change(s) = (last(s,max_ch(s)) - first(s,max_ch(s)))/abs(first(s,max_ch(s)));
end

change = (last-first);

% Find spike with large relative pre-spike power change
big_change = find(max_ch_change>thresh_change); 

%{
% Find the spikes with the largest change in 
[~,top_10_change] = sort(max_ch_change,'descend');
top_10_change = top_10_change(1:10);
%}

% Correlate biggest sp electrodes with biggest pre-spike electrodes on a
% spike by spike basis
rho = corr(change',sp_power'); % pairwise correlation between each pair of spikes
rho = diag(rho); % the diagonal is each spike with each corresponding spike

rho(big_change)
mean(rho(big_change))

%}

%{
figure
set(gcf,'position',[78 500 1330 305])

for sp = 1:2

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

    % Show mean across all spikes, show all channels, first vs last
    if 0
    plot(1+0.05*rand(size(first,2),1),mean(first,1),'o')
    hold on
    plot(2+0.05*rand(size(first,2),1),mean(last,1),'o')
    xlim([0 3])
    end

    pvals = nan(nchs,1);
    for ich = 1:nchs
        [~,ptemp] = ttest(first(:,ich),last(:,ich));
        pvals(ich) = ptemp;
    end

    % plot last-first for all electrodes, see if larger change in electrodes
    % with greater spike power
    change = last-first;
    %change = (last-first)./abs(first);
    mean_change = mean(change,1);
    subplot(1,2,sp)
    plot(mean_sp_power,mean_change,'o')
    [rho,pval] = corr(mean_sp_power',mean_change');
    title(sprintf('r = %1.2f, p = %1.3f',rho,pval))
end
%}
end