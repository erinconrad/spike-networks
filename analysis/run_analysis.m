clear

%% Parameters
metrics = {'abs_power','ers','ns'};
time_window = 0.1;

for m = 1:length(metrics)

    met = metrics{m};
    
    %% Get data and do basic processes
    % Get data
    orig = get_metric_data(met,time_window);

    % Find windows for each spike in which no early spike rise
    pre_spike = find_pre_spike_windows(orig.times);

    % We only need this once
    if m == 1
        % Get clinical info
        clinical = get_clinical_info(orig);
    end
    
    % Get the spike time powers
    spike_power = get_spike_time_powers(orig);
    
    % Get times of spikes
    timing = get_spike_timing;

    % Find bad spikes
    bad = identify_bad_spikes;

    % Remove bad spikes
    [clean,clean_spike,clean_pre,clean_timing] = remove_bad(bad,orig,spike_power,pre_spike,timing);

    all_metrics(m).name = met;
    all_metrics(m).data = clean;
    all_metrics(m).spike = clean_spike;
    all_metrics(m).pre = clean_pre;
    all_metrics(m).timing = clean_timing;
    
    % Get names
    all_names = {};
    for i = 1:length(orig.freq(1).pt)
        all_names = [all_names;orig.freq(1).pt(i).name];
    end
    

end


%% Find the patients with a significant pre-spike rise in absolute power
%{
This checks:
- Is there an overall pre-spike rise in absolute power in aggregate
- How many patients have this rise?
- Do these patients have a rise in the IED-free control periods?
- Clinical predictor of who has a rise?
- How many spikes have a rise?
%}
data = all_metrics(1).data.freq(1);
npts = length(data.pt);
[rise_pts,change] = count_sig_power_rise(data,all_metrics(1).spike,all_metrics(1).pre,clinical);

%% Visualize highest pre-rise spikes
if 0
p = 14;
[~,I] = sort(change{p},'descend');
show_specified_spike(p,I(1:10),all_names,bad)
end


%% Does something predict whether individual spikes have a pre-spike power rise?
%{
This is where I test:
- Does spike location (which electrode the spike is on) predict pre-spike
rise?
- More specifically, if the spike is in the SOZ, does that predict a rise?- NEED TO ADD *********
- Does spike timing predict pre-spike rise?
%}
rel_change = analyze_spikes_with_rise(data,all_metrics(1).spike,all_metrics(1).pre,all_metrics(1).timing,1:npts);

%% For spikes with a pre-spike rise, test additional features
%{
This is where I test, for spikes with a pre-spike rise:
- When does it begin?
- What frequencies does it involve?
- What electrodes does it involve?
- Is there a change in network connectivity?

%}

detailed_sp_analyses(rise_pts,all_metrics(1).spike,all_metrics(1).pre,...
    all_metrics(1).data,all_metrics(2).data,all_metrics(3).data);

%{
detailed_sp_analyses(rel_change,rise_pts,all_metrics(1).spike,all_metrics(1).pre,...
    all_metrics(1).data,all_metrics(2).data,all_metrics(3).data);
%}

%% Do analyses
if 0
    
    % Get the number of patients with a pre-spike rise and an aggregate
% statistic
count_sig_power_rise(clean,clean_spike,clean_pre,clinical,1,alpha)

    
    
    
% Plot the metric over time for spikes and non spikes
show_metric_over_time(clean,clean_spike,clean_pre,1,1)




analyze_spikes_with_rise(clean,clean_spike,clean_pre,clean_timing,1,1)


metric_rise_by_ch(clean,clean_spike,clean_pre,1,1,alpha)

agg_metric_rise_by_ch(clean,clean_spike,clean_pre,1,alpha)


end

