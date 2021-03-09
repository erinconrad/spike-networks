%function spike_comp

%{
why is df for ns_avg 8
%}

%% Clear
clear

%% Parameters
windows = 0.1; % can be 0.1 or 0.2
do_plot = 0;

% Don't change these
do_auto = 1; % should be 1
do_cumulative = 0; % should be zero
rm_rise = 1; % should be 1


%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));

bct_folder = locations.BCT;
addpath(genpath(bct_folder));
out_folder = [results_folder,'plots/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


for met_all = {'sd','ers','ns_big','ns_avg'}
met = met_all{1};

% this just says that I am using the automatically calculated peak IED
% channel for doing calculations of metric change
if do_auto
    if strcmp(met,'ns_big')
        met = 'ns_auto';
    elseif strcmp(met,'ns_avg')
        met = 'ns_avg';
    else
        met = [met,'_auto'];
    end
end

%% Find windows for each spike in which no early spike rise
pre_spike = multi_reviewer_pre_spike(windows);

%% Compare two reviewers
[earliest_rise,pt_rise] = compare_two_reviewers(pre_spike);
orig_pt_rise = pt_rise;

%% Get signal power deviation
sig_dev = get_sd;

% convert this to be similar to pre_spike
pre_spike = convert_sd(sig_dev,windows,pre_spike,met);

%% Get network metrics
[metrics,is_spike_soz,is_spike_depth] = get_specified_metrics(windows,pre_spike,met);



%% Remove bad spikes
[metrics,is_spike_soz,pre_spike,pt_rise,is_spike_depth] = ...
    remove_bad_spikes(metrics,is_spike_soz,pre_spike,pt_rise,is_spike_depth);

%% Remove the time windows with an early spike rise, get slopes, and do significance testing
include_times = include_which_times(metrics,met,pre_spike,nan);
metrics = generate_summary_stats(metrics,met,include_times,rm_rise,is_spike_soz,...
    do_cumulative,is_spike_depth,earliest_rise);

%% Build a classifier to predict spike vs not spike
%tbl = class_spike(metrics,met);

%% Sanity checks
% sanity_checks(metrics,met,pt_rise)

if strcmp(met,'sd_auto')
    metrics_sd = metrics;
elseif strcmp(met,'ers_auto')
    metrics_ers = metrics;
elseif strcmp(met,'ns_auto')
    metrics_ns_big = metrics;
elseif strcmp(met,'ns_avg')
    metrics_ns_avg = metrics;
end

%% Print some summary stats
if contains(met,'sd')
    nfreq = 1;
else
    nfreq = 3;
end
for f = 1:nfreq
    fprintf(['\n\nFor %s %s, the relative metric for IEDs was (M = %1.1e, SD = %1.1e)'...
        ' and for non-IEDs was (M = %1.1e, SD = %1.1e), t(%d) = %1.1f, p = %1.3f\n'],...
        met,metrics.time.freq(f).name,...
        mean(metrics.time.freq(f).(met).auc.data(:,1)),std(metrics.time.freq(f).(met).auc.data(:,1)),...
        mean(metrics.time.freq(f).(met).auc.data(:,2)),std(metrics.time.freq(f).(met).auc.data(:,2)),...
        metrics.time.freq(f).(met).auc.df,metrics.time.freq(f).(met).auc.tstat,...
        metrics.time.freq(f).(met).auc.pval);
    
    if metrics.time.freq(f).(met).auc.pval < 0.05/nfreq
        fprintf('\nThe earliest pre-IED change is %1.1f\n',metrics.time.freq(f).(met).auc.short.first_sig_time);
    end
    
    %{
    fprintf(['\nThe relative metric for SOZ was (M = %1.1e, SD = %1.1e)'...
        ' and for non-SOZ was (M = %1.1e, SD = %1.1e), t(%d) = %1.1f, p = %1.3f\n'],...
        nanmean(metrics.time.freq(f).(met).auc.soz.data(:,1)),nanstd(metrics.time.freq(f).(met).auc.soz.data(:,1)),...
        nanmean(metrics.time.freq(f).(met).auc.soz.data(:,2)),nanstd(metrics.time.freq(f).(met).auc.soz.data(:,2)),...
        metrics.time.freq(f).(met).auc.soz.df,metrics.time.freq(f).(met).auc.soz.tstat,...
        metrics.time.freq(f).(met).auc.soz.pval); 
   
    
    fprintf(['\nThe relative metric for lead IED was (M = %1.1e, SD = %1.1e)'...
        ' and for other sequence IED was (M = %1.1e, SD = %1.1e), t(%d) = %1.1f, p = %1.3f\n'],...
        nanmean(metrics.time.freq(f).(met).auc.first_v_other.data(:,1)),nanstd(metrics.time.freq(f).(met).auc.first_v_other.data(:,1)),...
        nanmean(metrics.time.freq(f).(met).auc.first_v_other.data(:,2)),nanstd(metrics.time.freq(f).(met).auc.first_v_other.data(:,2)),...
        metrics.time.freq(f).(met).auc.first_v_other.df,metrics.time.freq(f).(met).auc.first_v_other.tstat,...
        metrics.time.freq(f).(met).auc.first_v_other.pval); 

        %}
end
end

%% Make supplemental tables
meta_metrics.metrics_sd = metrics_sd;
meta_metrics.metrics_ers = metrics_ers;
meta_metrics.metrics_ns_big = metrics_ns_big;
meta_metrics.metrics_ns_avg = metrics_ns_avg;
do_save = 0;
make_supplemental_tables(meta_metrics,results_folder,do_save)

%% Figs
%methods_fig_2(metrics_sd,metrics_ers,earliest_rise,orig_pt_rise)
%methods_fig_3(metrics_ns_big,metrics_ns_avg,earliest_rise,orig_pt_rise)
%methods_fig_4(metrics_sd)

%plot_auc(metrics,met,out_folder,do_plot);
%plot_short_both(metrics,met,2,earliest_rise,out_folder,do_plot,rm_rise);
%soz_comparison(metrics,met,out_folder)
%count_sig_pts(metrics,met)

%plot_short(metrics,met,1,earliest_rise,out_folder,do_plot);




%end