function spike_stats(whichPts)

%% Parameters
% 1 = alpha/theta; 2 = beta, 3 = low gamma, 4 = high gamma, 5 = ultra high, 6 = broadband
freq_text = {'alpha/theta','beta','low\ngamma','high\ngamma','ultra high\ngamma','broadband'};
%freq_text = {'alpha/theta'};
t_text = {'t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11'};
n_f = length(freq_text);
n_times = 11;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
spike_times_file = [data_folder,'spike_times/times.mat'];
pt_file = [data_folder,'spike_structures/pt.mat'];

times = load(spike_times_file); % will result in a structure called "out"
times = times.out;
pt = load(pt_file); % will create a structure called "pt"
pt = pt.pt;

if isempty(whichPts) == 1
    whichPts = 1:length(times);
end

for whichPt = whichPts
    
    %% Prep patient
    
    % Skip it if name is empty
    if isempty(times(whichPt).name) == 1, continue; end
    name = times(whichPt).name;
    fprintf('\nDoing %s\n',name);
    
    pt_folder = [results_folder,name,'/'];
    stats_folder = [pt_folder,'stats/'];
    
    % Load the stats file
    out = load([stats_folder,'stats.mat']);
    stats = out.out;
    
    % Get variables
    z = stats.signal.z;
    dev = stats.signal.dev;
    bin_z = stats.signal.bin_z;
    bin_dev = stats.signal.bin_dev;
    ec = stats.network.ec;
    
    
    % Just look at alpha/theta
    z_ec = (((ec-mean(ec,3))./std(ec,0,3)));
    z_ec_at = squeeze(z_ec(1,:,:));
    
    
    %% Compare EC at spike time between spike chs and non spike chs
    
    
    %% Signal deviation, z score
    % Do a repeated measures ANOVA for signal deviation
    s_num = 1:size(bin_z,1);
    c = {};
    for i = 1:size(bin_z,2)
        c{i} = bin_z(:,i);
    end
    t_dev = table(s_num',c{[1:5,7:11]},...
        'VariableNames',{'spike',t_text{[1:5,7:11]}});
    
    rm_dev = fitrm(t_dev,'t1-t10 ~ spike');
    ranovatbl = ranova(rm_dev)
    
    % Try a friedman test instead
    [p,tbl,stats] = friedman(bin_z(:,[1:5,7:11]),1,'off');
    
    %% Eigenvector centrality
    % Do a repeated measures ANOVA for ec
    c = {};
    for i = 1:size(z_ec_at,2)
        c{i} = z_ec_at(:,i);
    end
    t_ec = table(s_num',c{[1:5,7:11]},...
        'VariableNames',{'spike',t_text{[1:5,7:11]}});
    rm_ec = fitrm(t_ec,'t1-t10 ~ spike');
    ranovatbl = ranova(rm_ec)
    
    % Friedman
    [p,tbl,stats] = friedman(z_ec_at(:,[1:5,7:11]),1,'off');
    
    % Paired T-tests to compare
    [h,p,ci,stats] = ttest(z_ec_at(:,1),z_ec_at(:,5))
    [h,p,ci,stats] = ttest(z_ec_at(:,1),z_ec_at(:,4))
    [h,p,ci,stats] = ttest(z_ec_at(:,1),z_ec_at(:,3))
    
end

end