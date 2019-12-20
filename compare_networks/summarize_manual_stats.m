function summarize_manual_stats(simple)

%% parameters
alpha = 0.05;

%% Get file locations, load spike times and pt structure
locations = spike_network_files;
main_folder = locations.main_folder;
results_folder = [main_folder,'results/'];
data_folder = [main_folder,'data/'];
script_folder = locations.script_folder;
addpath(genpath(script_folder));
pt_file = [data_folder,'spike_structures/pt.mat'];
bct_folder = locations.BCT;
addpath(genpath(bct_folder));
plot_folder = [results_folder,'plots/'];

if simple == 1
    out_folder = [results_folder,'perm_stats/simple/'];
elseif simple == 0
    out_folder = [results_folder,'perm_stats/coherence/'];
end

% Get the individual patient stat files
listing = dir([out_folder,'*_perm.mat']);

% Initialize table of significant values
sig_table = cell2table(cell(0,5),'VariableNames',{'Name','Freq','Time','F','p'});

% Loop through patients
for i = 1:length(listing)
    
    % Load the file
    name = listing(i).name;
    sim = load([out_folder,name]);
    sim = sim.sim;
    
    nfreq = length(sim);
    
    % Loop through frequencies
    for f = 1:nfreq
        
        % Loop through time points
        for t = 2:length(sim(f).p)
            
            % multiple comparisons (correct for number of time comparisons
            % and number of frequency comparisons)
            % 
            % Note I don't think this can ever be significant when testing
            % multiple frequencies ********
            if sim(f).p(t) < alpha/(nfreq)/(length(sim(f).p)-1)
                
                % Add info to the significance table
                sig_table = [sig_table;{name,f,t,sim(f).F(t),sim(f).p(t)}];
                
            end
            
        end
        
    end
    
    
end

% Store the significance table
save([out_folder,'sig_table.mat'],'sig_table')

% Display the significance table
sig_table

end