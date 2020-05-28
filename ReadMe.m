%{
This is the codebase for the project looking at pre- and post-spike
network changes. The following is the current working pipeline for the 
project:

1) Manually detect spikes
- These are stored in a google sheet in the spike_networks folder. I only
included spikes that didn't have any other spikes in the 3 seconds before
or 3 seconds after. 

2) Get EEG data for spikes
- choose_spikes/spike_processing: this script gets the EEG data surrounding
 manually detected spikes, and stores it in an output structure. It calls
 get_manual_spike_times_from_excel below to get the spike times.
- choose_spikes/get_manual_times_from_excel: this script takes manual
spikes from an excel file and outputs a structure with spike times
- choose_spikes/plot_manual_spike_examples: this script will take the EEG
data produced by the above functions and plot example spikes

3) Get adjacency matrices for spikes
- get_adj_matrices/manual_spike_networks: this script calculates the
adjacency matrices associated with the EEG data surrounding the manually
detected spikes. It calls the functions below to do the actual math in
calculating these networks.
- get_adj_matrices/get_adj_matrices: this is called by the above script to
calculate functional coherence across a range of frequency bands and across
all channel pairs to generate an adjacency matrix
- get_adj_matrices/get_simple_corr: this is called by the
manual_spike_networks script to calculate simple pair-wise channel by
channel correlations to generate an adjacency matrix.
- get_adj_matrices/pre_processing: this is called by manual_spike_networks
to do common average referencing and/or pre-whitening by removing the AR(1)
component of the signal.
- get_adj_matrices/plot_manual_networks: this takes the output adjacency
matrix and plot the average adjacency matrices

4) Do various cross-time network comparisons
- compare_networks/manual_sig_deviation: this script takes the EEG data
surrounding the spikes and finds time periods in which there is a
significant deviation from the baseline. This is essentially the control
for the other comparisons.
- network_metrics/manual_network_metrics: this script takes the adjacency
matrices produced in step 4 and compares a couple of specific network
metrics (e.g., node strength of the involved channels and global efficiency
of the full network) across time points
- compare_networks/manual_permanova: this script takes the adjacency
matrices produced in step 4 and does a permanova to test for any changes in
the adjacency matrix across time (using a correction for multiple
comparisons)
- compare_networks/manual_nbs: NEED TO MAKE. the goal will be to use
network based statistics to compare networks across times in a way that
takes into account the network structure, hopefully boosting the power to
find a difference relative to a permanova.

5) Visualize results of these comparisons
- compare_networks/summarize_manual_stats: This makes summary tables about
the statistics produced from step 4. It is only configured to do the simple
correlation case right now.

%}