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


%}