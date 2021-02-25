%{
This is the codebase for the project looking at pre- and post-spike
network changes. The following is the pipeline for the project:

1) Manually detect spikes
- These are stored in a google sheet in the spike_networks folder. I only
included spikes that didn't have any other spikes in the 2 seconds before
or 2 seconds after. 

2) Get EEG data for spikes
- choose_spikes/spike_processing: this script gets the EEG data surrounding
 manually detected spikes, and stores it in an output structure. It calls
 get_manual_spike_times_from_excel below to get the spike times.
- choose_spikes/get_manual_times_from_excel: this script takes manual
spikes from an excel file and outputs a structure with spike times
- choose_spikes/plot_manual_spike_examples: this script will take the EEG
data produced by the above functions and plot example spikes

4) Various spike prep things:

4) Get spike powers
- abs_power/get_abs_power: this script gets the absolute power for manually
detected spikes and saves it in a structure called sig_dev
- get_adj_matrices/ers_spike: this script gets frequency-specific powers
for manually detected spikes and saves them in a structure called ers
- math/get_power: this is called by ers_spike above and does the actual
calculation of frequency-specific power by taking the fft of the signal to
get it into the frequency domain, squaring the abs value, and then taking
just the component in the right frequency

5) Get spike network changes
- get_adj_matrices/manual_spike_networks: this script calculates the
adjacency matrices associated with the EEG data surrounding the manually
detected spikes. It calls the functions below to do the actual math in
calculating these networks.
- get_adj_matrices/get_adj_matrices: this is called by the above script to
calculate functional coherence across a range of frequency bands and across
all channel pairs to generate an adjacency matrix
- get_adj_matrices/pre_processing: this is called by manual_spike_networks
to do common average referencing and notch filtering
- network_metrics/get_ns: this takes adjacency matrices produced by the
above scripts and calculates node strength



%}