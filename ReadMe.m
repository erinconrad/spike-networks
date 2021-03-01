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
- choose_spikes/denote_bad_spikes: this script plots spike data; I used
this to double check spikes, looking for spikes that did not appear to be
true spikes or spikes that had other spikes in the 2 seconds surrounding
the spike that I had missed on ieeg. It outputs a file called "bad" with
the identities of the "bad" spikes.
- classifier/remove_bad_spikes: this takes the bad spikes I denoted and
removes them from further analysis
- choose_spikes/manually_find_rise: this script also plot spikes data and
allows the reviewer to manually select the earliest spike rise time. It
outputs a file called "early" with the earliest spike rise times for each
spike. Two independent board-certified epileptologists reviewed spikes
using this script, resulting in two sets of rise times for each spike.
- classifier/multi_reviewer_pre_spike: this script takes the "early"
structure containing the manually identified spike rise times from each of
the two reviewers. It loops through the spike and find the windows for
which the full time of the window occurs before the spike rise time for
both reviewers.
-classifier/include_which_times: this takes the time windows from the
multi_reviewer_pre_spike info and makes an array for each spike of which
time windows to include for analysis
- classifier/get_specified_metrics: ugly script that pulls info for power,
ers, node strength and puts them all into a uniform structure for further
analysis


4) Get spike powers
- abs_power/get_abs_power: this script gets the absolute power for manually
detected spikes and saves it in a structure called sig_dev
- get_adj_matrices/ers_spike: this script gets frequency-specific powers
for manually detected spikes and saves them in a structure called ers
- math/get_power: this is called by ers_spike above and does the actual
calculation of frequency-specific power by taking the fft of the signal to
get it into the frequency domain, squaring the abs value, and then taking
just the component in the right frequency
- classifier/get_sd: just pulls the absolute power data and puts it into a
structure
- classifier/convert_sd: restructuring of the absolute power data

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

6) Measure metric changes
- classifier/generate_summary_stats: this is the main analysis script. It
takes the metrics obtained above, and the pre-spike time windows to analyze
(from include_which_times) and it calculates the relative metric change and
does significance testing.

7) Make figures
- plots/methods_fig_2: Plots example spikes, a spike broken up into time
windows, the metric change for spike vs not a spike, and the metric change
at different times, and then the frequency specific power changes
- plots/methods_fig_3: Plots a single spike with time windows, a bunch of
adjacency matrices, the node strength for the biggest channel, and the
average node strength
- plots/methods_fig_4: plots an example spike sequence, how often the lead
and peak IED are in the SOZ, the power change for IEDs in SOZ vs not, and
the power change for lead vs non lead channels

%}