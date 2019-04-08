%{

Plan for spike networks project:
- the goal is to determine whether interictal spikes acutely change nodal,
local (close to the spike), or global network properties.
- I will examine the second prior to spikes, the second after, and then a
separate second removed from spikes. I will compare various network
measures in those 3 seconds. I will aggregate spikes within a patient and
do within-patient analyses. I will then aggregate some statistic across
patients to do a between-patient analysis.

Storage requirements:
- number of bytes = 
(N_spikes*2) * N_ch * N_s * N_freq * precision * N_patients 
= 1e3 * 2 * 100 * 10* 512 * 8 * 30 = 245 GB


Specific steps:
1) Identify spike and not-a-spike-times
    - for spike times, I want to be pretty sure I am getting real spikes,
    and so I will only take those spikes that are post-processing from the
    other spike projects. 
    - for not-a-spike, I want to be pretty sure they are not spikes. So I
    will take the pre-artifact removal spike times. For each one of these,
    I will look 5 seconds after the spike. If it is not within 5 seconds
    of another "spike", I will include it as "not-a-spike".
    - I will randomly take N spikes and N not-a-spikes, where N = 1,000
2) Download chunks of ieeg data
    - for the T seconds surrounding each spike and "not a spike", I will
    download the data. T = 10 seconds, centered around the spike or 
    not-a-spike.
3) Identify times that are before the spike and after the spike
    - Get the first spike in the sequence, find the peak, and then decide
    some T_begin that is long enough before the spike peak.
    - Get the last spike in the sequence, and decide some T_end that is
    long enough after the spike peak.
    - Take one second ending at the threshold start time and one second
    beginning at the threshold end time.
    - Do it again for another threshold time, maybe a second later.
3) Calculate adjacency matrices
    - For the seconds I decided, I will calculate adjacency matrices for 
    the major frequency bands (~6 matrices).
    - I will save the adjacency matrices in structures, 10 spikes at a
    time.
    - the size of each file will be 10 * 4 * 100 * 100 * 6 * 8 = 19 MB and
    there will be 200 files per patient
4) Calculate various network metrics
    - For each adjacency matrix, get various network metrics
        - global metrics: synchronizability, transitivity, global
        efficiency
        - nodal metrics of the electrodes involved in the spikes: node
        strength, CC, ???
5) Compare the network metrics in the pre-, post-, and non-spike times
    - Within a patient, could combine all of them and do a Friedman or
    repeated measures ANOVA
    - Between patients, I could combine the statistics and compare with a
    Friedman test
    


%}