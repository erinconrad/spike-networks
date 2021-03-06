function [pdc_out,dtf_out] = pdc(y,Fs,index_windows,freq_bands)

% Erin edited mvar to take out an error if the number of samples is too
% small relative to the number of electrodes...

L = size(y,1); % Number of samples
CH = size(y,2); % Number of channels
N_freq = round(Fs); % Number of frequency points
Fmax = Fs/2;      % Maximum frequency limit in the PDC and DTF plots

%% Time-invariant MVAR estimation -- to estimate the optimum moel order
[w, A_TI, C_TI, sbc, fpe, th] = arfit(y, 1, 20, 'sbc'); % ---> ARFIT toolbox
[tmp,p_opt] = min(sbc); % Optimum order for the MVAR model

[PDC_TI, DTF_TI] = PDC_DTF_matrix(A_TI,p_opt,Fs,Fmax,N_freq); % Compute time-varying PDC and DTF measures

%% Short-Time PDC and DTF of the original signal
Mode1 = 10; % ARFIT(10), Multichannel Yule-Walker(1) ---> BioSig toolbox 
s = 1;
for i = 1:size(index_windows,1)
    start_point = index_windows(i,1);
    end_point = index_windows(i,2);
    seg = y(start_point:end_point,:);
    seg = seg.*repmat(hamming(size(seg,1)),1,CH);
    [A_ST,RCF,C_ST] = mvar(seg, p_opt, Mode1); % ---> BioSig toolbox
    [PDC_ST(:,:,:,s), DTF_ST(:,:,:,s)] = PDC_DTF_matrix(A_ST,p_opt,Fs,Fmax,N_freq);
    s = s + 1;
    
end

finterp = 1:Fmax;

pdc_out = zeros(size(PDC_ST,1),size(PDC_ST,2),size(freq_bands,1),size(PDC_ST,4));
dtf_out = zeros(size(PDC_ST,1),size(PDC_ST,2),size(freq_bands,1),size(PDC_ST,4));
for f = 1:size(freq_bands,1)
    % DTF instead of pdc since this emphasizes sources rather than sinks
    pdc_out(:,:,f,:) = sum(PDC_ST(:,:,finterp>=freq_bands(f,1)& finterp<=freq_bands(f,2),:),3);
    dtf_out(:,:,f,:) = sum(DTF_ST(:,:,finterp>=freq_bands(f,1)& finterp<=freq_bands(f,2),:),3);
end

% pdc_out is nch x nch x nfreq x ntimes
% for a given frequency and time, the i,j component is the directed inflow
% from channel j to channel i, normalized by dividing by all the inflows to
% channel i. And so a large value indicates that j has a lot of inflow into
% i.

end