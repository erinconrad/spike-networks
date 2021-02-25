duration = 3; % 3 seconds
f = 60;
f1 = 60; 
f2 = 90;
f3 = 120;
fs = 512;
times = 0:1/fs:3;


if 0 
    
    Y = sin(2*pi*f*times);
    
    % Plot the signal
    figure
    plot(times,Y);

    % Plot the spectrogram
    figure
    spectrogram(Y,128,120,128,fs,'yaxis');

    % Get power in diff freq bands
    freq_bands = [59 61;61 63];
    powers = get_power(Y',fs,freq_bands);
end

% Make multichannel signal
Y = zeros(length(times),2);
Y(:,1) = 0.01*randn(length(times),1)+(sin(2*pi*f1*times)+sin(2*pi*f2*times))';
Y(:,2) = 0.01*randn(length(times),1)+(sin(2*pi*f1*times)+sin(2*pi*f3*times))';

% Plot the spectrogram
figure
subplot(2,1,1)
spectrogram(Y(:,1),128,120,128,fs,'yaxis');
set(gca,'fontsize',20)
subplot(2,1,2)
spectrogram(Y(:,2),128,120,128,fs,'yaxis');
set(gca,'fontsize',20)

% Ms cohere
freq_bands = [59 61; 89 91;119 121];
adj = get_adj_matrices(Y,fs,freq_bands);

range = [0 0];
for i = 1:3
    curr_max = max(max(adj(i).adj));
    curr_min = min(min(adj(i).adj));
    if curr_max > range(2)
        range(2) = curr_max;
    end
    
    if curr_min < range(1)
        range(1) = curr_min;
    end
end

% Plot the adj
figure
for i = 1:3
subplot(1,3,i)
imagesc(adj(i).adj)
caxis(range)
colorbar
title(sprintf('Coherence from %d-%d Hz',freq_bands(i,1),freq_bands(i,2)));
end

