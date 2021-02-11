function sanity_checks(metrics,met,pt_rise)

if 0
%% Plot an example spike and not spike
s = 30;
p = 1;
f = 1;
times = metrics.time.freq(f).(met).pt(p).times;
data_sp = metrics.time.freq(f).(met).pt(p).spike.data(s,:);
data_not = metrics.time.freq(f).(met).pt(p).not.data(s,:);
if size(pt_rise(p).both,1) ~= size(metrics.time.freq(f).(met).pt(p).spike.data,1)
    error('what')
end
rise_time = min(pt_rise(p).both(s,:));
figure
subplot(2,1,1)
plot(times,data_sp)
hold on
plot([rise_time rise_time],get(gca,'ylim'),'k--')
subplot(2,1,2)
plot(times,data_not)
end

if 0
%% Plot average across spikes for example pt
p = 1;
f = 1;
times = metrics.time.freq(f).(met).pt(p).times;
data_sp = mean(metrics.time.freq(f).(met).pt(p).spike.data,1);
data_not = mean(metrics.time.freq(f).(met).pt(p).not.data,1);
figure
subplot(2,1,1)
plot(times,data_sp)
subplot(2,1,2)
plot(times,data_not)
end

if 0
%% Relative change for example spike
s = 2;
p = 1;
f = 1;
times = metrics.time.freq(f).(met).pt(p).times;
data_sp = metrics.time.freq(f).(met).pt(p).spike.data(s,:);
data_not = metrics.time.freq(f).(met).pt(p).not.data(s,:);
if size(pt_rise(p).both,1) ~= size(metrics.time.freq(f).(met).pt(p).spike.data,1)
    error('what')
end
rise_time = min(pt_rise(p).both(s,:));
figure
subplot(2,1,1)
plot(times,(data_sp-data_sp(1))/abs(data_sp(1)))
hold on
plot([rise_time rise_time],get(gca,'ylim'),'k--')
subplot(2,1,2)
plot(times,(data_not-data_not(1))/abs(data_not(1)))
end

if 0
%% Relative change, median across spikes, example pt
for p = 1:15
f = 1;
times = metrics.time.freq(f).(met).pt(p).times;
data_sp = (metrics.time.freq(f).(met).pt(p).spike.data);
data_not = (metrics.time.freq(f).(met).pt(p).not.data);

data_sp = (data_sp-data_sp(:,1))./abs(data_sp(:,1));
data_not = (data_not-data_not(:,1))./abs(data_not(:,1));

data_sp = median(data_sp,1);
data_not = median(data_not,1);

figure
subplot(2,1,1)
plot(times,data_sp)
xlim([-2 -0.2])
subplot(2,1,2)
plot(times,data_not)
xlim([-2 -0.2])
pause
close(gcf)
end
    
end

%% Relative change, median across spikes, all pts
if 0
np = 15; 
nt = 21;
all_data_sp = zeros(np,nt);
all_data_not = zeros(np,nt);
for p = 1:np
    f = 1;
    times = metrics.time.freq(f).(met).pt(p).times;
    data_sp = (metrics.time.freq(f).(met).pt(p).spike.data);
    data_not = (metrics.time.freq(f).(met).pt(p).not.data);

    data_sp = (data_sp-data_sp(:,1))./abs(data_sp(:,1));
    data_not = (data_not-data_not(:,1))./abs(data_not(:,1));

    data_sp = median(data_sp,1);
    data_not = median(data_not,1);
    
    all_data_sp(p,:) = data_sp;
    all_data_not(p,:) = data_not;

end
figure
subplot(2,1,1)
errorbar(times,mean(all_data_sp,1),std(all_data_sp,[],1))
xlim([-2 -0.2])
subplot(2,1,2)
errorbar(times,mean(all_data_not,1),std(all_data_not,[],1))
xlim([-2 -0.2])
end

end
