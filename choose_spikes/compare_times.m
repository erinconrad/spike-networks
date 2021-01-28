function compare_times(erin,jim)

times1 = [];
times2 = [];

for i = 1:length(erin.spike)
    times1 = [times1;erin.spike(i).time];
    times2 = [times2;jim.spike(i).time];
    
end

histogram(times1,20)
hold on
histogram(times2,20)
legend('Erin','Jim')

end