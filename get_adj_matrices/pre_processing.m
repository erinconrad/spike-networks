function data = pre_processing(data,do_car,pre_whiten)

%% Common average reference
if do_car == 1
    avg = mean(data,2);
    old_data = data;
    data = data - avg;
end

%% Apply an AR(1) model
%{
The goal here is to fit an AR(1) model to data and then output the
residual, so that I am removing the autocorrelation component. I don't know
if I want to do this, because it tends to wash out the spikes.
%}
if pre_whiten == 1
    for j = 1:size(data,2)
        vals = data(:,j);
        mdl = regARIMA('ARLags',1);
        mdl = estimate(mdl,vals,'Display','Off');
        E = infer(mdl,vals);
        
        if 0
           figure
           plot(old_data(:,j))
           hold on
           plot(vals)
           hold on
           plot(E)
           legend('Original','CAR','whitened')
           pause
           close gcf
        end
        
        data(:,j) = E;
    end
end


end