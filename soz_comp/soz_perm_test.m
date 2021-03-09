clear

%% Parameters
nboot = 1e4;

%% Get data for soz 
[~,all,total_soz_n,all_main,all_soz_chs] = soz_info;
big = all(:,2);
lead = all(:,1);
main_big = all_main(:,2);
main_lead = all_main(:,1);
nsoz = total_soz_n(:,1);
nelecs = total_soz_n(:,2);
np = size(all,1);


%% For how many patients is the main spiking electrode in the SOZ?
n_big_in_soz = sum(big);
n_lead_in_soz = sum(lead); % same for lead electrode

%% Monte carlo
boot_successes = nan(nboot,1);
% Loop through MC iterations
for ib = 1:nboot
    %{
    if mod(ib,100) == 0
        fprintf('\nDoing iteration %d of %d\n',ib,nboot);
    end
    %}
        
        
    successes = nan(np,1);
    
    % Loop through patients
    for ip = 1:np
    
        % roll an n_elec sided die
        die_roll = randi(nelecs(ip));
        
        % count it as a success if it's less than or equal to the number of
        % soz electrode. This is equivalent to saying, under null
        % distribution, there is an nsoz/nelecs chance of the highest
        % spiking electrode being in the soz (assuming under the null
        % distribution that the highest spikeing electrode is chosen
        % randomly amongst all electrodes)
        if die_roll <= nsoz(ip)
            successes(ip) = 1;
        else
            successes(ip) = 0;
        end
        
    end
    
    n_successes = sum(successes == 1);
    boot_successes(ib) = n_successes;
    
end

%% Compare the actual number to the Monte Carlo number
boot_successes = sort(boot_successes);

p_big = (1+sum(boot_successes >= n_big_in_soz))/(nboot+1);
p_lead = (1+sum(boot_successes >= n_lead_in_soz))/(nboot+1);

fprintf('\nThe big and the lead electrode were the same for %d patients.\n',...
    sum(main_big == main_lead));
fprintf(['\nFor %d patients, the big electrode was a SOZ.\n'...
    'For %d patients, the lead electrode was a SOZ.\n'],...
    n_big_in_soz, n_lead_in_soz);
fprintf(['\nThe monte carlo probability is %1.3f for the big spiking electrode\n'...
    'and %1.3f for the lead electrode.\n'],...
    p_big,p_lead);

%% Plot
if 0
    figure
    plot(boot_successes,'o')
    hold on
    plot(xlim,[n_big_in_soz n_big_in_soz])
    plot(ylim,[n_lead_in_soz n_lead_in_soz])
    legend('Monte carlo','Big','Lead');
    
    
end

