function compare_spike_function(whichPts)

for whichPt = whichPts
    
    % Get functional network node strength
    [ns_fun,locs_fun] = pca_2(whichPt);
    
    % Get spike coactivation network node strength
    [ns_sp,locs_sp] = spike_coactivation_old(whichPt);
    
    % Throw an error if locs aren't the same
    if isequal(locs_fun,locs_sp) == 0
        error('what\n');
    end
    
    locs = locs_fun;
    
    figure
    subplot(1,2,1)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,ns_fun,'filled');
    title('Functional');
    
    subplot(1,2,2)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,ns_sp,'filled');
    title('Spike');
    
    [rho,pval] = corr(ns_fun',ns_sp','Type','Spearman');
    fprintf('Patient %d rho = %1.2f, p = %1.3f\n',whichPt,rho,pval);
    
    
end



end