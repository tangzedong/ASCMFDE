clear 
close all
ker = 1/5*ones(1,5);
titlename = {'CI+LS', 'PI+LS', 'NI+LS'};
nn = 1;
for t = [12,15,18]
    figure
    filename = sprintf('successrateASCMFDE%d-1.mat', t);
    load(filename);
    b = zeros(30,1000);
    for r = 1:30
        b(r,:) = data_result(r).crsuccessrate;
    end
    a1 = mean(b,1);
    std1 = 1/2*std(b,1);
    
    successrate1 = conv(a1,ker,'same');
    successrate1(1:5) = a1(1:5);
   
    lower = successrate1-std1;
    upper = successrate1 + std1;
    fill([1:1000,1000:-1:1],[lower, upper(end:-1:1)],[1, 0.8, 0.8],'facealpha',0.5);
    hold on
    plot(successrate1+std1, 'r-', 'LineWidth', 0.5,'Color', [1, 0.5, 0.5]);%,'MarkerIndices',1:50:length(successrate));
    hold on
    plot(successrate1-std1, 'r-', 'LineWidth', 0.5,'Color', [1, 0.5, 0.5]);%,'MarkerIndices',1:50:length(successrate));
    hold on
    
    filename = sprintf('successrateMFDE%d-1.mat', t);
    load(filename);
    b = zeros(30,1000);
    for r = 1:reps
        b(r,:) = data_result(r).crsuccessrate;
    end
    a2 = mean(b,1);
    std2 = 1/2*std(b,1);
    successrate2 = conv(a2,ker,'same');
    successrate2(1:5) = a2(1:5);
    lower = successrate2-std2;
    upper = successrate2 + std2;
    fill([1:1000,1000:-1:1],[lower, upper(end:-1:1)],[0.8, 0.8, 1],'facealpha',0.5);
    hold on
    plot(successrate2+std2, '-', 'LineWidth', 0.5, 'Color', [0.5, 0.5, 1]);%,'MarkerIndices',1:50:length(successrate));
    hold on
    plot(successrate2-std2, '-', 'LineWidth', 0.5, 'Color', [0.5, 0.5, 1]);%,'MarkerIndices',1:50:length(successrate));
    hold on
    
    
    b1 = plot(successrate1, 'r-', 'LineWidth', 2,'Color', [1, 0, 0]);%,'MarkerIndices',1:50:length(successrate));
    hold on
    b2 = plot(successrate2, '-', 'LineWidth', 2, 'Color', [0, 0, 1]);%,'MarkerIndices',1:50:length(successrate));
    hold on
    legend([b1,b2],{'ASConly','MFDE'});
    xlabel('Generation', 'FontName', 'Times New Roman', 'FontSize', 20);
    ylabel('Transferring Success Rate', 'FontName', 'Times New Roman', 'FontSize', 20);
    yl = ylim;
    ylim([0, yl(2)]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
    title(titlename{nn});
    grid on
    print(sprintf('ascfigureSucRate%d',t),'-depsc2');
    nn = nn + 1;
end