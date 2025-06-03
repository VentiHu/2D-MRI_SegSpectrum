clear
close all

load('..\Data_Spearman.mat')
load('..\D_T2_cluster_result.mat')
% CompartmentNames could be loaded from D_T2_cluster.mat, or defined here
%CompartmentNames = ["S", "R", "L", "PD"];

disp('Activity Index')
f1 = figure;
hold on
for j = 1:4
    mask = D_T2_cluster(:,:,j);
    pvalMatrix = pvalMatrix_AI;
    pval_list = pvalMatrix(mask == 1);
    
    disp(CompartmentNames{j})
    report_statistics(pval_list);
    %figure
    histogram(pval_list,'BinEdges', 0:0.05:1.0,'FaceAlpha', 0.5)
end
hold off
legend(CompartmentNames)
legend('Location', 'best');
xlabel('P-value with Activity Index'); 
ylabel('Number of meshpoint');
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'LineWidth', 1);
saveas(f1,'histogram_corrAI.png')

disp('Chronicity Index')
f2 = figure;
hold on
for j = 1:4
    mask = D_T2_cluster(:,:,j);
    pvalMatrix = pvalMatrix_CI;
    pval_list = pvalMatrix(mask == 1);

    disp(CompartmentNames{j})
    report_statistics(pval_list);
    
    histogram(pval_list,'BinEdges', 0:0.05:1.0,'FaceAlpha', 0.5)
end
hold off
legend(CompartmentNames)
legend('Location', 'best');
xlabel('P-value with Chronicity Index'); 
ylabel('Number of meshpoint');
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'LineWidth', 1);
saveas(f2,'histogram_corrCI.png')

function report_statistics(data)
    p05 = prctile(data, 5);
    disp(['05 percentile: ', num2str(p05)]);
    p25 = prctile(data, 25);
    disp(['25 percentile: ', num2str(p25)]);
    p50 = prctile(data, 50);
    disp(['50 percentile: ', num2str(p50)]);
    p75 = prctile(data, 75);
    disp(['75 percentile: ', num2str(p75)]);
    p95 = prctile(data, 95);
    disp(['95 percentile: ', num2str(p95)]);
end