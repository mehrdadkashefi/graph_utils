function scatter_SEM_1(D,y, F1, tick_labels)
%SCATTER_SEM Summary of this function goes here
%   Detailed explanation goes here
color_set = [27,158,119; 217,95,2; 117,112,179]/255;
alpha_val = 0.2;
rand_spread = 0.03;
n_F1 = unique(D.(F1));

group_mean = zeros(length(n_F1), 1);
group_sem = zeros(length(n_F1), 1);

for f1 =1:length(n_F1)
    idx = D.(F1) == n_F1(f1);
    group_mean(f1) = mean(D.(y)(idx));
    group_sem(f1) = std(D.(y)(idx));
end
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(group_mean);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    scatter(x, group_mean(:,i),'Marker',"_", 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor',color_set(i,:), 'LineWidth',2);
end

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    num_scatter_sample = length(D.(y));
    rand_val = rand_spread*random('Uniform',-1,1, num_scatter_sample,1);
    [C,~,ic] = unique(D.(F1));
    scatter(x(ic)+rand_val', D.(y), 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor',color_set(i,:), 'MarkerFaceAlpha',alpha_val, 'MarkerEdgeAlpha',alpha_val);
    errorbar(x, group_mean(:,i), group_sem(:,i), 'k', 'linestyle', 'none', 'LineWidth',1);
    scatter(x, group_mean(:,i),'Marker',"_", 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor','black', 'LineWidth',2);

end
xticks(x)
xticklabels(tick_labels)
end

