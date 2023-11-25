function scatter_SEM(D,y, F1, F2, F2_legend, varargin)
%SCATTER_SEM Summary of this function goes here
%   Detailed explanation goes here
 
try 
    color_set = varargin{1};
catch
    color_set = [27,158,119; 217,95,2; 117,112,179;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100;100 100 100]/255;
end

try 
    alpha_val = varargin{2};
catch
    alpha_val = 0.4;
end

try 
    rand_spread = varargin{3};
catch
    rand_spread = 0.03;
end
try
    plot_line = varargin{4};
catch
    plot_line = false;
end



[n_F1, ~, val_idx_F1] = unique(D.(F1));
[n_F2, ~, val_idx_F2] = unique(D.(F2));


group_mean = zeros(length(n_F1), length(n_F2));
group_sem = zeros(length(n_F1), length(n_F2));

for f1 =1:length(n_F1)
    for f2 =1:length(n_F2)
        idx = D.(F1) == n_F1(f1) & D.(F2) == n_F2(f2);
        group_mean(f1, f2) = mean(D.(y)(idx));
        group_sem(f1, f2) = std(D.(y)(idx));
    end
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
legend(F2_legend,'AutoUpdate','off','Location','best')
if plot_line
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        plot(x, group_mean(:, i), 'LineWidth',3, 'Color',color_set(i,:));
    end
end
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    num_scatter_sample = length(D.(y)(D.(F2) == n_F2(i)));
    rand_val = rand_spread*random('Uniform',-1,1, num_scatter_sample,1)';
    %[C,~,ic] = unique(D.(F1)(D.(F2) == n_F2(i)));
    scatter(x(val_idx_F1(D.(F2) == n_F2(i)))+rand_val, D.(y)(D.(F2) == n_F2(i)), 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor',color_set(i,:), 'MarkerFaceAlpha',alpha_val, 'MarkerEdgeAlpha',alpha_val);
    errorbar(x, group_mean(:,i), group_sem(:,i), 'k', 'linestyle', 'none', 'LineWidth',1);
    scatter(x, group_mean(:,i),'Marker',"_", 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor','black', 'LineWidth',3);

end
xticks([1:ngroups])
xticklabels(num2str(n_F1))

end

