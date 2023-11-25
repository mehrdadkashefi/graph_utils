function pointplot(D,y, F1, F2, F2_legend, varargin)
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
try
    error_type = varargin{5};
catch
    error_type = 'ci';
end


[n_F1, ~, val_idx_F1] = unique(D.(F1));
[n_F2, ~, val_idx_F2] = unique(D.(F2));
 

group_mean = zeros(length(n_F1), length(n_F2));
group_error = zeros(length(n_F1), length(n_F2));

for f1 =1:length(n_F1)
    for f2 =1:length(n_F2)
        idx = D.(F1) == n_F1(f1) & D.(F2) == n_F2(f2);
        group_mean(f1, f2) = mean(D.(y)(idx));
        if strcmp(error_type, 'sd')
            group_error(f1, f2) = std(D.(y)(idx));
        elseif strcmp(error_type, 'se')
            n_value = length(D.(y)(idx));
            group_error(f1, f2) = std(D.(y)(idx))/sqrt(n_value);
        elseif strcmp(error_type, 'ci')
            n_value = length(D.(y)(idx));
            group_error(f1, f2) =  1.96 * std(D.(y)(idx))/sqrt(n_value);
        end
    end
end
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(group_mean);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    scatter(x, group_mean(:,i),'Marker',"o", 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor',color_set(i,:), 'LineWidth',2);
end
legend(F2_legend,'AutoUpdate','off','Location','best')
if plot_line
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        plot(x, group_mean(:, i), 'LineWidth',4, 'Color',color_set(i,:));
    end
end
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, group_mean(:,i), group_error(:,i),'Color', color_set(i,:), 'linestyle', 'none', 'LineWidth',1);
    scatter(x, group_mean(:,i),100 ,'Marker',"o", 'MarkerFaceColor', color_set(i,:), 'MarkerEdgeColor',color_set(i,:), 'LineWidth',3);

end
xticks([1:ngroups])
xticklabels(num2str(n_F1))

end

