function bar_SEM(model_series,model_error,varargin)
%BAR_SEM Summary of this function goes here
%   Detailed explanation goes here

% default color set:
color_set = [215,25,28; 253,174,97; 255,255,191; 171,217,233; 44,123,182]/255;
% default xtick labels:
default_xtick_labels = [];   
% plot_type:
default_type = 'dashplot';
expected_types = {'dashplot','barplot'};
% line styles:
default_line_groups = 1;


% input parser:
p = inputParser;
addOptional(p,'xtick_labels',default_xtick_labels,@(x) length(x)==size(model_series,2));
addOptional(p,'type',default_type);
addOptional(p,'line_groups',default_line_groups)
parse(p,model_series,model_error,varargin{:});

% b = bar(model_series, 'grouped');
% hold on

% setting colors:
C = zeros(size(model_series,2),3);
if size(model_series,2)<=length(color_set)
    for i=1:size(model_series,2)
        C(i,:) = color_set(i,:);
    end
else

end


% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(model_series);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% create the bars:
if strcmp(p.Results.type,'dashplot')
    for i = 1:ngroups
        for j = 1:nbars
            % dash center:
            x = i - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
            width_dash = (-groupwidth/2 + (2*j-1) * groupwidth / (2*nbars))/2;
            x_fill = [x-width_dash , x+width_dash , x+width_dash, x-width_dash];
            y_fill = [model_series(i,j)-0.01 , model_series(i,j)-0.01 , ...
                 model_series(i,j)+0.01 , model_series(i,j)+0.01];
            fill(x_fill,y_fill,C(j,:),'FaceAlpha',0.6,'LineStyle','none')
            hold on
            % scatter(x,model_series(i,j),10,C(j,:),'filled')
        end
    end
    xticks(1:ngroups)
    
    if p.Results.line_groups
        linestyle = ':';
    else
        linestyle = 'none';
    end

    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', linestyle);
    end
    
    ylim([0,max(max(model_series))*120/100])

    xticklabels(p.Results.xtick_labels)
    box off

end

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
% end


