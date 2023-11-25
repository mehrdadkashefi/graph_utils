function bar_SEM(model_series,model_error)
%BAR_SEM Summary of this function goes here
%   Detailed explanation goes here
color_set = [215,25,28; 253,174,97; 255,255,191; 171,217,233; 44,123,182]/255;

b = bar(model_series, 'grouped');
hold on


if size(model_series,2)<=length(color_set)
    for i=1:length(b)
        b(i).FaceColor = color_set(i,:);
    end
else
    
end


% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(model_series);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

end
