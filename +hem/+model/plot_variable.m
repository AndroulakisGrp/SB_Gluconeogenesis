function plot_variable(t, y, data_t, data_y, t_pre, line_style, color)
% t: Time vector from the model
% y: Variable from the model
% data_t: Time vector for real data (Optional)
% data_y: Real data (Optional)
% t_pre: Where to start plotting t (Optional; default is 48)

if nargin<=2
    data_t = [];
end
if nargin<=3
    data_y = [];
end
if nargin<=4
    t_pre = 48;
end
if nargin<=5
    line_style = '';
end
if nargin<=6
    color = 'k';
end

hold on;

plot(t(t>t_pre)-t_pre, y(t>t_pre), sprintf('%s%s', color, line_style), 'LineWidth', 1)
% plot(t(t>t_pre)-t_pre, y(t>t_pre), sprintf('%s', line_style), 'Color',color,'LineWidth', 1)
% hold on
% line([10 10],[0 10]);
if ~isempty(data_t) && ~isempty(data_y)
%     plot(data_t, data_y, sprintf('%s.', color), 'MarkerSize', 20);
plot(data_t, data_y,'.', 'Color',color, 'MarkerSize', 20);
end

hold off;

% y-axis scaling
% ys = [data_y, y(t>t_pre)'];
% y0 = min(ys);
% y1 = max(ys);
% ylim([y0-0.1*(y1-y0), y1+0.1*(y1-y0)]);
