function [pc, score, latent, tsquare, mu, sigma] = pca(t_plot, y_plot, vars)
% Riboleukogram-like PCA analysis

concatenate_variables = 0; % Set to 0 to average variables together, like in riboleukogram paper
if concatenate_variables
    x = zeros(length(t_plot)*M, length(vars));
else
    x = zeros(length(t_plot), length(vars));
end

M = size(y_plot.P, 1); % Number of profiles

for i=1:M
    % x(k,j) = average over all simulations, the jth variable at the kth time point
    % or
    % x(k,j) = all simulations concatenated, the jth variable at the kth (concatenated) time point
    for j=1:length(vars)
        if concatenate_variables
            a = (i-1)*length(t_plot) + 1;
            b = a + 999;
            x(a:b, j) = eval(sprintf('y_plot.%s(i,:)', vars{j}))';
        else
            x(:,j) = mean(eval(sprintf('y_plot.%s', vars{j})))';
        end
    end
end
mu = mean(x);
sigma = std(x);
zx = bsxfun(@minus, x, mu);
zx = bsxfun(@rdivide, zx, sigma);

[pc, score, latent, tsquare] = princomp(zx);