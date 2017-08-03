function [t, y] = run(varargin)


%% Function parameters
p = inputParser;
p.StructExpand = true;


% INTEGRATION
% Enforce nonnegativity?
p.addParamValue('nonnegative', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% How long to integrate?
p.addParamValue('t_final', 48, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
% How long to run before starting everything? This is helpful if parameters
% or equations are changing, so initial conditions might also change.
p.addParamValue('t_pre', 48, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
% false: deterministic, one cell; true: stochastic, multiple cells
% Time range of exogenous LPS dosing
p.addParamValue('t_ex_LPS', [0, 0], @(x)validateattributes(x, {'numeric'}, {'vector', 'size', [1, 2]}));
% Stochastic?
p.addParamValue('stochastic', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));


% MODEL
% Include circadian components?
p.addParamValue('circadian', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Use LPS logistic growth term?
p.addParamValue('transition_logistic', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Use HPA model of cortisol?
p.addParamValue('HPA', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Load parameters from files if they are not passed explicitly
p.addParamValue('hh_circ', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_fen', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_all', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('new', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_jusko', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_sig', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('hh_opt', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('met', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('glu', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
% Period of scotoperiod
p.addParamValue('sp', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));
% Light stimulus minimum and maximum
p.addParamValue('lst1', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('lst2', 0, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('lst3', 6, @(x)validateattributes(x, {'numeric'}, {'vector'}));


% PLOTTING
% Show any plots?
p.addParamValue('make_plots', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Show data (for baseline conditions) on main plot?
p.addParamValue('plot_data', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% If make_plots is also true, show some additional plots (currently 2: all
% of the variables in subplots and individual vs average values of NFkBn)
p.addParamValue('extra_plots', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
p.addParamValue('ss_plots', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
% Color for the main plot
p.addParamValue('color', 'k', @(x)validateattributes(x, {'char'}, {}));
% Line style for the main plot
p.addParamValue('line_style', '-', @(x)validateattributes(x, {'char'}, {}));

p.addParamValue('runx', 0, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
p.addParamValue('runy', 0, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
p.addParamValue('fdpd', [], @(x)validateattributes(x, {'numeric'}, {'vector'}));

% Food stimulus minimum and maximum
p.addParamValue('fdhigh', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('fdlow', 0, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('fdphase', 6, @(x)validateattributes(x, {'numeric'}, {'vector'}));

p.addParamValue('cells', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('count', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));
p.addParamValue('paramtest', 1, @(x)validateattributes(x, {'numeric'}, {'vector'}));

p.parse(varargin{:});
s = p.Results;

s.t_final = s.t_final + s.t_pre;

%% Parameter and initial condition defaults

f = hem.util.get_hem_folder; % Path to +hem folder

% Initial conditions - should be recalculated if parameters or equations

y0_original = importdata(sprintf('%s/param/ICs.txt', f));
y0_met = importdata(sprintf('%s/param/ICmet.txt', f));
y0=[y0_original' y0_met']';




% First, load all default parameter values
params_to_test = {'hh_circ', 'hh', 'hh_fen', 'hh_all', 'new', 'met', 'hh_jusko', 'hh_sig', 'hh_opt','glu'};
load_param = zeros(size(params_to_test)); % Was parameter passed to function or loaded from file?
for i=1:length(params_to_test)
    if isempty(s.(params_to_test{i}))
        s.(params_to_test{i}) = importdata(sprintf('%s/param/%s.txt', f, params_to_test{i}));
        load_param(i) = true;
    else
        load_param(i) = false;
    end
end

% Is circadian enabled or not?
if ~s.circadian
    s.hh_circ = zeros(size(s.hh_circ));
end


%% Integration


options = odeset('Nonnegative', 1:27);

if ~s.stochastic
    %     keyboard;
    %tspan1  = [0 s.t_final];
    tspan1 = [0 1000];
   
    [t1, y1] = ode45(@hem.model.model2, tspan1, y0 , options, s);
    [tc, yc] = ode45(@hem.model.modelc, tspan1, y0 , options, s);    
    t = t1;
    y = y1;
    
    assignin('base', 't1', t1)
    
    
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc'};
    ids=1:length(vars);
   
    y_old = y;
    y_oldc = yc;
    y = struct;
    yc = struct;
    for j=1:length(ids)
        id = ids(j);
        y.(vars{j}) = y_old(:,id)';
        yc.(vars{j}) = y_oldc(:,id)';
    end
    
else
    % Stochastic constant noise integration
    y0=y0';
    mu = 0.03*ones(size(y0));
    M = length(y0);                    % num of variables
    dt = 0.1; N = s.t_final/dt;
    dW = sqrt(dt)*randn(N, M);         % Brownian increments
    y_temp = y0;
    y = [];
    y(1,:) = y_temp;
    %     keyboard;
    for j = 1:N
        t = dt*j;
        dy = hem.model.model(t, y_temp,s);
        y_temp = y_temp + dt*dy + mu.*dW(j,:);
        y(j+1,:) = y_temp;
        %        keyboard;
    end
    t = 0:dt:s.t_final;
    
    % keyboard;
    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc'};
    ids=1:length(vars);
    y_old = y;
    y = struct;
    for j=1:length(ids)
        id = ids(j);
        y.(vars{j}) = y_old(:,id)';
    end
end
% keyboard;



%% Plotting
if s.make_plots
    vars_to_plot = {'F','PERCRY','CLOCKBMAL','PGC1aN','Gluc'};
    titles = {'Cortisol','PER/CRY','CLOCK/BMAL1','PGC1-\alpha(N)','Pck1/G6pc'};

    
    n = length(vars_to_plot);
    
    %     Make appropriately sized subplot
    nx = ceil(sqrt(n));
    ny = ceil(sqrt(n));
    if nx*(ny-1) >= n
        ny = ny - 1;
    end
    
%     time=linspace(0,1000,10000);
%     lightreg=zeros(1,length(time));  lightb=zeros(1,length(time));  lightd=zeros(1,length(time));
%     for i=1:length(time)
%         if (mod(time(i)-6,24) < 12)
%             lightreg(i)=1;
%         else
%             lightreg(i)=0;
%         end
%     end
%     for i=1:length(time)
%         lightb(i)=1;
%         lightd(i)=0.1;
%     end

    for i=1:n
        var = vars_to_plot{i};
        hold on
        subplot(2, 3,i+1);
        hold on
        plot(tc, mean(yc.(var),1),'k','LineWidth',2);
        hold on
        plot(t, mean(y.(var), 1),'b--','LineWidth',2);
        %plot(t, mean(y.(var), 1),'r-.','LineWidth',2);
        
        xlim([840 840+24]);
        set(gca, 'XTick', 840:12:840+24);
        set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
        hold on;
        title(titles(i),'FontSize',18);
        ax=gca;
        ax.XAxis.FontSize=14;   ax.YAxis.FontSize=14;
    end
    
    %Restricted Feeding
    
    time=linspace(0,1000,10000);
    feedcontrol=zeros(1,length(time));
    feedres=zeros(1,length(time));
    for kk=1:length(time)
        if (mod(time(kk)-6,24)<12)
            feedcontrol(kk)=1;
        else feedcontrol(kk)=0;
        end
        if time(kk)<480
            if (mod(time(kk)-6,24) < 12)
                feedres(kk)=1;
            else
                feedres(kk)=0;
            end
        else
            if (mod(time(kk)-s.fdphase,24) < s.fdpd)
                feedres(kk)=s.fdhigh;
            else
                feedres(kk)=s.fdlow;
            end
        end
    end
    subplot(2,3,1)
    hold on
    plot(time,feedres,'b--','LineWidth',2)
    hold on
    plot(time,feedcontrol,'k:','LineWidth',2)
    xlim([480 480+24]);
    set(gca, 'XTick', 480:12:480+24);
    ylim([0 2.2])
    title('Feed','FontSize',18);
    set(gca, 'XTickLabel',{'12am', '12pm', '12am'})

end

if s.extra_plots
    % Plot all variables along with their max and min values

    vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
        'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
        'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
        'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
        'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc'};
    n = length(vars);
    
    figure
    hold on
    for i=1:n
        var = vars{i};
        subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
        plot(t, mean(y.(var), 1),'b','LineWidth',2);
        hold on
        plot(tc, mean(yc.(var),1),'k:','LineWidth',2);
        xlim([840 840+24]);
        set(gca, 'XTick', 840:12:840+24);
        set(gca, 'XTickLabel',{'12am', '12pm', '12am'})
        hold on;
        title(var,'FontSize',10);
    end
    text(0.5, 1,'\bf Reversed Feeding','HorizontalAlignment',...
        'center','VerticalAlignment', 'top')

end
%%
if s.ss_plots
     vars = {'F', 'mRNA_R_F', 'R_F', 'FR' 'FRn', 'CRH', 'ACTH', 'F_per',...
         'percry', 'PERCRY', 'nucPERCRY', 'bmal', 'BMAL', 'nucBMAL',...
         'CLOCKBMAL','R_F_per', 'FR_per', 'FRn_per', 'M_F_per', 'FM_per',...
         'FMn_per', 'mRNA_P', 'P', 'mRNA_R_P','R_P', 'PR','P_cen',...
         'NAD','NAM','NMN','SIRT1','CLOCKBMALSIRT1','NAMPT','feed2','feed3','EntF',...
        'pgc1a','PGC1a','PGC1aN','PGC1aNa','FOXO1','Gluc'};

n = length(vars);
    figure
    hold on
    for i=1:n
        var = vars{i};
        subplot(ceil(sqrt(n)), ceil(sqrt(n)), i);
        plot(t, mean(y.(var), 1),'k-');

        xlim([400 700]);

        title(vars(i),'FontSize',10);
    end

end
end
