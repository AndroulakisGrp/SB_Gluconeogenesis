tic
close all;
clear all;
clc;

s.t_pre = 960;
s.t_final=100;
s.make_plots = true;
s.plot_data = true;


f = hem.util.get_hem_folder;
hh_store = importdata(sprintf('%s/param/hh_all.txt', f));


% lb=[0.05033,0.7163,0.03018,0.1573,0.01182,0.38,0.50111,0.01029,0.012802,0.008261,0.11205,0.00579,0.004135,0.003685,0.167228,0.037359,0.01156,0.001523,0.006875,0.002041,0.06702,0.009156,0.002918,0.011393,0.001932,0.0014,0.012688,0.00453,0.09,0.01,0.0056,0.0012,0.003,0.0005,0.0012,0.0024,0.0002,0.036,0.0216,0.0075,0.0024,0.0006,0.0012,0.0045,0.0006,0.0009,3e-05,0.0009,0.0001,0.08,0.02,0.03,1.9e-05,0.049605,0.5717,0.406033,1.042588,0.5,1,1,0.01];
% ub=[1.3,1.3,2.18,14.73,1.82,2,71.11,102.92,128.02,82.61,3.05,0.56,41.35,36.85,1672.28,373.59,115.6,15.23,68.75,20.410,670.2,91.56,29.18,113.93,19.32,14.01,126.88,45.3,900,100,56,12,30,5,12,24,2,360,216,75,24,6,12,45,6,9,0.3,9,1,800,200,300,0.19,496.05,5.17,1060.33,4.88,100,100,5,1];
lb=[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,5,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,.01];
ub=[10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,13,10,10,10,10,10,10,10,10,10,10,10,10,10];


v = [69,70,72,73,80,81]; % Parameters of the vector hh_circ(:,1) that want to change
vars = {'P'}; % Functions which parameters are goint to be optimized.


% options = optimset('LargeScale','on','UseParallel','always','FunvalCheck','on','Display','iter','MaxFunEvals',55000);

hh_change=hh_store;
[err,fval,exitflag,output] = fminsearch(@(hh_change) hem.model.obj_fun_any(hh_change,hh_store,vars,v, s),hh_change(v));
% opts = optimset('Algorithm','sqp');
% problem = createOptimProblem('fmincon','objective',...
%  @(hh_change) hem.model.obj_fun_any(hh_change,hh_store,vars,v, s),'x0',hh_change(v),'lb',lb(v),'ub',ub(v),'options',opts);
% ms = MultiStart;
% [x,fval,exitflag,output,solutions] = run(ms,problem,10);


%% Plot results to compare

s.make_plots = true;
% s.hh_opt = hh_store;
figure(1);
clf;
hem.model.run2(s);
% figure(2);
% clf;
% run_model(.1, 8, s);
% s.hh_circ = hh_circ;
% figure(2);
% clf;
% hem.model.run(0, 9, s);
% figure(4);
% clf;
% run_model(.1, 8, s);
% hh_circ
toc