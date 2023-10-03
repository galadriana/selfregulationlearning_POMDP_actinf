%% CASE 1 DISRUPTION
%% GENERATE PLOT TO VISUALIZE THE ANTICORRELATION BETWEEN POLICIES


clear all
close all


% Path

cd 'C:\Users\gadri\spm12\toolbox\DEM\GVG_and_tutorial_pomdp\Supplementary code\GVG_DEM\hierarchical\Thesis\Model_1\reproducibility_test\case1b_1test'

% LOAD

% Performance
load('test_reproducibility_100rnf.mat')

% Cases
load('index_subebaja_1b.mat')


% DIRICHlet 
load('1b_learn_cosines_matrices_disruption_down.mat')
load('1b_learn_cosines_matrices_disruption_up.mat')


% MDP_Q
load('CQ_allDOWN_mean_table_reproducibility_100rnf.mat')
load('CQ_allDOWN_std_table_reproducibility_100rnf.mat')
load('CQ_allUP_mean_table_reproducibility_100rnf.mat')
load('CQ_allUP_std_table_reproducibility_100rnf.mat')

CQ_allUP_mean_table = CQ_allUP_mean_table(index_subebaja,:);
CQ_allDOWN_mean_table = CQ_allDOWN_mean_table(index_subebaja,:);
CQ_allUP_std_table = CQ_allUP_std_table(index_subebaja,:);
CQ_allDOWN_std_table = CQ_allDOWN_std_table(index_subebaja,:);

%EFE
load('EFEpolup_table_reproducibility_100rnf.mat');
load('EFEpoldown_table_reproducibility_100rnf.mat')

EFE_up_mean = mean(EFEpolup_table(index_subebaja,:));
EFE_down_mean = mean(EFEpoldown_table(index_subebaja,:));
EFE_up_std = std(EFEpolup_table(index_subebaja,:));
EFE_down_std = std(EFEpoldown_table(index_subebaja,:));

% VFE

load('F_table_reproducibility_100rnf.mat')
VFE_up_mean = mean(F_table(index_subebaja,:));
VFE_up_std = std(F_table(index_subebaja,:));


%% Color set
color = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0 0 0], [0 0 0 0.5], [0 0.4470 0.7410 0.5], [0.8500 0.3250 0.0980 0.8], [0.4660 0.6740 0.1880], [0.954 0.345 0.543], [0.4 0.6 0.1 0.5], [0.954 0.345 0.543 0.5], [0.4660 0.6740 0.1880 0.5], [0.92 0.69 0.1125], [0.6350 0.0780 0.1840], [0.2081 0.1663 0.5292] }; % Cell array of colros.

 
%% VISUALIZATION UP
mdpq_range = int16(linspace(1,600,600));
t = size(cosines_matrices_up_only,1);

figure()

subplot(4,1,1)
plot(mean(thermometer_table(index_subebaja,1:t)',2), LineWidth=2., Color=[0 0 0]);
%xlabel('Time step');
ylabel('Thermometer position');
title(' Thermometer behavior')

subplot(4,1,2)
s1EFE = shadedErrorBar(mdpq_range(1:t), EFE_up_mean(1:t), EFE_up_std(1:t), 'lineprops', '-blue','patchSaturation',0.085);
s1EFE.mainLine.Color = color{1};
s1EFE.mainLine.LineWidth = 1.5;
hold on
s2EFE=shadedErrorBar(mdpq_range(1:t), EFE_down_mean(1:t), EFE_down_std(1:t), 'lineprops', '-red','patchSaturation',0.085);
s2EFE.mainLine.Color = color{2};
s2EFE.mainLine.LineWidth = 1.5;
%xlabel('Time step');
ylabel('EFE');
legend('Policies causing up ' , 'Policies causing down', Location='northeast', fontsize=12);
title('Expected Free Eenergy')

subplot(4,1,3)
cosup=plot(mean(cosines_matrices_up_only,2)); % 'ColorLimits',[0 1]
cosup.Color = color{1};
cosup.LineWidth=2;
hold on
cosdown=plot(mean(cosines_matrices_down_only,2)); % 'ColorLimits',[0 1]
cosdown.Color= color{2};
cosdown.LineWidth=2;
title('Transition matrix memory formation');
%xlabel('Time (belief update)')
ylabel('Similarity with learning')

subplot(4,1,4)
s1 =shadedErrorBar(mdpq_range(1:t), CQ_allUP_mean_table(1,1:t), CQ_allUP_std_table(1,1:t), 'lineprops', '-blue','patchSaturation',0.085);
s1.mainLine.Color = color{1};
s1.mainLine.LineWidth = 1.5;
hold on
s2=shadedErrorBar(mdpq_range(1:t), CQ_allDOWN_mean_table(1,1:t), CQ_allDOWN_std_table(1,1:t), 'lineprops', '-red','patchSaturation',0.085);
s2.mainLine.Color = color{2};
s2.mainLine.LineWidth = 1.5;
xlabel('Time (belief update)')

ylabel('Belief');

% xline(2);
% text(2,0.4,'Initial belief update')
% xline(9);
% text(7,0.4,'Final belief update')
%legend('Policies causing up ' , 'Policies causing down', Location='northwest', fontsize=12);
%title('Posterior predictive distribution in policies leading down')
title(' Policy Prediction Optimization', fontsize=21);



sgtitle('Single-Layer Architectures with Disruption cases', fontsize=15)
fontsize(gcf,15,"points");



