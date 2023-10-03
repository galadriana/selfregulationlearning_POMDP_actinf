%% CASE 1
%% GENERATE PLOT TO VISUALIZE THE 3 subcases of behavior: Early regulator, early downregulator, disruption


clear all
close all


% Path

cd 'C:\Users\gadri\spm12\toolbox\DEM\GVG_and_tutorial_pomdp\Supplementary code\GVG_DEM\hierarchical\Thesis\Model_1\reproducibility_test\case1b_1test'
% LOAD



% Performance
load('test_reproducibility_100rnf.mat')

%INDEX



load('index_positive_100_prepro.mat')
load('index_down_100_prepro.mat')
load('index_disruption.mat')

% STATS
earlyup_mean = mean(thermometer_table(index_positive_100,:));
earlyup_std = std(thermometer_table(index_positive_100,:));

earlydown_mean = mean(thermometer_table(index_down_100,:));
earlydown_std = std(thermometer_table(index_down_100,:));

disrup_mean = mean(thermometer_table(index_disruption,:));
disrup_std = std(thermometer_table(index_disruption,:));

%% % Color set
color = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0 0 0], [0 0 0 0.5], [0 0.4470 0.7410 0.5], [0.8500 0.3250 0.0980 0.8], [0.4660 0.6740 0.1880], [0.954 0.345 0.543], [0.4 0.6 0.1 0.5], [0.954 0.345 0.543 0.5], [0.4660 0.6740 0.1880 0.5], [0.92 0.69 0.1125], [0.6350 0.0780 0.1840], [0.2081 0.1663 0.5292] }; % Cell array of colros.
mdpq_range = int16(linspace(1,600,600));

figure()
subplot(3,1,1)
s1up = shadedErrorBar(mdpq_range, earlyup_mean, earlyup_std, 'lineprops', '-blue','patchSaturation',0.085);
s1up.mainLine.Color = color{1};
s1up.mainLine.LineWidth = 1.5;
ylim([-15, 15])
ylabel('Thermometer ', fontsize=13);
title('Early up-regulators', fontsize=13);

subplot(3,1,2)
s1down = shadedErrorBar(mdpq_range, earlydown_mean, earlydown_std, 'lineprops', '-red','patchSaturation',0.085);
s1down.mainLine.Color = color{2};
s1down.mainLine.LineWidth = 1.5;
ylim([-15, 15])
ylabel('Thermometer ', fontsize=13);
title('Early down-regulators', fontsize=13)

subplot(3,1,3)
s1disrup = shadedErrorBar(mdpq_range, disrup_mean, disrup_std, 'lineprops', '-black','patchSaturation',0.085);
s1disrup.mainLine.Color = color{3};
s1disrup.mainLine.LineWidth = 1.5;
ylim([-15, 15])
ylabel('Thermometer ', fontsize=13);
xlabel('Time steps', fontsize=13);
title('Disruption dynamic', fontsize=13);


sgtitle('Behavior profiles of the single layer model', fontsize=22)
