%% COMPARISON OF TIME TO REACH UP REGULATION as BOX PLOT


clear all
close all

%% BINARIZE AND EXTRACT TIME TO REACH UP


% CASE 1
case1_thermo = load('test_reproducibility_100rnf.mat');
case1_thermo.thermometer_table;

index_time_zero_case1 = [];



for i = 1:size(case1_thermo.thermometer_table,1);

    arr_zero = min(find(any(case1_thermo.thermometer_table(i,:),1 >0)));
            index_time_zero_case1 = [index_time_zero_case1; arr_zero];

end


% CASE 2

load('test_reproducibility_case2.mat')

case2_thermo =  load('test_reproducibility_case2.mat');
case2_thermo.thermometer_table;

index_time_zero_case2 = [];


for i = 1:size(case2_thermo.thermometer_table,1)

    arr_zero = min(find(any(case2_thermo.thermometer_table(i,:),1 >0)));
            index_time_zero_case2 = [index_time_zero_case2; arr_zero];

end

% CASE 3


load('test_reproducibility_case3_3.mat')

case3_thermo =  load('test_reproducibility_case3_3.mat');
case3_thermo.thermometer_table_3;

index_time_zero_case3 = [];


for i = 1:size(case3_thermo.thermometer_table_3,1)

    arr_zero = min(find(any(case3_thermo.thermometer_table_3(i,:),1 >0)));
            index_time_zero_case3 = [index_time_zero_case3; arr_zero];

end

%% VISUALIZATION Three CASES
 

three_cases = [index_time_zero_case1 index_time_zero_case2 index_time_zero_case3 ];
h=boxplot( three_cases);
set(h,'linewidth',2) 
%legend('Time to reach up regulation')

set(gca,'xtick',1:3, 'xticklabel',{'Single Layer', 'Single layer with Policy priors', 'Hierarchical'},'Fontsize',15)

color = [[0 0.4470 0.7410],  [0.4 0.6 0.1],[0.8500 0 0.0980], 'b', 'b', 'b']; 
h = findobj(gca,'Tag','Box'); 
for j=1:length(h) 
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5); 
end 

ylabel('Time steps');
ylim([-5 101])
xlabel('Generative Model')
title('Time to reach up regulation')
%set(gca,'xtick',1:2, 'xticklabel',{'Single Layer','Hierarchical'})

