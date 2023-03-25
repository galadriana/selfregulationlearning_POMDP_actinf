

%%
%
% GVG
% script to plot the actions and polcies
%
% 


% Latex style
set(0,'defaulttextinterpreter','latex');

%INTERACTIVO: DEFINE NOMBRE DEL EXPERIMENTO PARA LOS PLOTS
%--------------------------------------------------------------------------
if ~isa(prompt, 'char')
prompt = 'Define name of experiment:';
type_experiment= input(prompt, 's');
if isempty(type_experiment)
    txt = 'no_definido';
end
end 


%% EXTRACT factors and outcomes 
%--------------------------------------------------------------------------
if size(MDP(1).D,2) ==1

    if ~iscell(MDP(1).A)
    Ns = 1;
    No = 1;
    Nu = 1;
    X  = {MDP.X};
    C  = {MDP.C};
    else 
    Ns = size(MDP(1).A{:},2); %numel(MDP(1).B);                 
    No = size(MDP(1).A{:},1); % numel(MDP(1).A);                 
    X  = MDP.X;
    C  = MDP.C;  
    Nu = size(MDP(1).B{1},3) ;
    end
else
    Ns = numel(MDP(1).B);                 % number of hidden state factors
    No = numel(MDP(1).A);                 % number of outcome factors
    X  = [MDP.X];
    C  = MDP.C;
    for f = 1:Ns
        Nu(f) = size(MDP(1).B{f},3) > 1;
    end
end


% Nu = size(MDP(1).B{1},3) ;


%% Factors and outcomes to plot
%--------------------------------------------------------------------------

maxg  = 3;
gf  = 1:min(Ns,maxg); 
gg  = 1:min(No,maxg); 
nf  = numel(gf);
ng  = numel(gg);
Nf = numel(MDP(1).B);  

% Numero de trials 
%--------------------------------------------------------------------------
% Varia el tipo de plot para MDP.Q

numtrials = size(MDP,2);

%% EXTRACT VARIABLES
%--------------------------------------------------------------------------
if size(MDP(1).D,2) ==1
    if ~isa(MDP(1).P, 'double')
        fprintf('No Actions')
    else
        prob_action = [MDP.P].' ; % MDP.P(M,T) - probability of emitting action 1,...,M at time 1,...,T
        post_policy = [MDP.R].' ; % MDP.R      - conditional expectations over policies
        numeros_pol= size(MDP(1).R,1);
    end
else

    prob_action = [MDP.P];
    post_policy = [MDP.R];
    numeros_pol= size(MDP(1).R,1);

end

% MDP.P: Probability of the actions - posterior beliefs about control states
% Probability of emitting action
% The probability of emitting each particular action, 
% expressed as a softmax function of a vector 
% containing the probability of each action summed over each policy.

%% VFE
% Negative total variational free energy over time.
%--------------------------------------------------------------------------
if size(MDP(1).D,2) ==1
    F_total = [MDP.H];
else
    total_vfe = [MDP.H];
    F_total = [];
    for i = 1:2:numtrials
        F_total = [F_total, total_vfe(i)];
    end
end



%% Vector of actions
if size(MDP(1).D,2) ==1
    vect_actions= [MDP(:).u];
else
    actions = [MDP.u];
    vect_actions = actions(1,:);
end

%% EXTRACT FREE ENERGIES of the policies
%--------------------------------------------------------------------------

if size(MDP(1).D,2) ==1
    %%% VFE of policies (F)
    vfe_pol= [MDP.F].';
    VFE_rows= (vfe_pol)';
    %%% EFE of policies (G)
    efe_pol= [MDP.G].';
    efe_rows = (efe_pol)';
else
    %%% VFE of policies (F)
    %%% EFE of policies (G)
    VFE_rows= [];
     efe_rows = [];
    for i = 1:2:numtrials
        VFE_rows = [VFE_rows, MDP(i).F(:,1)];
        efe_rows = [efe_rows, MDP(i).G(:,1) ];
    end


end

%% Color set
color = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0 0 0], [0 0 0 0.5], [0 0.4470 0.7410 0.5], [0.8500 0.3250 0.0980 0.8], [0.4660 0.6740 0.1880], [0.954 0.345 0.543], [0.4 0.6 0.1 0.5], [0.954 0.345 0.543 0.5], [0.4660 0.6740 0.1880 0.5], [0.92 0.69 0.1125], [0.6350 0.0780 0.1840], [0.2081 0.1663 0.5292] }; % Cell array of colros.

% 1. blue
% 2. red
% 3. black
% 4. black-transparency
% 5. blue-transparency
% 6. red-transparency
% 7. green
% 8. pink
% 9. green-transparency 
% 10. pink-transparency
% 11. green-transparency
% 12. yellow
% 13. deep red 
% 14. deep blue

%% PLOT OF THE ACTIONS

act_up = (vect_actions(1,:)==1);
act_down = (vect_actions(1,:)==2);
act_left = (vect_actions(1,:)==3);
act_rigth = (vect_actions(1,:)==4);

if any(vect_actions>=5)
    act_lost = (vect_actions(1,:)>=5);
end

figure()

p1=plot(act_up, 'o','Color', color{1},'MarkerSize', 15, 'LineWidth',0.2);
p1.MarkerFaceColor = color{1};

hold on

p2=plot(act_down, 'o', 'Color', color{2}, 'MarkerSize', 15, 'LineWidth',0.2);
p2.MarkerFaceColor = color{2};


p2=plot(act_left, 'o', 'Color', color{7}, 'MarkerSize', 10, 'LineWidth',0.2);
p2.MarkerFaceColor = color{7};

p2=plot(act_rigth, 'o', 'Color', color{12}, 'MarkerSize', 10, 'LineWidth',0.2);
p2.MarkerFaceColor = color{12};

if exist('act_lost')
    plot(act_lost, '.', 'Color', color{14}, 'Markersize',20);
else
    ;
end

if size(MDP(1).X{gf(Nf)},1) > 12
    spm_spy(MDP(1).X{gf(Nf)}, 18, 1)
end

ylim([0.9 1.1]);
set(gca,'XTickLabel',{});
set(gca,'XTick', 1:size(MDP(1).X,2));
set(gca,'YTick', 1:numel(MDP(1).label.action{1}));
ylabel('Action in the first level');
legend('action up', 'action down', 'action left', 'action rigth');
xlabel('Time');
%%%%set(gca,'YTickLabel', MDP(1).label.action{1});
sgtitle(sprintf([type_experiment, ': Selected action during the simulation '])) ;


%% PLOTS probability of emitting action
%--------------------------------------------------------------------------
figure(); 

subplot(2,2,1);
p1=plot(prob_action(:,1), '-o','Color', color{1},'MarkerSize', 10, 'LineWidth',2);
p1.MarkerFaceColor = color{1};
ylabel('Probability');
xlabel('Time steps');
ylim([-0.1 1.1]);
legend('Action up');

subplot(2,2,2); 

p2=plot(prob_action(:,2), '-o','Color', color{2}, 'MarkerSize', 10, 'LineWidth',2);
p2.MarkerFaceColor = color{2};
ylabel('Probability');
xlabel('Time steps');
ylim([-0.1 1.1]);
legend('Action down');

subplot(2,2,3); 

p2=plot(prob_action(:,3), '-o','Color', color{7}, 'MarkerSize', 10, 'LineWidth',2);
p2.MarkerFaceColor = color{7};
ylabel('Probability');
xlabel('Time steps');
ylim([-0.1 1.1]);
legend('Action left');

subplot(2,2,4); 

p2=plot(prob_action(:,4), '-o','Color', color{12}, 'MarkerSize', 10, 'LineWidth',2);
p2.MarkerFaceColor = color{12};
ylabel('Probability');
xlabel('Time steps');
ylim([-0.1 1.1]);
legend('Action rigth');

sgtitle(sprintf([type_experiment,': Posterior belief of emiting an action'])) ;


%% PLOT OF THE VFE POLICIES
%--------------------------------------------------------------------------

figure();
plot(VFE_rows', 'LineWidth',2, 'Color',color{4});
ylabel('VFE');
xlabel('Time steps');
legend_names = "Policy" + string(1:numeros_pol);
legend(legend_names);
sgtitle(sprintf([type_experiment, ': \n Variational Free Energy of ', num2str(numeros_pol), ' policies']));

%% PLOT OF THE EFE POLICIES
%--------------------------------------------------------------------------

figure();
plot(efe_rows', 'LineWidth', 2, 'Color',color{4});
hold on
plot(efe_rows(1,:)', '.','Color', color{1}, MarkerSize=20);
plot(efe_rows(2,:)', '.','Color', color{2}, MarkerSize=20);
plot(efe_rows(3,:)', '.','Color', color{7}, MarkerSize=20);
plot(efe_rows(4,:)', '.','Color', color{12}, MarkerSize=20);
ylabel('EFE');
xlabel('Time steps');
legend_names = "Policy " + string(1:numeros_pol);
legend(legend_names);

sgtitle(sprintf([type_experiment, ': \n Expected Free Energy of ', num2str(numeros_pol), ' policies']));


%% Beliefs about each policy for each time step
%--------------------------------------------------------------------------
% Expectations over policies

pUN = [MDP.un];

pol1_un = pUN(1,:);
pol2_un = pUN(2,:);
pol3_un = pUN(3,:);
pol4_un = pUN(4,:);


figure(); 

subplot(2,2,1);

p1=plot(pol1_un, '-o','Color', color{1},'MarkerSize', 4, 'LineWidth',2);
p1.MarkerFaceColor = color{1};
ylabel('Probability');
xlabel('Updates');
ylim([-0.1 1.1]);
legend('policy up');

subplot(2,2,2); 

p2=plot(pol2_un , '-o','Color', color{2}, 'MarkerSize', 4, 'LineWidth',2);
p2.MarkerFaceColor = color{2};
ylabel('Probability');
xlabel('Updates');
ylim([-0.1 1.1]);
legend('policy down');


subplot(2,2,3); 

p2=plot(pol3_un , '-o','Color', color{7}, 'MarkerSize', 4, 'LineWidth',2);
p2.MarkerFaceColor = color{7};
ylabel('Probability');
xlabel('Updates');
ylim([-0.1 1.1]);
legend('policy left');

subplot(2,2,4); 

p2=plot(pol4_un , '-o','Color', color{12}, 'MarkerSize', 4, 'LineWidth',2);
p2.MarkerFaceColor = color{12};
ylabel('Probability');
xlabel('Updates');
ylim([-0.1 1.1]);
legend('policy rigth');


sgtitle(sprintf([type_experiment,': Beliefs about each policy for each time step'])) ;


%% PLOT posterior over policies
%--------------------------------------------------------------------------

% figure(3); 
% hold on
% for i = 1:size(MDP.R,1)
%     p{i} = post_policy(:,i)';
%     plot(p{i}, 'Color', color{i});
% end
% title('Posterior over policies');
% %subtitle('50 test');
% ylabel('Posterior');
% xlabel('Time steps');
% hold off
% legend_names = "Policy" + string(1:numeros_pol);
% legend(legend_names);


% SI FALLA

maxpol= find(post_policy == max(max(post_policy)));
minpol= find(post_policy == min(min(post_policy)));
[rowmin,colmin] = ind2sub(size(post_policy),minpol);
[rowmax,colmax] = ind2sub(size(post_policy),maxpol);

figure();
plot(post_policy, LineWidth=2);
hold on
%plot(rowmin,colmin, 'o', 'Color', [1 0.5 0], 'MarkerSize', 5, 'LineWidth', 5)
% % str1 = {'Minimal value policy:' num2str(post_policy(minpol))};
% % text(rowmin,colmin,str1, "FontSize", 10);
% % %plot(rowmax,colmax, 'o', 'Color', [0.6350 0.0780 0.1840], 'MarkerSize', 5, 'LineWidth', 5)
% % str2 = {'Maximal value policy:' num2str(post_policy(maxpol))};
% % text(0.250,10, str2, "FontSize", 15);
p= plot(post_policy(:,1), 'o-k');
p.MarkerFaceColor = 'k';
ylabel('Posterior');
xlabel('Time steps');
legend_names = "Policy" + string(1:numeros_pol);
legend(legend_names);
hold off
sgtitle(sprintf([type_experiment,': Posterior of the ', num2str(numeros_pol), ' policies']));


%% Posteriors over states under each policy at the end of the trial.
% PARA MULTIPLES TRIALS 

% Posterior probability of each state conditioned on each policy 
% at the end of the trial after successive rounds of updating at each time point.

%1 cell per state factor. 
% Rows = states. 
% Columns = tau. 
% Third dimension = policy number.

if numtrials > 2

for i = 1:numtrials
    for j = 1: size(MDP(1).P,1)
    pqi{i,j} = MDP(i).Q{1,1}(:,:,j); % MDP.Q
    end
end

policy1= [pqi{:,1}];
policy2= [pqi{:,2}];
policy3= [pqi{:,3}];
policy4= [pqi{:,4}];

figure();

subplot(2,2,1)
plot(policy1', '-o');
xlabel('time');
legend('Up', 'Down', 'Left', 'Rigth'); 
title('Policy Up')

subplot(2,2,2)
plot(policy2', '-o')
xlabel('time');
legend('Up', 'Down', 'Left', 'Rigth'); 
title('Policy Down')

subplot(2,2,3)
plot(policy3', '-o')
xlabel('time');
legend('Up', 'Down', 'Left', 'Rigth'); 
title('Policy Left')

subplot(2,2,4)
plot(policy4', '-o')
xlabel('time');
legend('Up', 'Down', 'Left', 'Rigth'); 
title('Policy  Rigth')


sgtitle(sprintf([type_experiment, ': \n Posterior probability of each state \n conditioned on each policy ']));


end







