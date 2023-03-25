%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% GVG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%auxiliary script for hierarchical model

% MDP    - structure (see spm_MDP_VB)
% FACTOR - the hidden factors (at the second alevel) to plot
% T      - flag to return cell of expectations (at time T; usually 1)
% 
% [o1, s1, u1]  - observations, states, and actions 
% at level 1 (lower level) with the prescence of level 2 (higher level)
% 
%
% Learning lower level
% Learning higher level
%
%
% MDP(n).mdp(n) = 1st level
% MDP(n).o = 2nd level

%==========================================================================


%% Latex style
set(0,'defaulttextinterpreter','latex');



%% EXTRACT factors,  outcomes and dimentions
%--------------------------------------------------------------------------

% NUMBERS OF FACTORS 2nd LEVEL
% F1 = Sequence states
% F2 = Time
% F3 = Feedback
for i = 1:size(MDP(1).D,2)
  eval(sprintf('F%d = [1:i]', i)); 
end


% NUMBER OF TRIALS higherlevel
trials_l2 = numel(MDP);
nepochs_l2 = size(MDP(1).F,2);

for m = 1:numel(MDP)
    Ns = size(MDP(m).s(F1,:),2);
end




%% EXTRACT VARIABLES
%--------------------------------------------------------------------------

%% STATES


% FIRST LEVEL 

% TRUE STATES -  s
s_1h = [];
for i = 1:trials_l2
   s_1h = [ s_1h, MDP(i).mdp.s]; % true states states of the 1st layer hierarchical model
end
s1_1 = s_1h(F1,:);


% BELIEF STATES - X
% held about each of the hidden state factors (retrospectively) 

vbelief_s_L2 = [];

for i = 1: trials_l2
    for j = 1:2: nepochs_l2
    
        vbelief_s_L2 = [  vbelief_s_L2, MDP(i).mdp(j).X{F1}  ]; % belief states of the 1st layer with hierarchical model

    end
end


% SECOND LEVEL

% TRUE STATES -  s
s_2 = []; 
for i = 1:trials_l2
    s_2= [s_2, MDP(i).s];
end
s_2_F1 = s_2(F1,:);

% BELIEF STATES - X

vbelief_sF1 = [];

for i = 1:trials_l2
    vbelief_sF1= [vbelief_sF1, MDP(i).X{F1}];
end



%% OBSERVATIONS

% FIRST LEVEL - o

% MDP(1).mdp.o

% MDP(1).mdp.O - specified as the sampled outcome



% SECOND LEVEL
n_obsL2   = size(MDP(1).o, 2);
total_oL2 = (trials_l2 * n_obsL2);

vect_oL2 = [];

for i = 1:n_trials
        vect_oL2 = [ vect_oL2,  MDP(i).o]; % observations of the 2nd layer with hierarchical model
end

vect_oL2_thermo = vect_oL2(F1,:); % 2nd layer Factor 1 observations

%% VECTORS

% TRUE
s_1_L2 = (s1_1(1,:)==1);
s_2_L2 = (s1_1(1,:)==2);
s_3_L2= ( s1_1(1,:)==3);
s_4_L2 = (s1_1(1,:)==4);

if any(s1_1 >= 5)
    s_5_L = (s1_1(1,:) >= 5);
end

% BELIEF
s_belief_1_L2 = vbelief_s_L2(1,:);
s_belief_2_L2 = vbelief_s_L2(2,:);
s_belief_3_L2 = vbelief_s_L2(3,:);
s_belief_4_L2 = vbelief_s_L2(4,:);

if any(size(vbelief_s_L2(:,1),1) >= 5)
    s_belief_5_L2 = (vbelief_s_L2(5:end,:) >= 5);
end

% OBSERVATION
o_up_L2 = (vect_oL2_thermo(1,:)==1);
o_down_L2 = (vect_oL2_thermo(1,:)==2);
o_baseline_L2 = (vect_oL2_thermo(1,:)>=3);




%% ACTION SELECTION




%% THERMOMETER as 2nd LEVEL 
% Obtain the sequence of observations
%--------------------------------------------------------------------------
% if 1  =   +1
% if 2  =   -1
% if 3  =   +0

move = 0;
iter = 1;
vect_moveL2= zeros(1, total_oL2);
vect_thermoL2 = zeros(1, total_oL2);


% while loop observations 
while iter < total_oL2+1
   %fprintf('value of a: %d\n', vect_o( iter ) );

% IF rules over the timestep

   if vect_oL2_thermo(iter) == 1
       vect_moveL2(iter) = move + 1;

   elseif vect_oL2_thermo(iter) == 2
       vect_moveL2(iter) = move - 1;

   elseif vect_oL2_thermo(iter) == 3
      vect_moveL2(iter) =  move + 0; 
   end
      % The general rule output
   %fprintf('value of move: %f\n', vect_move(iter) );  

   S = sum(vect_moveL2(:));
  % fprintf('value of move: %f\n', vect_move(iter) );  
  % fprintf('value of S: %f\n', S ); 

     % RULE TO RESTRICT values between -5 and 5
   if S > 10
       S = 10;
   elseif S < -10
       S = -10;
   else
       S = S;
   end

    vect_thermoL2 (iter) = S;

   %tf = isequal(vect_move(iter),vect_o(iter)) % here we check if are equal
   % The general iterable
   iter = iter + 1;

end

%% LEARNING

% FIRST LEVEL VFE   B
Fb_1 = [];
for i = 1:numel(MDP)
   Fb_1 = [ Fb_1, MDP(i).mdp.Fb];
end
Fb_1 = nonzeros(Fb_1');
Fb1_1 = [];
for i = 1:2:size(Fb_1,1)
    Fb1_1 = [Fb1_1,  Fb_1(i) ];
end

% SECOND LEVEL

Fa_f1 = [];
Fa_f2 = [];
for i = 1:trials_l2

    Fa_f1 = [Fa_f1, MDP(i).Fa(1,1)];
    Fa_f2 = [Fa_f2, MDP(i).Fa(1,2)];

end




%% %%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Color PLOT
%--------------------------------------------------------------------------
%color = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0 0 0], [0 0 0 0.5], [0.4660 0.6740 0.1880], [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.954 0.345 0.543] }; % Cell array of colros.

color = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0 0 0], [0 0 0 0.5], [0 0.4470 0.7410 0.5], [0.8500 0.3250 0.0980 0.8], [0.4660 0.6740 0.1880], [0.954 0.345 0.543], [0.4 0.6 0.1 0.5], [0.954 0.345 0.543 0.5], [0.4660 0.6740 0.1880 0.5], [0.92 0.69 0.1125], [0.6350 0.0780 0.1840] }; % Cell array of colros.

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

%% PLOT THERMOMETER SIMULATION with 2nd LEVEL

figure();

p=plot(vect_thermoL2, '-ok');
p.MarkerFaceColor = 'k';
hold on
legend('Thermometer');
ylim([-10 10]);
title(['Thermometer with 2nd level']);
ylabel('Position in thermometer');
xlabel('Time step');
hold off



%% PLOT SEQUENCE OF STATES, BELIEFS, and OBSERVATIONS

figure()


subplot(3,1,1);

% TRUE STATES
p1=plot(s_1_L2,'.', 'Color', color{1},'MarkerSize',30);
p1.MarkerFaceColor = color{1};
hold on
plot(s_2_L2,'.', 'Color', color{2}, 'Markersize', 30);
plot(s_3_L2, '.', 'Color', color{7}, 'Markersize', 20);
plot(s_4_L2, '.', 'Color', color{12}, 'Markersize',20);
ylim([0.8 1.2]);
title(sprintf('True hidden states during time '));
yticks(1);
yticklabels('');
legend('sate 1', 'state 2', 'state 3', 'state 4')
xlabel('Time');
ylabel('States');

subplot(3,1,2); 
% INFERED STATES
p2=plot(s_belief_1_L2, '.', 'Color', color{1},'MarkerSize',30);
p2.MarkerFaceColor = color{1};
hold on
plot(s_belief_2_L2,'.', 'MarkerSize',30, "Color",color{2});
plot(s_belief_3_L2,'.', 'MarkerSize',20, "Color",color{7});
plot(s_belief_4_L2,'.', 'MarkerSize',20, "Color",color{12});
ylim([-0.2 1.2]);
title(sprintf('Belief about hidden states during time '));
%yticks(0 1);
%yticklabels('');
xlabel('Time');
ylabel('Posterior probability');


subplot(3,1,3); 
% OBSERVATIONS
p3= plot(o_up_L2,'.', 'MarkerSize',30);
p3.MarkerFaceColor= color{1};
hold on
plot(o_down_L2,'.', 'MarkerSize',30, "Color", color{2});
hold on
plot(o_baseline_L2,'.', 'MarkerSize',20, "Color", color{3});


ylim([0.8 1.2]);
title(sprintf('Observations during time '));
yticks(1);
yticklabels('');
xlabel('Time');
ylabel('Thermometer observation');

sgtitle(sprintf(['Hierarchical: \n Comparison between the true hidden states, \n infered hidden states \n and observations during time '])) ;




%% PLOT LEARNING

figure()
subplot(2,1,1)
plot(Fa_f1, LineWidth=2 );
hold on
plot(Fa_f2, LineWidth=2 );
legend('States', 'Feedback');
ylabel('Variational Free energy');
xlabel('Time');
title('2nd level Learning parameter of the likelihood ')

subplot(2,1,2)
plot(Fb1_1, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
ylabel('Variational Free energy');
xlabel('Time');
title('1sr level Learning parameter of the transition ')

sgtitle(sprintf([  '  Learning parameters  ']));

