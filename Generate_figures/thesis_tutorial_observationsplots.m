
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% GVG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTS AND ANALYSIS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

set(0,'defaulttextinterpreter','latex');

%INTERACTIVO: DEFINE NOMBRE DEL EXPERIMENTO PARA LOS PLOTS
%--------------------------------------------------------------------------

prompt = 'Define name of experiment:';
type_experiment= input(prompt, 's');
if isempty(type_experiment)
    txt = 'no_definido';
end

%% POMDP

%--------------------------------------------------------------------------
% if ~iscell(MDP(1).A)
%     Ns = 1;
%     No = 1;
%     Nu = 1;
%     X  = {MDP.X};
%     C  = {MDP.C};
% else s
%     Ns = size(MDP(1).A{:},2); %numel(MDP(1).B);                 
%     No = size(MDP(1).A{:},1); % numel(MDP(1).A);                 
%     X  = MDP.X;
%     C  = MDP.C;  
% end

nf = numel(MDP(1).D);

if nf == 1

    Ns = numel(MDP(1).B);                 % number of hidden state factors
    No = numel(MDP(1).A);                 % number of outcome factors
    X  = [MDP.X];
    C  = MDP.C;
    for f = 1:Ns
        Nu(f) = size(MDP(1).B{f},3) > 1;
    end
else
    ;
end

if nf > 1
    Ns = numel(MDP(1).B);                 % number of hidden state factors
    No = numel(MDP(1).A);                 % number of outcome factors
    X  = [MDP.X];
    C  = MDP.C;
    for f = 1:Ns
        Nu(f) = size(MDP(1).B{f},3) > 1;
    end
else

    Ns = 1;
    No = 1;
    Nu = 1;
    X  = {MDP.X};
    C  = {MDP.C};
end

% Negative total variational free energy over time.
%--------------------------------------------------------------------------
F_total = [MDP.H];

% factors and outcomes to plot
%--------------------------------------------------------------------------
maxg  = 3;
gf  = 1:min(Ns,maxg); 
gg  = 1:min(No,maxg); 
nf  = numel(gf);
ng  = numel(gg);




%% %%%%%%%%% PLOTS of the MDP 
%--------------------------------------------------------------------------
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

%% STATES DYNAMICS of the SIMULATION
%--------------------------------------------------------------------------
% Overall posteriors over states at the end of the trial (MDP.X).
% 1 modality observation: visual - location {up, down, baseline}


% Obtain the serie of states of the 'brain ' 
%--------------------------------------------------------------------------
if nf == 1

    n_trials    = size(MDP,2);
    n_states   = size(MDP(1).s,2);
    total_s = (n_trials * n_states);
    vect_s = zeros(n_trials, n_states);
    
    for i = 1:n_trials
        for j = 1:n_states
            vect_s(i,j) = [MDP(i).s(1,j)];
        end
    end

else

    n_trials = size(MDP,2);
    n_states = size(MDP(1).s,2);
    total_s = (n_trials * n_states);
    vect_s = [];
    for i = 1:n_trials
        vect_s = [vect_s, MDP(i).s(1,:)];
    end

end


vect_s = reshape(vect_s', 1, total_s );

%% Posterior beliefs about hidden states - Beliefs held about each of the hidden state factors (retrospectively) 
% Or CONDIITONED by the OBSERVATIONS
%--------------------------------------------------------------------------
% (i.e., after all observations have been made). 
% Black shading indicates a probability of one; white of zero.

if nf == 1
    vbelief_s = cell2mat(X);
else 
   vbelief_s = [];
    for i = 1:2:size(F_total,2)
     vbelief_s = [  vbelief_s, cell2mat(X(i)) ];
    end
end



level = 500;
n = ceil(level/2);
cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
cmap = [cmap1; cmap2(2:end, :)];

if size(MDP(1).X,2) <= 1
    minval= min(min(abs(MDP(1).X{:,:}))) ;
    maxval= max(max(abs(MDP(1).X{:,:})));
    largoX= linspace(minval, maxval, 10);
    lab1 = num2str(largoX(1));
else
    for n = 1: size(MDP(1).X,2)
        minval= min(min(abs(MDP(1).X{:,n}))) ;
        maxval= max(max(abs(MDP(1).X{:,n})));
        largoX= linspace(minval, maxval, 10);
        lab1 = num2str(largoX(1));
    end
end

% Obtain the serie of conditional expectations over hidden states 
%--------------------------------------------------------------------------

dob_s = zeros(Ns, n_states, n_trials );

if nf == 1
    for i = 1:n_trials
        for j = 1: Ns
            for k = 1: n_states
                dob_s(j, k, i) = [MDP(i).X{1,1}(j,k)];
            end
        end
    end
    dob_s = reshape(dob_s, [Ns, n_trials * n_states] );
else
    dob_s = vbelief_s;
end




%%  DYNAMICS of the SIMULATION - OBSERVATIONS
%--------------------------------------------------------------------------
% 1 modality observation

% Obtain the serie of observation 
%--------------------------------------------------------------------------
n_trials    = size(MDP,2);
n_obs   = size(MDP(1).o,2);
total_o = (n_trials*n_obs);
vect_o = zeros(n_trials, n_obs);

for i = 1:n_trials
    for j = 1:n_obs
        vect_o(i,j) = [MDP(i).o(1,j)];
    end
end

vect_o= reshape(vect_o', 1, total_o );

% Obtain the sequence of observations
%--------------------------------------------------------------------------
% if 1  =   +1
% if 2  =   -1
% if 3  =   +0

move = 0;
iter = 1;
vect_move= zeros(1, total_o);
vect_thermo = zeros(1, total_o);


% while loop observations 
while iter < total_o+1
   %fprintf('value of a: %d\n', vect_o( iter ) );

% IF rules over the timestep

   if vect_o(iter) == 1
       vect_move(iter) = move + 1;

   elseif vect_o(iter) == 2
       vect_move(iter) = move - 1;

   elseif vect_o(iter) == 3
      vect_move(iter) =  move + 0; 
   end
      % The general rule output
   %fprintf('value of move: %f\n', vect_move(iter) );  

   S = sum(vect_move(:));
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

   vect_thermo(iter) = S;

   %tf = isequal(vect_move(iter),vect_o(iter)) % here we check if are equal
   % The general iterable
   iter = iter + 1;

end


%% VECTOR OF ERRORS
% Errors between all observations

vs = vect_s;
vs(vs==4)=3;
S = vs ~= vect_o;
meanabsoluteerror =  mean(abs(vs - vect_o));
meansqerr = mean((vs-vect_o).^2);



%% Vectors True hidden states, infered hidden states and observations 
%                       divided by COLOR

s_1 = (vect_s(1,:)==1);
s_2 = (vect_s(1,:)==2);
s_3 = (vect_s(1,:)==3);
s_4 = (vect_s(1,:)==4);

if any(vect_s >= 5)
    s_5 = (vect_s(1,:) >= 5);
end

s_belief_1 = vbelief_s(1,:);
s_belief_2 = vbelief_s(2,:);
s_belief_3 = vbelief_s(3,:);
s_belief_4 = vbelief_s(4,:);

if any( size(vbelief_s(:,1),1)  >= 5)
    s_belief_5 = (vbelief_s(5:end,:) >= 5);
end


o_up = (vect_o(1,:)==1);
o_down = (vect_o(1,:)==2);
o_baseline = (vect_o(1,:)==3);

% if any(vect_o==3)
%     o_green = (vect_o(1,:)==3);
% end


%% PLOT THERMOMETER SIMULATION 
%--------------------------------------------------------------------------
figure();
p=plot(vect_thermo, '-ok');
p.MarkerFaceColor = 'k';
hold on
legend('Thermometer');
ylim([-10 10]);
title(['Thermometer observation  - ', type_experiment]);
ylabel('Position in thermometer');
xlabel('Time step');
hold off



%% PLOT SEQUENCE OF STATE INFERENCE 
figure()

subplot(3,1,1);

% TRUE STATES
p1=plot(s_1,'.', 'Color', color{1},'MarkerSize',30);
p1.MarkerFaceColor = color{1};
hold on
plot(s_2,'.', 'Color', color{2}, 'Markersize', 30);
plot(s_3, '.', 'Color', color{7}, 'Markersize', 20);
plot(s_4, '.', 'Color', color{12}, 'Markersize',20);

if exist('s_5')
    plot(s_5, '.', 'Color', color{14}, 'Markersize',20);
else
    ;
end

ylim([0.8 1.2]);
title(sprintf('True hidden states during time '));
yticks(1);
yticklabels('');
legend('sate 1', 'state 2', 'state 3', 'state 4')
xlabel('Time');
ylabel('States');

subplot(3,1,2); 
% INFERED STATES
p2=plot(s_belief_1, '.', 'Color', color{1},'MarkerSize',30);
p2.MarkerFaceColor = color{1};
hold on
plot(s_belief_2,'.', 'MarkerSize',30, "Color",color{2});
plot(s_belief_3,'.', 'MarkerSize',20, "Color",color{7});
plot(s_belief_4,'.', 'MarkerSize',20, "Color",color{12});

if exist('s_belief_5')
    plot(s_belief_5, '.', 'Color', color{14}, 'Markersize',20);
else
    ;
end

ylim([-0.2 1.2]);
title(sprintf('Belief about hidden states during time '));
%yticks(0 1);
%yticklabels('');
xlabel('Time');
ylabel('Posterior probability');


subplot(3,1,3); 
% OBSERVATIONS
p3= plot(o_up,'.', 'MarkerSize',30);
p3.MarkerFaceColor= color{1};
hold on
plot(o_down,'.', 'MarkerSize',30, "Color", color{2});
hold on
plot(o_baseline,'.', 'MarkerSize',20, "Color", color{3});
% if any(vect_o==3)
% 
% plot(o_green, '.', 'MarkerSize',30, "MarkerFaceColor", color{7});
% 
% end

ylim([0.8 1.2]);
title(sprintf('Observations during time '));
yticks(1);
yticklabels('');
xlabel('Time');
ylabel('Thermometer observation');

sgtitle(sprintf([type_experiment, ': \n Comparison between the true hidden states, \n infered hidden states \n and observations during time '])) ;



%% VECTOR OF ERRORS

% % Errors between all observations
% vs = vect_s;
% vs(vs==4)=3;
% S = vs ~= vect_o;
% meanabsoluteerror =  mean(abs(vs - vect_o));
% meansqerr = mean((vs-vect_o).^2);
% 
% figure()
% plot(S, '*', 'Color', color{13}, 'MarkerSize',30);
% str = {'Mean square error of:' num2str(meansqerr)};
% text(size(S,2)/2, 1.2, str, "FontSize", 15);
% ylim([0.5 1.5]);
% %yticks([]);
% ylabel('Error', "FontSize", 15);
% xlabel('Time', "FontSize", 15);
% title ([type_experiment,': Errors between state and observation'], "FontSize", 15);

%% CONDITIONAL EXPECTATION


%%%%%%%%%%%% MEJOR ACA PONER EL ERROR 
%%%  veces en que era esperado el estado y no fue acertado

ce_1 = (dob_s(1,:) > 0.9);
ce_2 = (dob_s(2,:) > 0.9);
ce_3 = (dob_s(3,:) > 0.9);
ce_4 = (dob_s(4,:) > 0.9);

%%% DIFFERENCE betweeb expected state and true state 
error_s1 = (((sum((s_1==ce_1)==0))*100)/total_s);
error_s2 = (((sum((s_2==ce_2)==0))*100)/total_s);
error_s3 = (((sum((s_3==ce_3)==0))*100)/total_s);
error_s4 = (((sum((s_4==ce_4)==0))*100)/total_s);


figure();
subplot(2,2,1);
%State 1 - UP inference 
plot(s_1,'.', 'MarkerSize',30, "MarkerFaceColor", color{1});
hold on
plot(ce_1, '.k', 'Marker','+', MarkerSize=10);
xlabel('Time');
ylabel('1 inference');
ylim([-0.1 1.1]);
str = {'Percentage error of:' num2str(error_s1) '%'};
text(size(error_s1,2)/2, 1, str, "FontSize", 15);
set(gca,'YTick',0:numel(MDP(1).label.name{gf(1)}));
set(gca,'YTickLabel',{'Abscent', 'Present'});

subplot(2,2,2); 
plot(s_2,'.', 'MarkerSize',30, "MarkerFaceColor", color{1});
hold on
plot(ce_2, '.k', 'Marker','+', MarkerSize=10);
str = {'Percentage error of:' num2str(error_s2) '%'};
text(size(error_s2,2)/2, 1, str, "FontSize", 15);
xlabel('Time');
ylim([-0.1 1.1]);
ylabel('2 inference');
set(gca,'YTick',0:numel(MDP(1).label.name{gf(1)}));
set(gca,'YTickLabel',{'Abscent', 'Present'});

subplot(2,2,3); 
plot(s_3,'.', 'MarkerSize',30, "MarkerFaceColor", color{1});
hold on
plot(ce_3, '.k', 'Marker','+', MarkerSize=10);
str = {'Percentage error of:' num2str(error_s3) '%'};
text(size(error_s3,2)/2, 1, str, "FontSize", 15);
xlabel('Time');
ylim([-0.1 1.1]);
ylabel('3 inference');
set(gca,'YTick',0:numel(MDP(1).label.name{gf(1)}));
set(gca,'YTickLabel',{'Abscent', 'Present'});

subplot(2,2,4); 
plot(s_4,'.', 'MarkerSize',30, "MarkerFaceColor", color{1});
hold on
plot(ce_4, '.k', 'Marker','+', MarkerSize=10);
str = {'Percentage error of:' num2str(error_s4) '%'};
text(size(error_s4,2)/2, 1, str, "FontSize", 15);
xlabel('Time');
ylim([-0.1 1.1]);
ylabel('4 inference');
set(gca,'YTick',0:numel(MDP(1).label.name{gf(1)}));
set(gca,'YTickLabel',{'Abscent', 'Present'});

sgtitle(sprintf([type_experiment,': Belief states and errors \n confidence with a 90/% confidence '])) ;


%% VFE
figure();
plot(F_total, LineWidth=2, Color=color{8} );
ylabel('Variational Free energy');
xlabel('Time');
subtitle ('averaged across states and policies at each time point');
title ([type_experiment, ': Total energy over time'], "FontSize", 15);


%%%


