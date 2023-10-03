

function [thermovector, thermomovement ] = thermometer(MDP)

%% Rename
%MPD = MDP;

%% Obtain the serie of observation 
%--------------------------------------------------------------------------

% TYPE OF MDP

if isfield(MDP, 'MDP')
% Is a hierarchical structure
    % MDP.mdp.o
    los_f = size(MDP(1).D,2);
    for i = 1:los_f
      eval(sprintf('F%d = [1:i]', i)); 
    end

    trials_l2 = numel(MDP);


    if numel(MDP) ==1
        vect_oL1 = [MDP.mdp.o];
        vect_oL1_thermo = vect_oL1(1,:);
        vect_oL1_NOTthermo = vect_oL1(2,:);
        
        n_mdpL1 = size([MDP(1).mdp.o],2);
        total_oL1 =  trials_l2 * n_mdpL1 ;
    else
        n_obsL1   = numel(MDP(1).o(1,:));
        total_oL1 = (trials_l2 * n_obsL1);
        vect_oL1 = [];
        for i=1:trials_l2
            vect_oL1 = [vect_oL1, MDP(i).mdp.o];
        end
        vect_oL1_thermo = vect_oL1(F1,:);
    end
    
    
        % SECOND LEVEL
        % MDP(1).mdp.o
        n_obsL2   = size(MDP(1).o, 2);
        total_oL2 = (trials_l2 * n_obsL2);
        
        vect_oL2 = [];
        
        for i = 1:trials_l2
                vect_oL2 = [ vect_oL2,  MDP(i).o]; % observations of the 2nd layer with hierarchical model
        end
        
        vect_oL2_thermo = vect_oL2(F1,:); % 2nd layer Factor 1 observations
        vect_oL2_NOTthermo = vect_oL2(size(F2,2),:); % 2nd layer Factor 2 observations
    
        vect_o = vect_oL2_thermo;
        total_o = total_oL2;



else

% Is a 1 layer structure
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

end



%% THERMOMETER 
% Obtain the sequence of observations
%--------------------------------------------------------------------------
% if 1  =   +1
% if 2  =   -1
% if 3  =   +0

move = 0;
iter = 1;
vect_move= zeros(1, total_o);
vect_thermo = zeros(1, total_o);



% Init loop 

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




%% PLOT THERMOMETER SIMULATION 
%--------------------------------------------------------------------------
figure();
p=plot(vect_thermo, '-');
p.MarkerFaceColor = 'k';
hold on
legend('Thermometer');
ylim([-11 11]);
yline(0,'-', 'Baseline','LineWidth',2);
title(['Thermometer observation ']);
ylabel('Position in thermometer');
xlabel('Time step');
hold off

thermomovement = vect_move;
thermovector = vect_thermo;

end
