

function [Q_allUP_mean, Q_allUP_std, Q_allDOWN_mean, Q_allDOWN_std, Q_allBASELINE_mean, Q_allBASELINE_std  ] = mdpq_posteriorconsequences(MDP)

% PARA MULTIPLES TRIALS 
% MDP.Q

% Posterior probability of each state conditioned on each policy 
% at the end of the trial after successive rounds of updating at each time point.

%1 cell per state factor. 
% Rows = states. 
% Columns = tau. 
% Third dimension = policy number.

numtrials = size(MDP,2);



for i = 1:numtrials
    for j = 1: size(MDP(1).P,1)
    pqi{i,j} = MDP(i).Q{1,1}(:,:,j); % MDP.Q
    end
end



Q1= [pqi{:,1}];
Q2= [pqi{:,2}];
Q3= [pqi{:,3}];


Q_allUP = [Q1(1,:); Q2(1,:); Q3(1,:)];
Q_allDOWN = [Q1(2,:); Q2(2,:); Q3(2,:)];
Q_allBASELINE = [Q1(3,:); Q2(3,:); Q3(3,:)];

Q_allUP_mean = mean(Q_allUP);
Q_allUP_std = std(Q_allUP);

Q_allDOWN_mean = mean(Q_allDOWN);
Q_allDOWN_std = std(Q_allDOWN);

Q_allBASELINE_mean = mean(Q_allBASELINE);
Q_allBASELINE_std = std(Q_allBASELINE);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure();
% % 
% % subplot(2,2,1)
% % plot(policy1', '-o');
% % xlabel('time');
% % legend('Up', 'Down', 'Left', 'Rigth'); 
% % title('Policy Up')
% % 
% % subplot(2,2,2)
% % plot(policy2', '-o')
% % xlabel('time');
% % legend('Up', 'Down', 'Left', 'Rigth'); 
% % title('Policy Down')
% % 
% % subplot(2,2,3)
% % plot(policy3', '-o')
% % xlabel('time');
% % legend('Up', 'Down', 'Left', 'Rigth'); 
% % title('Policy Left')
% % 
% % subplot(2,2,4)
% % plot(policy4', '-o')
% % xlabel('time');
% % legend('Up', 'Down', 'Left', 'Rigth'); 
% % title('Policy  Rigth')
% % 
% % 
% % sgtitle(sprintf(['  Posterior probability of each state \n conditioned on each policy ']));
% 
% %%%%%%%%%%%%%%%%%%%%
% policie_q = ([pqi{:,1}]');
% 
% policyup = policie_q(:,1)';
% policydown = policie_q(:,2)';
% policybaseline = policie_q(:,3:end)';




end