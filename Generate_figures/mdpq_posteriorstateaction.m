

function [policie_q, policyup, policydown, policybaseline] = mdpq_posteriorstateaction(MDP)

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

% policy1= [pqi{:,1}];
% policy2= [pqi{:,2}];
% policy3= [pqi{:,3}];
% policy4= [pqi{:,4}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure();
% 
% subplot(2,2,1)
% plot(policy1', '-o');
% xlabel('time');
% legend('Up', 'Down', 'Left', 'Rigth'); 
% title('Policy Up')
% 
% subplot(2,2,2)
% plot(policy2', '-o')
% xlabel('time');
% legend('Up', 'Down', 'Left', 'Rigth'); 
% title('Policy Down')
% 
% subplot(2,2,3)
% plot(policy3', '-o')
% xlabel('time');
% legend('Up', 'Down', 'Left', 'Rigth'); 
% title('Policy Left')
% 
% subplot(2,2,4)
% plot(policy4', '-o')
% xlabel('time');
% legend('Up', 'Down', 'Left', 'Rigth'); 
% title('Policy  Rigth')
% 
% 
% sgtitle(sprintf(['  Posterior probability of each state \n conditioned on each policy ']));

%%%%%%%%%%%%%%%%%%%%
policie_q = ([pqi{:,1}]');

policyup = policie_q(:,1)';
policydown = policie_q(:,2)';
policybaseline = policie_q(:,3:end)';




end