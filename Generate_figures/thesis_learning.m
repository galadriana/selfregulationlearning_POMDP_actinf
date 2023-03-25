
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% GVG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxiliary script for learning dynamics of the model
% MDP - structure (see spm_MDP_VB)
% Fx - the learned parameters : Fa: Likelihood; Fb:transition, Fc: prefered
% outcom. 


%% PLOTS AND ANALYSIS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

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

%% SIMULATION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% EXTRACT VFE of learning parameters

%--------------------------------------------------------------------------
names = fieldnames(MDP(1)); 
match = {'Fa' 'Fb' 'Fc' 'Fd' };
TF= matches(names, match);

parameters= names(find(TF==1));

Fx = [];
for i = size(parameters, 2)
    a= parameters{i};
    Fx.(a) = MDP(1).(parameters{1});
end

%% %%%%%%%%% PLOTS of the  VFE LEARNING PARAMETERS
%--------------------------------------------------------------------------

%% A
if ~isfield(MDP(1), 'Fa')
       ;
else
    if size(MDP(1).D,2) ==1
        Fa = [MDP.Fa].'; % calculated as : MDP(m).Fb(f) = - spm_KL_dir(MDP(m).b{f},pB{m,f});
        figure();
        plot(Fa, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
        ylabel('Variational Free energy');
        xlabel('Time');
        sgtitle(sprintf([ type_experiment, '\n Learning parameter of the likelihood matrix ']));
    else
        Fa = [];
        for i = 1: size(MDP,2)
            Fa = [Fa, MDP(i).Fa(1,1) ];
        end
        figure();
        plot(Fa, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
        ylabel('Variational Free energy');
        xlabel('Time');
        sgtitle(sprintf([ type_experiment, '\n Learning parameter of the likelihood matrix ']));
    end
end




%% B
if ~isfield(MDP(1), 'Fb')
       ;
else
    if size(MDP(1).D,2) ==1
        Fb = [MDP.Fb].'; % calculated as : MDP(m).Fb(f) = - spm_KL_dir(MDP(m).b{f},pB{m,f});
        figure();
        plot(Fb, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
        ylabel('Variational Free energy');
        xlabel('Time');
    else
        Fb = [];
        for i = 1: size(MDP,2)
            Fb = [Fb, MDP(i).Fb(1,1) ];
        end
        figure();
        plot(Fb, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
        ylabel('Variational Free energy');
        xlabel('Time');
        sgtitle(sprintf([ type_experiment, '\n Learning parameter of the Transition matrix ']));

    end
end







%% C

if ~isfield(MDP(1), 'Fc')
       ;
else
    if size(MDP(1).D,2) ==1
        Fc = [MDP.Fc].'; % calculated as : MDP(m).Fb(f) = - spm_KL_dir(MDP(m).b{f},pB{m,f});
        figure();
        plot(Fc, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
        ylabel('Variational Free energy');
        xlabel('Time');
        sgtitle(sprintf([ type_experiment, '\n Learning parameter of the preference over the observations ']));
    else
        Fc = [];
        for i = 1: size(MDP,2)
            Fc = [Fc, MDP(i).Fc(1,1) ];
        end
        figure();
        plot(Fc, Color=[0.4940 0.1840 0.5560 0.9], LineWidth=2 );
        ylabel('Variational Free energy');
        xlabel('Time');
        sgtitle(sprintf([ type_experiment, '\n Learning parameter of the preference over the observations ']));

    end
end










