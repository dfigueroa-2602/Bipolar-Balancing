Variables;
% Discrete Plant
[Ad,Bd,~,~] = ssdata(c2d(ss_Parallel_gamma,Ts,'zoh'));

% Delayed Plant
Ad_delay = [Ad Bd; zeros(nu,nx) zeros(nu)]; Bd_delay = [zeros(nx,nu);eye(nu)];
C_delay = [C_gamma zeros(nx,nu)];

% Discrete Integrative Plant (The one that the gains will be designed for)
Ha = [1 0 0];
Ad_aug = [ Ad_delay             zeros(nx+nu,nu);
          [-Ha*Ts*C_gamma zeros(nu)]     eye(nu)    ];
Bd_aug = [ Bd_delay; zeros(nu)];

Tsim = 0.3; T_settling = 40e-3; A_ref = 5;

search = 1;

beta_c = 3e-5; n_part = 100; iter = 50;
c1 = 2.05; c2 = 2.05; R_mode = 0;
Qmax = 6; Rmax = 2; Qmax_vel = 1; Rmax_vel = 0.1; rang_coef = 0.6;
nQ = nx + 2*nu; nR = size(Bd_aug,2); dim = nQ + (R_mode==1)*(nu/2);

[vel_clamp,Kap,swarm,space_range] = Swarm_Init(n_part,R_mode,0,nQ,nR,c1,c2,Qmax,Rmax,0,rang_coef,Qmax_vel,Rmax_vel,0);
xmax_val = 6; xmin = repmat(-xmax_val,dim,1); xmax = repmat(xmax_val,dim,1);

rng(1,'twister');
if search == 1
    swarms = {};
    % Randomize the seed, then open a pool of workers
    rng('shuffle'); if isempty(gcp('nocreate')), parpool; end

    % Initialization of important matrices
    b_fitness = zeros(iter,1); fitness = inf(1,n_part);
    gbest_vec = zeros(iter,1);
    % Initialize particles
    for i = 1:iter
        parfor n = 1:n_part
            try
                sw = swarm(n,:,:);
                sw = squeeze(sw(1,1,:));
                fitness(n) = LQR_Search(0,Ts,0,[],R_mode,sw,Ad_aug,Bd_aug, ...
                    Ad,Bd,[],[],Ha,beta_c,Vdc,w,A_ref,Tsim,T_settling);
            catch
                fitness(n) = 1e6;
                disp(['Evaluation for particle no. ' num2str(n) ' was aborted']);
            end
        end
        for n = 1:n_part
            if fitness(n) < swarm(n,4,1)
                swarm(n,3,:) = swarm(n,1,:);
                swarm(n,4,1) = fitness(n);
            end
        end
        PSO_Algorithm;
    end

    [Kx, Kr, Ku, Ki] = Final_Value(0,b_swarm,w,[],Ad_aug,Bd_aug,nx,nu,0,0,0);
    delete(gcp('nocreate'));
end

% currentFolder = pwd;
% C_Folder = extractBefore(currentFolder,'MATLAB');
% C_folder = append(C_Folder,'PLECS/');
% cd(C_folder)
% save HDT_Gains.mat Kx Ku Kr
% cd(currentFolder);
% 
% Kgain2Ccode({Kx, Ku, Kr, Ard, Brd, Ad_aug, Bd_aug},{'Kx_val', 'Ku_val', 'Kr_val', 'Ard_val', 'Brd_val', 'Ad_aug_val', 'Bd_aug_val'},'Matrices')