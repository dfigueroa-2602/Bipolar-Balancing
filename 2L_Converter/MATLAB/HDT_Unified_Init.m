Variables_HDT;

ss_HDT = ss(A,B,eye(nx),[]); [Ad,Bd,~,~] = ssdata(c2d(ss_HDT,Ts,'zoh'));

Ad_delay = [Ad Bd; zeros(nu,nx) zeros(nu)]; Bd_delay = [zeros(nx,nu);eye(nu)];

Ar1 = [0 w; -w 0]; Br1 = [0;1]; Ar3 = [0 (3*w); -(3*w) 0]; Ar5 = [0 (5*w); -(5*w) 0]; Ar7 = [0 (7*w); -(7*w) 0];
[Ard1,Brd1,~,~] = ssdata(c2d(ss(Ar1,Br1,eye(2),[]),Ts,'tustin'));
[Ard3,Brd3,~,~] = ssdata(c2d(ss(Ar3,Br1,eye(2),[]),Ts,'tustin'));
[Ard5,Brd5,~,~] = ssdata(c2d(ss(Ar5,Br1,eye(2),[]),Ts,'tustin'));
[Ard7,Brd7,~,~] = ssdata(c2d(ss(Ar7,Br1,eye(2),[]),Ts,'tustin'));
Ard1 = blkdiag(Ard1,Ard1,Ard1,Ard1); Brd1 = blkdiag(Brd1,Brd1,Brd1,Brd1);
Ard3 = blkdiag(Ard3,Ard3,Ard3,Ard3); Brd3 = blkdiag(Brd3,Brd3,Brd3,Brd3);
Ard5 = blkdiag(Ard5,Ard5,Ard5,Ard5); Brd5 = blkdiag(Brd5,Brd5,Brd5,Brd5);
Ard7 = blkdiag(Ard7,Ard7,Ard7,Ard7); Brd7 = blkdiag(Brd7,Brd7,Brd7,Brd7);

w_vec = [1 3 5 7]*w; n_h = length(w_vec);

Ard_HDT = blkdiag(Ard1,Ard3,Ard5,Ard7); Brd_HDT = [Brd1;Brd3;Brd5;Brd7];
nxr = size(Brd_HDT,1); nur = size(Brd_HDT,2);

C_delay = [C zeros(nx,nu)];

currentFolder = pwd;
C_Folder = extractBefore(currentFolder,'MATLAB');
C_folder = append(C_Folder,'PLECS/');
cd(C_folder)
save HDT_Resonant.mat Ard_HDT Brd_HDT
cd(currentFolder)

Ha = [0 0 1 0 0 0 0 0 0 0 ;
      0 0 0 1 0 0 0 0 0 0 ;
      0 0 0 0 1 0 0 0 0 0 ;
      0 0 0 0 0 1 0 0 0 0];

Hx = C_delay(3:6,:);

Ad_aug = [Ad_delay  zeros((nx + nu),nxr) ;
          -Brd_HDT*Hx   Ard_HDT                 ];
Bd_aug = [Bd_delay; zeros(nxr,nu)];

Tsim = 0.3; V_ref = 100; i_ref = 1; T_settling = 40e-3;
A_ref = [V_ref; i_ref];

search = 1;

beta_c = 3e-5; n_part = 100; iter = 50;
c1 = 2.05; c2 = 2.05; R_mode = 0;
Qmax = 4; Rmax = 2; Qmax_vel = 1; Rmax_vel = 0.1; rang_coef = 0.6;
nQ = nx/2 + nu/2 + n_h; nR = size(Bd_aug,2)/2; dim = nQ + (R_mode==1)*(nu/2);

[vel_clamp,Kap,swarm,space_range] = Swarm_Init(n_part,R_mode,0,nQ,nR,c1,c2,Qmax,Rmax,0,rang_coef,Qmax_vel,Rmax_vel,0);
xmax_val = 4; xmin = repmat(-xmax_val,dim,1); xmax = repmat(xmax_val,dim,1);

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
                fitness(n) = LQR_Search(1,Ts,n_h,w_vec,R_mode,sw,Ad_aug,Bd_aug, ...
                    Ad,Bd,Ard_HDT,Brd_HDT,Ha,beta_c,Vdc,w,A_ref,Tsim,T_settling);
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
    % [Kx, Kr, Ku, Ki] = Final_Value(Res_mode,b_swarm,w,w_vec,A_aug,B_aug,nx,nu,nxr,n_h,R_mode)
    [Kx_HDT, Kr_HDT, Ku_HDT, Ki] = Final_Value(1,b_swarm,w,w_vec,Ad_aug,Bd_aug,nx,nu,nxr,n_h,0);
    delete(gcp('nocreate'));
end

currentFolder = pwd;
C_Folder = extractBefore(currentFolder,'MATLAB');
C_folder = append(C_Folder,'PLECS/');
cd(C_folder)
save HDT_Gains.mat Kx_HDT Ku_HDT Kr_HDT
cd(currentFolder);

% Kgain2Ccode({Kx, Ku, Kr, Ard, Brd, Ad_aug, Bd_aug},{'Kx_val', 'Ku_val', 'Kr_val', 'Ard_val', 'Brd_val', 'Ad_aug_val', 'Bd_aug_val'},'Matrices')