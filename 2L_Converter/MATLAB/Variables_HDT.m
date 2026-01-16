clc
fsw = 20e3; Ts = 1/(2*fsw); Vg = 10e3; Vs = 400; Vdc = 700;
NDT = Vg/Vs; NCT = 5;
Lfs = 200e-6; Rfs = 100e-3; Cfs = 12e-6;
Cdc = 5000e-6;
Lfp = 200e-6; Ly = 100e-6; Ry = 5e-3; Rfp = 100e-3; Cfp = 20e-6;

w = 2*pi*50;

T = [2/3 -1/3 -1/3; 0 1/sqrt(3) -1/sqrt(3); 1/3 1/3 1/3];
KDT = [1 0 -1; -1 1 0; 0 -1 1];

%%%%%%%%%%%%%%%%%%%%% SERIES CONVERTER %%%%%%%%%%%%%%%%%%%%%%

As = [-Rfs/Lfs   0         -1/Lfs      0     ;
       0        -Rfs/Lfs    0         -1/Lfs ;
       1/Cfs     0          0          0     ;
       0         1/Cfs      0          0    ];

Bs = [ 1/Lfs     0     ;
       0         1/Lfs ;
       0         0     ;
       0         0    ];

Ps = [ 0         0     ;
       0         0     ;
      -1/Cfs     0     ;
       0        -1/Cfs];

Cs = eye(size(As, 1));
Ds = zeros(size(Cs, 1), size(Bs, 2));
s_conv = ss(As,Bs,Cs,Ds);

KT_ab = T*KDT/(T); KT_ab = KT_ab(1:2, 1:2);

Piy = Ps*NCT*(1/(NDT*sqrt(3)))*KT_ab;

%%%%%%%%%%%%%%%%%%%% PARALLEL CONVERTER %%%%%%%%%%%%%%%%%%%%%
Ap = [-(Rfp/Lfp)*eye(2)  0*eye(2)       -(1/Lfp)*eye(2)  ;
       zeros(2)         -(Ry/Ly)*eye(2) -(1/Ly)*eye(2)   ;
       (1/Cfp)*eye(2)    (1/Cfp)*eye(2)  zeros(2)       ];

Bp = [(1/Lfp)*eye(2)  ;
      zeros(2)        ;
      zeros(2)       ];

Mp = eye(size(Ap, 1));
Dp = zeros(size(Mp, 1), size(Bp, 2));
p_conv = ss(Ap,Bp,Mp,Dp);

%%%%%%%%%%%%%%%%%%%%% FULL CONVERTER %%%%%%%%%%%%%%%%%%%%%%
Pvy = [zeros(2)      ;
      (1/Ly)*eye(2)  ;
      zeros(2)      ];

PiL = [ zeros(2)        ;
        zeros(2)        ;
       -(1/Cfp)*eye(2) ];

Pvg = (1/(NDT*sqrt(3)))*Pvy*KT_ab';
Pvcs = (1/(NDT*sqrt(3)))*NCT*Pvy*KT_ab';

M_a = zeros(3,9);
M_a(1:3,4:6) = eye(3);

Mp = zeros(2,6);
Mp(1:2,3:4) = eye(2);

M_b = zeros(3,6);
M_b(1:3,4:6) = eye(3);

Ms_ = zeros(2,4);
Ms_alphabeta(1:2,3:4) = eye(2);

A_complete = [As                  Piy*Mp    ;
              Pvcs*Ms_alphabeta   Ap       ];

B_complete = [Bs          zeros(4,2)    ;
              zeros(6,2)  Bp ];

C_complete = eye(size(A_complete, 1));
D_complete = zeros(size(C_complete, 1), size(B_complete, 2));

s_HDT = ss(A_complete,B_complete,C_complete,D_complete);

A = A_complete;
B = B_complete;
C = C_complete;

nx = size(A,1);
nu = size(B,2);

currentFolder = pwd;
C_Folder = extractBefore(currentFolder,'MATLAB');
C_folder = append(C_Folder,'PLECS/');
cd(C_folder)
save HDT_Variables.mat fsw Ts Vg Vs Vdc NDT NCT Lfs Rfs Cfs Cdc Lfp Rfp Cfp Ry Ly w
cd(currentFolder)