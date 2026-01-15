clc
fsw = 20e3; Ts = 1/(2*fsw); Vg = 10e3; Vs = 400; Vdc = 700;
NDT = Vg/Vs; NCT = 5;
Lfs = 200e-6; Rfs = 100e-3; Cfs = 12e-6;
Cdc = 5000e-6;
Lfp = 200e-6; Ly = 100e-6; Ry = 5e-3; Rfp = 100e-3; Cfp = 20e-6;

w = 2*pi*50;

T = [2/3 -1/3 -1/3; 0 1/sqrt(3) -1/sqrt(3); 1/3 1/3 1/3];
KDT = [1 0 -1; -1 1 0; 0 -1 1];

%%%%%%%%%% PARALLEL CONVERTER GAMMA COMPONENT %%%%%%%%%%%%%

A_gamma = [-Rfp/Lfp  0           1/Lfp ;
            0       -Ry/Ly      -1/Ly  ;
            1/Cfp    1/Cfp       0    ];

B_gamma = [1/Lfp; 0; 0];

C_gamma = eye(3);

ss_Parallel_gamma = ss(A_gamma,B_gamma,C_gamma,[]);

nx = size(A_gamma,1);
nu = size(B_gamma,2);

currentFolder = pwd;
C_Folder = extractBefore(currentFolder,'MATLAB');
C_folder = append(C_Folder,'PLECS/');
cd(C_folder)
save HDT_Variables.mat fsw Ts Vg Vs Vdc NDT NCT Lfs Rfs Cfs Cdc Lfp Rfp Cfp Ry Ly w
cd(currentFolder)