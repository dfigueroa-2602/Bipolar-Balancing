function [vel_clamp,Kap,swarm,space_range] = Swarm_Init(n_part,R_mode,xi_Mode,nQ,nR,c1,c2,Qmax,Rmax,xi_Max,rang_coef,Qvel_max,Rvel_max,xivel_max)
    % Initializes the particles
    % Definive Version
    % Added case for nx not par
    % Dave Figueroa - September 2025
    phi = c1 + c2;
    
    % Constriction factor by Clerc and Kennedy
    Kap = 2/abs(2 - phi - sqrt(phi^2 - 4*phi));
    assert(isfinite(Kap) && Kap>0,'Bad Kap');
    
    % Bounds for Q and R
    Q_blk = [-Qmax Qmax];
    Q_blks = repmat(Q_blk, nQ, 1);
    Qvel_clamp = repmat(Qvel_max,nQ,1);
    if xi_Mode == 1
        xi_blk = [0 xi_Max];
        nxi = 1;
        xivel_clamp = xivel_max;
    else
        xi_blk = [];
        nxi = 0;
        xivel_clamp = [];
    end
    if R_mode == 1
        R_blk = [-Rmax Rmax];
        R_blks = repmat(R_blk, nR, 1);
        Rvel_clamp = repmat(Rvel_max,nR,1);
    else
        nR = 0;
        R_blks = [];
        Rvel_clamp = [];
    end
    vel_clamp = [Qvel_clamp;Rvel_clamp;xivel_clamp];
    
    space_range = [Q_blks;R_blks;xi_blk];
    dim = nQ + nR + nxi;

    p_center = mean(space_range,2);
    p_range = range(space_range,2)*rang_coef;

    % Initialization of the swarm
    % (:,1,:) Position
    % (:,2,:) Velocity
    % (:,3,:) pbest
    % (:,4,:) pbest fitness
    swarm = zeros(n_part, 4, dim);

    % Compute the box center and a shrunken span to initialize positions
    for index = 1:n_part
        swarm(index,1,:) = p_center+(p_range.*(rand(dim,1) - 0.5));  %Initialization of particles positions
        swarm(index,2,:) = vel_clamp.*(rand(dim,1) - 0.5)*2;   %Initialization of particles velocities
    end
    
    % Initialization of the pbest fitness
    swarm(:,4,:) = 1e40;
end