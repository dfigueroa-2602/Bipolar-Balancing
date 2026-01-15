function fitness = LQR_Search(Res_mode,Ts,n_h,w_vec,R_mode,particle,...
                              A_aug,B_aug,Ad,Bd,Ard,Brd,Ha,...
                              beta_c,VDC,w,A_ref,Tsim,T_settling)
% Function that injects the given particle into the lqr function, then call
% the simulation function.
% Version 0.4
%   Added the option for odd number of states
% Dave Figueroa
% January 2026

    % Plant dimensions
    nx = size(Bd,1); nu = size(Bd,2);

    if Res_mode == 1
        % Resonant case
        nxr = size(Brd,1);
        states_per_h = nxr / n_h;
        if mod(nx,2) == 0
            % Split particle p
            idx_x_end = nx/2; idx_u_end = idx_x_end + nu/2;
        else
            % Use the complete particle
            idx_x_end = nx; idx_u_end = idx_x_end + nu;
        end

        p_x  = particle(1:idx_x_end);                  % Plant states (half)
        p_u  = particle(idx_x_end+1:idx_u_end);        % Plant inputs (half)
        p_hr = particle(idx_u_end+1:idx_u_end+n_h);    % 1 scalar per harmonic
        
        if mod(nx,2) == 0
            % Duplication of each particle to construct the Q matrix
            Qx = repelem(p_x,2); Qu = repelem(p_u,2);
            if isrow(Qu), Qu = Qu'; end
        else
            Qx = p_x; Qu = p_u;
            if isrow(Qu), Qu = Qu'; end
        end

        % Construction of the resonant Q with Bryson scaling ----
        Qr_vec = zeros(nxr,1);
        offset = 0;
        for i = 1:n_h
            qi  = 10.^p_hr(i);
            wi2 = w_vec(i)^2;
            Qr_block = qi * wi2 * ones(states_per_h,1);
            Qr_vec(offset+1:offset+states_per_h) = Qr_block;
            offset = offset + states_per_h;
        end

        Q_diag = [10.^Qx; 10.^Qu; Qr_vec]; % Assembly of the Q matrix
        Q = diag(Q_diag);

        % LQR on full augmented system
        A_lqr = A_aug; B_lqr = B_aug;

    else
        % Non resonant case. With additional states for integrators
        nL = size(A_aug,1);
        assert(nL == nx + 2*nu, 'Expected augmented size nx+2*nu');

        idx_x_end = nx; idx_u_end = nx + nu; idx_i_end = nx + 2*nu;

        p_x = particle(1:idx_x_end);
        p_u = particle(idx_x_end+1:idx_u_end);
        p_i = particle(idx_u_end+1:idx_i_end);

        if isrow(p_u), p_u = p_u'; end
    
        Q_diag = [10.^p_x;10.^p_u;10.^p_i]; % Assembly of the Q matrix
        Q = diag(Q_diag);
    
        A_lqr = A_aug; B_lqr = B_aug; Ard = []; Brd = [];
    end

    if R_mode == 1
        % Case in which also the R weights are searched
        Ru = repelem(particle(end - nu/2 + 1:end),2);
        R  = diag((10.^Ru)');
    else
        % R is identity
        R = eye(nu);
    end

    if any(R(:) < 0) || any(Q(:) < 0)
        disp('One of the elements of the Q and R matrixes is negative!')
        fitness = 1e6; return;
    end

    try
        [K,~,~] = dlqr(A_lqr, B_lqr, Q, R);
    catch
        disp('Error in the LQR')
        fitness = 1e6; return;
    end

    Kx = K(:,1:idx_x_end); Ku = K(:,idx_x_end+1:idx_u_end);

    if Res_mode == 1
        Kr = K(:,idx_u_end+1:end);
    else
        Kr = zeros(size(K,1),0); Ki = K(:,idx_u_end+1:idx_i_end);
    end

    % ----------------------------
    % Stability check (on the system actually used in LQR)
    % ----------------------------
    Acl  = A_lqr - B_lqr*K;
    spec = max(abs(eig(Acl)));
    if spec >= 1
        fitness = 1e6;
        disp('One of the particles is too close of the instable area!')
        return;
    end

    try
        if Res_mode == 1
            fitness = Sim_ResSystem(Tsim,Ts,VDC,w,A_ref,T_settling,beta_c,...
                                Ad,Bd,Ard,Brd,Ha,Kx,Kr,Ku);
        else
            fitness = Sim_AugSystem(Tsim,Ts,VDC,w,A_ref,T_settling,beta_c,...
                         Ad,Bd,Ha,nx,nu,Kx,Ku,Ki);
        end
    catch ME
        fprintf('\nSimulation failed: %s\n', ME.message)
        disp(getReport(ME, 'extended'))
        fitness = 1e6;
        return;
    end
end
