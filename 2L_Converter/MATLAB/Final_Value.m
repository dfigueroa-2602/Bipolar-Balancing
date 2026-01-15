function [Kx, Kr, Ku, Ki] = Final_Value(Res_mode,b_swarm,w,w_vec,A_aug,B_aug,nx,nu,nxr,n_h,R_mode)
% Final_Value
% Returns gains from the best p_best for EITHER:
%   (A) Resonant augmented system  x_aug = [x; xu; xr]      size = nx + nu + nxr
%   (B) Non-resonant augmented     x_aug = [x; xu; xi]      size = nx + nu + nu = nx + 2*nu
%
% Assumes the same swarm packing you used in LQR_Search():
%   - Resonant:  [ sw_x( nx/2 ); sw_u( nu/2 ); sw_hr( n_h ); (optional R terms) ]
%               (duplicated for +-)
%   - Non-res:   [ sw_x( nx );   sw_u( 2*nu ); (optional R terms) ]
%
% Version 0.4
% Dave Figueroa
% January 2026

    p_best = squeeze(b_swarm(:));  % force column vector

    if Res_mode == 1
        % Resonant case
        states_per_h = nxr / n_h;
        if mod(nx,2) == 0
            % Split best particle
            idx_x_end = nx/2; idx_u_end = idx_x_end + nu/2;
        else
            % Use the complete particle
            idx_x_end = nx; idx_u_end = idx_x_end + nu;
        end

        p_x  = p_best(1:idx_x_end);                  % Plant states (half)
        p_u  = p_best(idx_x_end+1:idx_u_end);        % Plant inputs (half)
        p_hr = p_best(idx_u_end+1:idx_u_end+n_h);    % 1 scalar per harmonic

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

        p_x = p_best(1:idx_x_end);
        p_u = p_best(idx_x_end+1:idx_u_end);
        p_i = p_best(idx_u_end+1:idx_i_end);

        if isrow(p_u), p_u = p_u'; end
    
        Q_diag = [10.^p_x;10.^p_u;10.^p_i]; % Assembly of the Q matrix
        Q = diag(Q_diag);
    
        A_lqr = A_aug; B_lqr = B_aug;
    end

    if R_mode == 1
        % Case in which also the R weights are searched
        Ru = repelem(p_best(end - nu/2 + 1:end),2);
        R  = diag((10.^Ru)');
    else
        % R is identity
        R = eye(nu);
    end

    [K,~,~] = dlqr(A_lqr, B_lqr, Q, R);

    % Split K = [Kx Ku Ki]
    Kx = K(:, 1:nx);
    Ku = K(:, nx+1 : nx+nu);     % gain on xu (input-state)
    if Res_mode == 1
        Kr = K(:, nx+nu+1 : end);
        Ki = [];
    else
        Ki = K(:, nx+nu+1 : nx+2*nu);
        Kr = [];
    end
end
