function cost = Sim_AugSystem(Tsim,Ts,VDC,w,A_ref,t_settling,beta_c,...
                               Ad,Bd,Ha,nx,nu,Kx,Ku,Ki)
% Sim_System_aug
% Simulates the plant:
%   x_k_1 = Ad*x_k + Bd*u_k
% In which, delays are considered.
% Control law:
%   u_k = -Kx*x_k -Ku*xu_k -Ki*xi_k
% where:
% xi_k_1 = xi_k + Ts*(r_k - Ha*x_k)
% Version 0.1
% Dave Figueroa
% January 2026

    V_ref_amp = A_ref(1);

    % Samples
    N = round(Tsim / Ts);

    % Time vector
    t = (0:N-1).' * Ts;

    % States
    x_k  = zeros(nx,1); xu_k = zeros(nu,1); xi_k = zeros(nu,1);
    u_prev = zeros(nu,1);

    % Logs
    e_log  = zeros(N,nu); du_log = zeros(N,nu);

    % Reference generator (supports scalar or nu-vector reference)
    t_step = 0.1; tau = t_settling/4; env = 1 - exp(-t/tau);

    if isscalar(A_ref)
        Av_vec = zeros(N,1);
        Av_vec(t >= t_step) = V_ref_amp;
        r = (Av_vec .* env .* sin(w*t));           % Nx1
    else
        % if you ever pass per-channel amplitudes in A_ref (nu x 1)
        Av = A_ref(:).';
        Av_mat = repmat(Av, N, 1);                % Nxnu
        r = (env .* sin(w*t)) .* Av_mat;          % Nxnu
    end

    for k = 1:N
        % Output to be tracked
        z_k = Ha*x_k;                 % scalar or vector
        if isscalar(z_k), z_k = repmat(z_k, nu, 1); else, z_k = z_k(:); end

        % Reference at k
        if size(r,2) == 1, r_k = repmat(r(k), nu, 1); else, r_k = r(k,:).'; end

        % Error calculation
        e_k = r_k - z_k;
        e_log(k,:) = e_k.';

        % Integration of the error
        xi_k_1 = xi_k + Ts*e_k;

        % Control law
        u_k = -Kx*x_k - Ku*xu_k - Ki*xi_k;

        % Saturation
        usat_k = max(min(u_k, VDC), -VDC);

        % State update according to plant dynamics
        x_k_1 = Ad*x_k + Bd*xu_k;

        % du (use saturated)
        du_k = u_k - u_prev;
        du_log(k,:) = du_k.';
        u_prev = usat_k;
        
        % States update
        x_k = x_k_1;
        xu_k = usat_k;
        xi_k = xi_k_1;
    end

    term_e  = sum(e_log.^2, 2);
    term_du = sum(du_log.^2, 2);
    cost = (1/N) * sum(term_e + beta_c*term_du);
end
