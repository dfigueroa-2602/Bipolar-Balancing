function Resonant2Ccode(Ard,Brd,mode)
% Prints C-style code for the resonant system:
%   rho_k_1 = Ard*rho_k + Brd*e_k
% And optionally prints update lines:
%   var_Control_struct.States.rho_k[i] = var_Control_struct.States.rho_k_1[i];

    if nargin < 3, mode = 0; end

    % --- Validate inputs ---
    if ~isnumeric(Ard) || ndims(Ard) ~= 2
        error('Ard must be a numeric 2-D matrix.');
    end
    if ~isnumeric(Brd) || ndims(Brd) ~= 2
        error('Brd must be a numeric 2-D matrix.');
    end

    [nRows_Ard, nCols_Ard] = size(Ard);
    [nRows_Brd, nCols_Brd] = size(Brd);

    if nRows_Ard ~= nCols_Ard
        error('Ard must be square. Got %dx%d.', nRowsA, nColsA);
    end
    if nRows_Brd ~= nRows_Ard
        error('Brd must have same number of rows as Ard. Ard is %d rows, Brd is %d rows.', nRowsA, nRowsB);
    end

    clc;
    
    [nRows_Ard, nCols_Ard] = size(Ard);
    [~, nCols_Brd] = size(Brd);
    ResonantSystem = '';
    Update = '';

    for r = 1:nRows_Ard
        line = sprintf('%s[%d] = ', 'rho_k_1', r-1);
        for c = 1:nCols_Ard
            term_r = sprintf('%s[%d][%d] * %s[%d]', 'Ard', r-1, c-1, 'rho_k', c-1);
            if c == 1
                line = [line, term_r]; %#ok<AGROW>
            else
                line = [line, ' + ', term_r]; %#ok<AGROW>
            end
        end
        for c = 1:nCols_Brd
            term_e = sprintf('%s[%d][%d] * %s[%d]', 'Brd', r-1, c-1, 'e_k', c-1);
            line = [line, ' + ', term_e]; %#ok<AGROW>
        end
        line = [line, ';', newline]; %#ok<AGROW>
        ResonantSystem  = [ResonantSystem, line]; %#ok<AGROW>
    end

    for r = 1:nRows_Ard
        line = sprintf('var_Control_struct.States.rho_k[%d] = var_Control_struct.States.rho_k_1[%d];%s',r-1, r-1, newline);
        Update = Update + string(line);
    end

    if mode == 0
        fprintf('%s', ResonantSystem);
    else
        fprintf('%s', Update);
    end
end