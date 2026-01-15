
function Kgain2Ccode(Kcell, mat_names, h_filename)
% Kgain2Ccode({K1, K2, ...}, {'MACRO1','MACRO2',...}, 'Gain')
% Writes one header Gain.h with multiple #define blocks.

    % --- sanity checks ---
    assert(iscell(Kcell),      'First argument must be a cell array of matrices');
    assert(iscell(mat_names),  'Second argument must be a cell array of names');
    assert(numel(Kcell) == numel(mat_names), ...
        'Kcell and mat_names must have the same number of elements');

    nMat = numel(Kcell);

    % --- build text for all matrices ---
    all_blocks = '';

    for idx = 1:nMat
        Ki    = Kcell{idx};
        namei = mat_names{idx};

        if ~ismatrix(Ki) || ~isnumeric(Ki)
            error('Element %d of Kcell is not a numeric 2-D matrix.', idx);
        end

        % Convert matrix to row strings
        Ki_str = mat2str(Ki);
        Ki_str = extractBetween(Ki_str, '[', ']');
        Ki_str = split(Ki_str, ';');   % one cell per row

        nRows = size(Ki,1);

        rows = cell(nRows,1);
        for r = 1:nRows
            row_str    = Ki_str{r};
            row_str    = strrep(row_str, ' ', ', ');
            rows{r}    = ['{ ' row_str ' }'];
        end

        % Build #define block for this matrix
        block = ['#define ' namei ' {' '   \\' '\n'];
        for r = 1:nRows
            if r < nRows
                block = [block '    ' rows{r} ', ' '\\' '\n'];  %#ok<AGROW>
            else
                block = [block '    ' rows{r} ' \\' '\n'];  %#ok<AGROW>
            end
        end
        block = [block '};' '\n\n'];    %#ok<AGROW>

        all_blocks = [all_blocks block];   %#ok<AGROW>
    end
        CurrentFolder = pwd;
        C_Folder = extractBefore(CurrentFolder,'MATLAB');
        C_folder = append(C_Folder,'C_code/src/');
        C_file = fopen([C_folder h_filename '.h'], 'wt');
        aux = ['#pragma once\n#ifndef ' h_filename '_H_\n#define ' h_filename '_H_\n\n'];
        fprintf(C_file, aux);
        fprintf(C_file, [all_blocks '\n\n']);
        aux = ['#endif /* ' h_filename '_H_ */'];
        fprintf(C_file, aux); fclose(C_file);
        disp(['Your(s) matrix/matrices gain(s) is/are ready in the ' C_folder ' folder'])
end