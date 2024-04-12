function table2latex(T, filename, column_names, sigfigs, ignore_col)
    if nargin < 2
        filename = 'table.tex';
        fprintf(['Output path is not defined. The table will be ' ...
            'written in %s.\n'], filename);
    elseif ~ischar(filename)
        error('The output file name must be a string.');
    else
        if ~strcmp(filename(end-3:end), '.tex')
            filename = [filename '.tex'];
        end
    end
    if nargin < 1, error('Not enough parameters.'); end
    if ~istable(T), error('Input must be a table.'); end
    
    % Parameters
    n_col = size(T,2);
    col_spec = [];
    for c = 1:n_col, col_spec = [col_spec 'c']; end
    col_names = strjoin(column_names, ' & ');
    row_names = T.Properties.RowNames;
    if ~isempty(row_names)
        col_spec = ['l' col_spec]; 
        col_names = ['& ' col_names];
    end
    
    % Writing header
    fileID = fopen(filename, 'w');
    fprintf(fileID, '\\begin{tabular}{%s}\n', col_spec);
    fprintf(fileID, '\\toprule\n');
    fprintf(fileID, '%s \\\\\n', col_names);
    fprintf(fileID, '\\midrule\n');
    
    % Writing the data
    for row = 1:size(T,1)
        temp{1,n_col} = [];
        for col = 1:n_col
            value = T{row,col};
            if isstruct(value)
                error('Table must not contain structs.');
            end
            while iscell(value), value = value{1,1}; end
            if isinf(value), value = '$\infty$'; end
            if ismember(col,ignore_col)
                temp{1,col} = '';
            else
                temp{1,col} = convertStringsToChars("\num{" ...
                    + sigfig(value,sigfigs(col)) + "}");
            end
        end
        if ~isempty(row_names)
            temp = [row_names{row}, temp];
        end
        fprintf(fileID, '%s \\\\ \n', strjoin(temp, ' & '));
        clear temp;
    end
    
    % Closing the file
    fprintf(fileID, '\\bottomrule\n');
    fprintf(fileID, '\\end{tabular}');
    fclose(fileID);
end