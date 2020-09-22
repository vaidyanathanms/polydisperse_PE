function n_lines = linecount(filename)
% https://www.mathworks.com/matlabcentral/answers/81137-pre-determining-the-number-of-lines-in-a-text-file
% MAY NOT BE EXACT NUMBER OF LINES. MAY BE OFF BY ONE. BUT GIVES A DECENT
% ESTIMATOR OF THE NUMBER OF LINES FOR ZEROING AN ARRAY
  fid_line = fopen(filename);
  if fid_line < 0
    error('Failed to open file %s\n', filename);
  end
    n_lines = 0;
    while true
        t = fgetl(fid_line);
        if ~ischar(t)
            break;
        else
            n_lines = n_lines + 1;
        end
    end
    fclose(fid_line);