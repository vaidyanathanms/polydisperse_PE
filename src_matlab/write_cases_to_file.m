function write_cases_to_file(fylename)

fin_cons = fopen(fylename,'r');

if fin_cons <= 0 % check for average list
    fprintf('%s does not exist', fylename);
    continue;
end