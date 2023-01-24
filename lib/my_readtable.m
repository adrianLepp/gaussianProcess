function T = my_readtable(filename,numRecords) %#codegen
    f = fopen(filename,"r");
    % THIS is os ugly, but for the start it works...
    
% • The size of the character vector line changes as the function reads each line of the initial
%   description. The coder.varsize directive instructs the code generator to produce code that
%   dynamically allocates memory for the variable line.
% • The function hardcodes the column header names 'time' and 'temp' instead of reading them
%   from the CSV file at run time. This is because code generation requires table variable names to be
%   compile-time constants.

    % Scan and ignore a variable number of description lines at the beginning
    % of the CSV file.
    line = fgetl(f);
    coder.varsize("line");
    while(~ismember(';',line))
        line = fgetl(f);
    end
    
    % Table variable names.
    names = {'time','x1','x2','x3','dx1','dx2','dx3','y1','y2','y3'};
    i = 1;
    
    time = zeros(1,numRecords);
    x1 = zeros(1,numRecords);
    x2 = zeros(1,numRecords);
    x3 = zeros(1,numRecords);
    dx1 = zeros(1,numRecords);
    dx2 = zeros(1,numRecords);
    dx3 = zeros(1,numRecords);
    y1 = zeros(1,numRecords);
    y2 = zeros(1,numRecords);
    y3 = zeros(1,numRecords);
    
    % Read each line in the CSV file till you reach EOF
    while(~feof(f) )
        %%line = fgetl(f);
        % read until semicolon%[^;
        
        [result,count] = fscanf(f,'%f;%f;%f;%f;%f;%f;%f;%f;%f;%f',10);
        assert(count == 10) 
        time(i) = result(1);
        x1(i) = result(2);
        x2(i) = result(3);
        x3(i) = result(4);
        dx1(i) = result(5);
        dx2(i) = result(6);
        dx3(i) = result(7);
        y1(i) = result(8);
        y2(i) = result(9);
        y3(i) = result(10);
      
        if (i == 180) 
            break;
        end
        i = i + 1;
    end
    
    % Construct the table from the values read in the previous code block.
    T = table(time',x1',x2',x3',dx1',dx2',dx3',y1',y2',y3','VariableNames',names);
    fclose(f);
end
