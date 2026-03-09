% This script tests if a matrix A is a universal key

% Set up the dimensions
n = 3;
d = 4;
D = 7;
% Load data from a text file
%fid0 = fopen('inputdata.txt','r');
%fgetl(fid0); % Ignore the first line
%invec = fscanf(fid0,'%d',3);
%n = invec(1);
%d = invec(2);
%D = invec(3);
%fclose(fid0);

% Save data in text file
fid = fopen(sprintf('outputresults_n%d_d%d_D%d.txt',n,d,D),'w');


fprintf(1,'n=%d , d=%d , D=%d+%d=%d\n',n,d,d,D-d,D);
fprintf(fid,'n=%d , d=%d , D=%d+%d=%d\n',n,d,d,D-d,D);

% Generate the part-two key
rng(42);
A2 = randn(d,D-d);
A = [eye(d),A2];
fprintf(1,'Key: A = \n');
disp(A)

fprintf(fid,'Key: A = \n');
for k=1:d
    for j=1:D
        fprintf(fid,'%14.12f  ',A(k,j));
    end
    fprintf(fid,'\n');
end

% Call the main function
solutionstruct = testKeyTwo(A2,n,d,D,fid);


% Report results
if (~solutionstruct.failtoincrement)
    fprintf(1,'Stoped before all cases were explored.\n');
    fprintf(1,'Combination Number %d out of %d\n',solutionstruct.i_xipi,(solutionstruct.nfact)^(D));
    fprintf(fid,'Stoped before all cases were explored.\n');
    fprintf(fid,'Combination Number %d out of %d\n',solutionstruct.i_xipi,(solutionstruct.nfact)^(D));
end

if (solutionstruct.counterexample == 0)
    fprintf(1,'The key is universal!');
    fprintf(fid,'The key is universal!');
else
    fprintf(1,'Counterexample found!');
    fprintf(1,'X:\n');

    disp(solutionstruct.badX);
    
    fprintf(1,'W:\n');
    disp(solutionstruct.badW);

    fprintf(1,'XA:\n');
    disp(solutionstruct.badX * A);
    
    fprintf(1,'WA:\n');
    disp(solutionstruct.badW * A);

    % Save to outputfile
    fprintf(fid,'Counterexample found!\n');
    fprintf(fid,'X:\n');
    for k=1:n
        for j=1:d
            fprintf(fid,'%f ',solutionstruct.badX(k,j));
        end
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'\n W:\n');
    for k=1:n
        for j=1:d
            fprintf(fid,'%f ',solutionstruct.badW(k,j));
        end
        fprintf(fid,'\n');
    end

    fprintf(fid,'\n XA:\n');
    XA = solutionstruct.badX * A;
    for k=1:n
        for j=1:D
            fprintf(fid,'%f ',XA(k,j));
        end
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'\n WA:\n');
    WA = solutionstruct.badW * A;
    for k=1:n
        for j=1:D
            fprintf(fid,'%f ',WA(k,j));
        end
        fprintf(fid,'\n');
    end

end

fclose(fid);
