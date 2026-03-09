function solutionstruct = testKeyTwo(A,n,d,D,fid)
%testKey Test if a specified key matrix A is admissible for all data
%        matrices
% Use: solutionstruct = testKeyTwo(A,n,d,D,fid)
% Input:
%   A is a dx(D-d) matrix for the key: R = [I_d ; A]
%   n,d,D : integers so that data matrices are n x d, and full key is d x D
%   fid : optional: File ID where to save tracing information
% Output:
%   solutionstruct : structire with the following fields:  
%     solutionstruct.nfact = n!
%     solutionstruct.failtoincrement = true when all (n!)^D permutations
%            have been checked; false if function finds a bad example
%            before reaching the end
%     solutionstruct.counterexample = number of counterexamples
%     solutionstruct.i_xipi = index of permutations (Xi1,...,Xi_(D-d),Pi1,...,Pi_d))
%     solutionstruct.badZ = basis of null space of (D-d)n x dn matrix, reshaped as n x d x r tensor
%     solutionstruct.badpiindex = PI permutation index;
%     solutionstruct.badPi = PI permutation;
%     solutionstruct.badIz = column index that achieves the maximum;
%     solutionstruct.badtestnorm = achieved norm ||PI*X-[Pi1 x1,...,Pid xd]||
%     solutionstruct.badX = bad X : n x d matrix that disproves universality;
%     solutionstruct.badW = [Pi1 x1,...,Pid xd];
% Radu Balan, 25 May 2021

if nargin < 5
    fid = 0;
end

Pi_list = perms(1:n); % Last one is the identity; it is a n! x n matrix
nfact = size(Pi_list,1); % = n!
xipi = ones(1,D); % xi are the first D-d entries; pi are the last d entries
i_xipi = 0;
counterexample = 0;
failtoincrement = false;
solutionstruct.nfact = nfact;

%parpool(40)

while ((i_xipi < nfact^(D)) && (counterexample == 0))
    qq = max(abs(xipi(D-d+1:D)-xipi(D-d+1)));
    if (qq>0)
        % non-diagonal element
        
        % Create Permutation matrices
        Xi = zeros(n,n,D-d);
        Pi = zeros(n,n,d);
        for i=1:D-d
            xivect = Pi_list(xipi(i),1:n);
            for k=1:n
                Xi(k,xivect(k),i) = 1;
            end
        end
        for i=1:d
            pivect = Pi_list(xipi(i+D-d),1:n);
            for k=1:n
                Pi(k,pivect(k),i) = 1;
            end
        end
        % Generate linear operator L : 
        % L_(i,j) block, where 1<=i<=D-d , 1<=j<=d is given by
        % L_(i,j) = A(j,i) *(Xi_i-Pi_j)
        L = zeros(n*(D-d),n*d);
        for i=1:D-d
            Xii = reshape(Xi(1:n,1:n,i),n,n);
            for j=1:d
                Pij = reshape(Pi(1:n,1:n,j),n,n);
                L((i-1)*n+1:i*n,(j-1)*n+1:j*n) = A(j,i)*(Xii-Pij);
            end
        end
        % Find null space and matrices Z1,...,Zr
        V = null(L);
        r = size(V,2);
        if (r>0)
            % Construct r matrices Z
            ZZ = zeros(n,d,r);
            for i=1:r
                ZZ(1:n,1:d,i) = reshape(V(1:n*d,i),n,d);
            end
            testnorm = zeros(nfact,1);
            iZ = zeros(nfact,1);
            parfor (i=1:nfact , 16)
                % Current permutation:
                PI = zeros(n,n);
                pivect = Pi_list(i,:);
                for j=1:n
                    PI(j,pivect(j)) = 1;
                end
                difnorm = zeros(r,1);
                for k=1:r
                    % Construct N(PI,Zk)
                    Zk = reshape(ZZ(1:n,1:d,k),n,d);
                    W = zeros(n,d);
                    for s=1:d
                        Pis = reshape(Pi(1:n,1:n,s),n,n);
                        W(1:n,s) = Pis*Zk(1:n,s);
                    end
                    difnorm(k) = norm(PI*Zk-W);
                end
                [testnorm(i),iZ(i)] = max(difnorm);
            end        
            [minnorm,imin] = min(testnorm);
            if (minnorm > 1e-6)
                % Bad matrix found
                counterexample = counterexample +1;
                bad.xipi = xipi;
                bad.piindex = imin;
                bad.Pi = Pi;
                bad.Z = ZZ;
                bad.testnorm = testnorm;
                bad.iZ = iZ;
                bad.X = reshape(ZZ(1:n,1:d,iZ(imin)),n,d);
                W = zeros(n,d);
                for s=1:d
                    Pis = reshape(Pi(1:n,1:n,s),n,n);
                    W(1:n,s) = Pis*bad.X(1:n,s);
                end
                bad.W = W;
            end
        end
    end
    % Go to next combination of permutations
    [xipi,failtoincrement] = incrementvect(xipi,D,nfact);
    i_xipi = i_xipi + 1;  
    
    % Trace: Print every 1000000 iterations of i_xipi if fid ~= 0
    if (fid ~= 0 )
        if mod(i_xipi,1000000) == 0
            fprintf(fid,'Iteration: %d\n',i_xipi);
        end
    end
end

solutionstruct.failtoincrement = failtoincrement;
solutionstruct.counterexample = counterexample;
solutionstruct.i_xipi = i_xipi;
if (counterexample > 0)
    solutionstruct.badZ = bad.Z;
    solutionstruct.badpiindex = bad.piindex;
    solutionstruct.badPi = bad.Pi;
    solutionstruct.badIz = bad.iZ;
    solutionstruct.badtestnorm = bad.testnorm;
    solutionstruct.badX = bad.X;
    solutionstruct.badW = bad.W;
else
    solutionstruct.badZ = [];
    solutionstruct.badpiindex = 0;
    solutionstruct.badPi = [];
    solutionstruct.badIz = 0;
    solutionstruct.badtestnorm = [];
    solutionstruct.badX = [];
    solutionstruct.badW = [];
end

end

