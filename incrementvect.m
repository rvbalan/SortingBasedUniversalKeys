function [newx,fail] = incrementvect(x,q,m)
%incrementvect Increment a Vector of integers modulo base
%   Detailed explanation goes here

if (q>1)
    if (x(q)<m)
        newx = x;
        newx(q) = newx(q)+1;    
        fail = false;
    else
        [subvect,fail] = incrementvect(x(1:q-1),q-1,m);
        newx = [subvect,1];
    end
else
    if (x(q)<m)
        newx = x;
        newx(q) = newx(q)+1;    
        fail = false;
    else
        newx = 1;
        fail = true;
    end    
end


end

