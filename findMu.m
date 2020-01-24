function mu = findMu(t, x, mu)
%Finds the index mu such that x is in [t_mu t_mu+1).
%
%t is a row-vector of length n+p+1 and x is in [t_1, t_n+p+1)
%mu as input is a guess that may minimize searching, for instance during
%plotting.
%assumes that the interval exists

% initial tests to minimize searching during plotting 
if nargin == 3 && mu<length(t)
    if t(mu)<=x && x<t(mu+1)
        return; % mu provided in input was correct
    end
    if x<t(mu); %mu provided was too big - start over
        mu=1;
    end  %continue search from this mu
        
else % mu not provided - start from 1
    mu=1;
end

while x >= t(mu+1) && mu<length(t)-1 
    %possibility: print statement if x corresponds to last element in tvec
    mu=mu+1;
end

end %end function

