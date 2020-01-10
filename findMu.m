function mu = findMu(tvec, x, mu)
%Finds the index mu such that x is in [tvec_mu tvec_mu+1].
%(Compendium: Exercise 2.11)
%tvec is a row-vector of length n+p+1 and x is in [tvec_1, tvec_n+p+1)
%mu as input is a guess that may minimize searching, for instance during
%plotting.
%assumes that the interval exists

% initial tests to minimize searching during plotting 
if nargin == 3
    if tvec(mu)<=x && x<tvec(mu+1)
        return; % mu provided in input was correct
    end
    if x<tvec(mu); %mu provided was too big - start over
        mu=1;
    end  %continue search from this mu
        
else % mu not provided - start from 1
    mu=1;
end

while x >= tvec(mu+1) && mu<length(tvec)-1 
    %possibility: print statement if x corresponds to last element in tvec
    mu=mu+1;
end

end %end function

