function [out]=asstr(inp,n)
%% n is the denominator used to compute the percentages
if nargin <2
    n=80;
end

if inp < 0 | inp > 1
        error('arcsin sqrt transform only defined for numbers between 0 and 1');
else
out=asin(sqrt(inp));
out(inp==0)=asin(sqrt(1/(4*n)));
out(inp==1)=asin(sqrt(1-1/(4*n)));
end