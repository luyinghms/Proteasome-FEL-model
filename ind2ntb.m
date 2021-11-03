function ret = ind2ntb(base,N,ind)
% ind can be single Nt binding state index or 
% a list of Nt binding states' index

if nargin<3
    ret = dec2base((1:1:base^N)-1,base,N);
else
    ret = dec2base(ind-1,base,N);

end