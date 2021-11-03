function ret = ntb2ind(base,ntb)
% ntb can be single Nt binding state or 
% a list of Nt binding states

if strcmp(class(ntb),'char')
    ret = base2dec(ntb,base) + 1;
else
    ret = base2dec(char(ntb+48),base) +1;
end

end