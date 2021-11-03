function [en_v, en_i]=minEn(en_vec,N)
%N: how many lowest states to output
len=length(en_vec);
D=[en_vec;1:len]';
D1=sortrows(D,1)';
en_v=D1(1,1:N);
en_i=D1(2,1:N);
end