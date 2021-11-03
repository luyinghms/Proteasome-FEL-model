function ret=bdmap(proj)
for i=1:4^proj.N
    proj.map.ind2ntb4(i)=ind2ntb(4,6,i);
end
for i=1:3^proj.N
    proj.map.ind2ntb3(i)=ind2ntb(3,6,i);
end
a=proj.map.ind2if;
for i=1:2^proj.N
    a(i).if=ind2ntb(2,6,i)-48;
    if_vec=a(i).if;
    % General rules
% # of open interface >=2, <=4
    temp = nnz(if_vec==0);  %number of open interfaces
    if temp<2 || temp>4
        a(i).flag=0;
        continue;
    end            
 % calculating the attachment state given an interface 
    N=proj.N;
    att_f1=ones(1,N);
    for j=1:(N-1)
        if if_vec(j)==1
            att_f1(j+1)=att_f1(j);
        else
            att_f1(j+1)=att_f1(j)+1;
        end
    end
    if if_vec(N)==1
        att_f1(att_f1==att_f1(N))=1;
    end
    att_f1_f=hist(att_f1,1:max(att_f1));
    [att_f1_f_max, att_f1_f_maxi]=max(att_f1_f);
    if att_f1_f_max<3  %if the largest domain is <3 exclude
        a(i).flag=0;
        continue;
    else
        if (att_f1_f_max == 3) && (temp == 2) % if 3 + 3, exclude
            a(i).flag=0;
            continue;
        else            
            attach1(1:N)=0;
            attach1(att_f1==att_f1_f_maxi)=1;
            attach(i,:)=attach1;
            a(i).at=attach1;
            a(i).apo=~attach1;
            %height calculation
            lowk=0;
            for k=1:N
                k1=k+1;
                if k1>N 
                    k1=1;
                end
                if (attach1(k)==1)&&(attach1(k1)==0)
                    lowk=k;
                    break;
                end
            end
            hi1=ones(1,N)*(N+1); 
            for k=1:N
                if attach1(lowk)==0
                    break;
                end
                hi1(lowk)=k;
                lowk=lowk-1;
                if lowk<1
                    lowk=N;
                end
            end
            a(i).hi = hi1;            
            for k=1:N  %this allows close non-apo states in detached ATPases if the interface is closed.
                if (attach1(k)==0)&&(if_vec(k)==1)
                    a(i).apo(k)=0;
                end
            end
            a(i).flag=1;
        end
    end    
end
proj.map.ind2if=a;

ret=proj;
end