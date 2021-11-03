function ret=evolN(proj,cstate,Nsteps,todisp)
%evolve the PTSM from the starting cstate for Nsteps (chemical + conformational)
%in this version, the two target conforamtions with the lowest energy is
%considered
Nstates=2; %2
if nargin<4
    todisp=0;
end
for i=1:Nsteps
    cstate.ntTr=chemtrvec(proj,cstate.ntind,cstate.ifind);
    cftren_vec=ones(1,2^proj.N)*999999;
    cftrdis_vec=zeros(1,2^proj.N);
    for j=1:2^proj.N
        if proj.map.ind2if(j).flag==1
            [cftren_vec(j),cftrdis_vec(j)]=tren(proj,cstate.ntind,cstate.ifind,j);
        end
    end
    
    [cftren_k,ife]=minEn(cftren_vec,Nstates);
    cstate.cfTr=[];
    cstate.cfTrind=[];
    kk=0;
    for j=1:Nstates
        if cstate.ifind~=ife(j)
            cstate.cfTr=[cstate.cfTr proj.para_k.tr*exp((proj.para_e.trenref-cftren_k(j))/2.0/proj.kT)];
            cstate.cfTrind=[cstate.cfTrind ife(j)];
            kk=kk+1;
        end
    end
    vec=[cstate.ntTr cstate.cfTr];
    
    cs1=[0 cumsum(vec)];
    r=rand()*cs1(end);
    k=1;
    while ~((cs1(k)<r)&&(cs1(k+1)>r))
        k=k+1;
        if k>(4^proj.N+kk)
            break;
        end
    end
    
    if k>=(4^proj.N+1) %conformation transition
        k1=k-4^proj.N;   
        if todisp
            disp(sprintf('cf: %s(%d): %s(%d)-->%s(%d), k=%f, dis=%d',ind2ntb(4,6,cstate.ntind),cstate.ntind,ind2ntb(2,6,cstate.ifind), ...
                cstate.ifind,ind2ntb(2,6,cstate.cfTrind(k1)),cstate.cfTrind(k1),cstate.cfTr,cftrdis_vec(cstate.cfTrind(k1))));
        end
        cstate.cfM(cstate.cfTrind(k1),cstate.ifind)=cstate.cfM(cstate.cfTrind(k1),cstate.ifind)+1;
        cstate.cfV(cstate.ifind)=cstate.cfV(cstate.ifind)+1/cs1(end);
        cstate.ifind=cstate.cfTrind(k1);
        cstate.dis=cstate.dis+cftrdis_vec(cstate.cfTrind(k1));
        cstate.time=cstate.time+1/cs1(end);        
    else  %chemical transition
        if todisp
            disp(sprintf('chem: %s(%d)-->(%s)%d: %s(%d), k=%f',ind2ntb(4,6,cstate.ntind),cstate.ntind,ind2ntb(4,6,k),k, ...
                ind2ntb(2,6,cstate.ifind),cstate.ifind,cstate.ntTr(k)));
        end
        s1=char(ind2ntb(4,6,cstate.ntind))-48;
        s2=char(ind2ntb(4,6,k))-48;
        ij=find(s1~=s2);
        if (s1(ij)==3)&&(s2(ij)==1) %if atp hydrolysis
            cstate.ATPhy=cstate.ATPhy+1;
        end
        cstate.ntM(k,cstate.ntind)=cstate.ntM(k,cstate.ntind)+1;
        cstate.ntV(cstate.ntind)=cstate.ntV(cstate.ntind)+1/cs1(end);
        cstate.ntind=k;
        cstate.time=cstate.time+1/cs1(end);
    end
end
proj.cstate=cstate;
ret=proj;
end
