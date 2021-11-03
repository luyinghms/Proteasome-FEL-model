function ret=outputpajek(proj,fname,para)
if nargin<3
    para(1)=1/3;
end

drawcutoff=0;
m1=sum(proj.cstate.cfM,1);
m2=sum(proj.cstate.cfM,2);
indmap=zeros(1,2^proj.N);
ind=0;
for i=1:2^proj.N
    if (m1(i)>drawcutoff)||(m2(i)>drawcutoff)
        ind=ind+1;
        indmap(i)=ind;
    end
end
    
    
fid=fopen(fname,'w');
fprintf(fid,'*Vertices  %d\n',ind);

ii=0;
for i=1:2^proj.N
    if indmap(i)>0
        ii=ii+1;
        fprintf(fid,'%d   "%d"\n',ii,i);
    end
end

%normlization factor
ML=proj.cstate.cfM(:);
MLmin=min(ML(ML>0));


fprintf(fid,'*Arcs\n');
for i=1:2^proj.N
    for j=1:2^proj.N
        if (proj.cstate.cfM(i,j)==0)||(i==j)
            continue;
        end
        sif=j;
        eif=i;
        at_j=proj.map.ind2if(sif).at;
        hi_j=proj.map.ind2if(sif).hi;
        at_i=proj.map.ind2if(eif).at;
        hi_i=proj.map.ind2if(eif).hi;
        at_all = at_j & at_i;
        dhi= mean(hi_j(at_all))-mean(hi_i(at_all)); % >0 forward; <0 backward
        if dhi==0
            cr='Green';
        else
            if dhi>0
                cr='Red';
            else
                cr='Blue';
            end
        end
          
        fprintf(fid,'%d    %d    %f   c %s\n',indmap(j),indmap(i),(proj.cstate.cfM(i,j)/MLmin)^(para(1)),cr);
    end
end

fclose(fid);
ret=proj;
end