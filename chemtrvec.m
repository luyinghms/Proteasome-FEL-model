function ret=chemtrvec(proj,inds,ifs)
%calculate the chemmical transition vector based on the current Nt state
%and interface state
nts=char(proj.map.ind2ntb4(inds))-48;
apo=proj.map.ind2if(ifs).apo;
att=proj.map.ind2if(ifs).at;
hi=proj.map.ind2if(ifs).hi;
ifvec=proj.map.ind2if(ifs).if;
vec=zeros(1,4^proj.N);
for i=1:proj.N
    switch nts(i)
        case 0 %can bind ADP, gs or ATP
            nts1=nts; nts1(i)=1; ind1=ntb2ind(4,char(nts1+48));
            nts2=nts; nts2(i)=2; ind2=ntb2ind(4,char(nts2+48));
            nts3=nts; nts3(i)=3; ind3=ntb2ind(4,char(nts3+48));
            if apo(i)==1
                vec(ind1)=proj.conADP/1e6*proj.para_k.kon_adp_apo(i);
                vec(ind2)=proj.conATPs/1e6*proj.para_k.kon_atps_apo(i);
                vec(ind3)=proj.conATP/1e6*proj.para_k.kon_atp_apo(i);
            else
                vec(ind1)=proj.conADP/1e6*proj.para_k.kon_adp(i);
                vec(ind2)=proj.conATPs/1e6*proj.para_k.kon_atps(i);
                vec(ind3)=proj.conATP/1e6*proj.para_k.kon_atp(i);                
            end         
        case 1 %ADP off
            nts1=nts; nts1(i)=0; ind1=ntb2ind(4,char(nts1+48));
            if apo(i)==1
                vec(ind1)=proj.para_k.koff_adp_apo(i);
            else
                if ifvec(i)==1  %if closed interface
                    vec(ind1)=proj.para_k.koff_adp(i)*exp(proj.para_e.eif_adp(i)/proj.kT);
                else
                    vec(ind1)=proj.para_k.koff_adp(i);
                end
            end
        case 2 %ATPgs off
            nts1=nts; nts1(i)=0; ind1=ntb2ind(4,char(nts1+48));
            if apo(i)==1
                vec(ind1)=proj.para_k.koff_atps_apo(i);
            else
                if ifvec(i)==1                    
                    vec(ind1)=proj.para_k.koff_atps(i)*exp(proj.para_e.eif_atp(i)/proj.kT);
                else
                    vec(ind1)=proj.para_k.koff_atps(i);
                end
            end
        case 3 %ATP off or hygrolysis (attached state only)
            nts1=nts; nts1(i)=0; ind1=ntb2ind(4,char(nts1+48));
            nts2=nts; nts2(i)=1; ind2=ntb2ind(4,char(nts2+48));
            if att(i)==0
                vec(ind1)=proj.para_k.koff_atp_apo(i);
            else
                if ifvec(i)==1
                    vec(ind1)=proj.para_k.koff_atp(i)*exp(proj.para_e.eif_atp(i)/proj.kT);
                else
                    vec(ind1)=proj.para_k.koff_atp(i);
                end
                
                if hi(i)<=proj.N
                    vec(ind2)=proj.para_k.kh_base(hi(i)); %kh may depend on distance to the 20S proteasome
                    if ifvec(i)==0    %if open interface, kh=0 due to absence of Arg-finger 
                        vec(ind2)=0;
                    end
                    if isfield(proj,'WBM')
                        if (i==proj.WBM)
                            vec(ind2)=0;
                        end
                    end
                end
            end
    end
    
end
ret=vec;

end