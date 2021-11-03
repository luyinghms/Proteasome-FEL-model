function [ret,dhi]=tren(proj,ntbind,sif,eif)
%calulating the energy accompanying the conformational change from sif to
%eif, given the nucleotide-binding index ntbind, and parameters in proj

%starting energy
if proj.map.ind2if(sif).flag * proj.map.ind2if(eif).flag ==0
     ret=99999;
else 
    if_j=proj.map.ind2if(sif).if;
    apo_j=proj.map.ind2if(sif).apo;
    i2_vec=char(proj.map.ind2ntb4(ntbind))-48;
    Ebs = sum(if_j.*proj.para_e.eb);
    Ents = sum((i2_vec==1 & ~apo_j).*proj.para_e.ent_adp) + sum((i2_vec>1 & ~apo_j).*proj.para_e.ent_atp) - sum((i2_vec==1 & apo_j).*proj.para_e.ent_adp_apo) - sum((i2_vec>1 & apo_j).*proj.para_e.ent_atp_apo);
    Einfs = sum((i2_vec==1 & if_j).*proj.para_e.eif_adp) + sum((i2_vec>1 & if_j).*proj.para_e.eif_atp);
    ens=Ebs+Ents+Einfs; %<0: can drive translocation, >=0, cannot drive translocation
    
    if_j=proj.map.ind2if(eif).if;
    apo_j=proj.map.ind2if(eif).apo;
    Ebe = sum(if_j.*proj.para_e.eb);
    Ente = sum((i2_vec==1 & ~apo_j).*proj.para_e.ent_adp) + sum((i2_vec>1 & ~apo_j).*proj.para_e.ent_atp) - sum((i2_vec==1 & apo_j).*proj.para_e.ent_adp_apo) - sum((i2_vec>1 & apo_j).*proj.para_e.ent_atp_apo);
    Einfe = sum((i2_vec==1 & if_j).*proj.para_e.eif_adp) + sum((i2_vec>1 & if_j).*proj.para_e.eif_atp);
    ene=Ebe+Ente+Einfe; %<0: can drive translocation, >=0, cannot drive translocation
    
    %if Lid is taken into consideration
    if isfield(proj,'LidID')  %LidID is the interface id for association with Lid
        if proj.LidID==sif
            ens=ens+proj.LidEN;
        end
        if proj.LidID==eif
            ene=ene+proj.LidEN;
        end
    end
    
    at_j=proj.map.ind2if(sif).at;
    hi_j=proj.map.ind2if(sif).hi;
    at_i=proj.map.ind2if(eif).at;
    hi_i=proj.map.ind2if(eif).hi;
    at_all = at_j & at_i;
    dhi= mean(hi_j(at_all))-mean(hi_i(at_all)); % >0 forward; <0 backward
    if isempty(dhi)
        ret=99999;
    else
        if dhi>=0
            ret=ene-ens+dhi*proj.para_e.transEn; %<0: favorable; >0: unfavorable
        else
            ret=ene-ens+dhi*proj.para_e.transEn*proj.bktrf;
        end
    end
end

end