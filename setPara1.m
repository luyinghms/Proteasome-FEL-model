function ret=setPara1(cons,p1)
%p1(1): koff_base; p1(2): bktrf; p1(3): transEn; p1(4):tr; p1(5): kh_base

N = 6;
kon_base=1e5;%2.7e4;
kon_apo_base=1e5;
koff_base=p1(1); %**
koff_apo_base=200; %Kd_APO=2mM
conADP = cons(1);
conATPs = cons(2);
conATP = cons(3);

proj.N = N;
proj.conATP = conATP;
proj.conADP = conADP;
proj.conATPs = conATPs;
proj.con = cons;
proj.kT = 0.59; % kt=0.59kcal/mol
proj.paraVar=0.2; % give a small variation of the parameters to avoid multiple minimums
proj.bktrf=p1(2); % **energy coefficient for backward translocation. 1: conservative force, 2: frictional force

para_e.eb = [0.45 0.15 0.83 1.0 0.68 0.88];
para_e.eb = (para_e.eb-mean(para_e.eb))*proj.paraVar;% normalize
para_e.ent_atp = - 7.4 * ones(1,proj.N); %kd=exp(-dG/kT)
para_e.ent_adp = - 7.4 * ones(1,proj.N);
para_e.eif_adp = - 0.52 * ones(1,proj.N);
para_e.eif_atp = - 2.1 * ones(1,proj.N);
para_e.ent_atp_apo = - 3.7; %remove an apo Nt, no need to break the interface. 
para_e.ent_adp_apo = - 3.7;
para_e.trenref=0.0; %translocation energy reference
para_e.transEn=p1(3);% **translocation per step (2AA) causes 1.5kcal/mol energy

proj.para_e = para_e;

para_k.kon_adp = kon_base * ones(1,proj.N); %/M/sec
para_k.kon_atp = kon_base * ones(1,proj.N);
para_k.kon_atps = kon_base * ones(1,proj.N);
para_k.kon_adp_apo = kon_apo_base * ones(1,proj.N);
para_k.kon_atp_apo = kon_apo_base * ones(1,proj.N);
para_k.kon_atps_apo = kon_apo_base * ones(1,proj.N);

para_k.koff_adp = koff_base * ones(1,proj.N); %/sec
para_k.koff_atp = koff_base * ones(1,proj.N);
para_k.koff_atps = koff_base * ones(1,proj.N);
para_k.koff_adp_apo = koff_apo_base * ones(1,proj.N);
para_k.koff_atp_apo = koff_apo_base * ones(1,proj.N);
para_k.koff_atps_apo = koff_apo_base * ones(1,proj.N);
para_k.tr=p1(4); %** conformational transition rate for the energy difference in para_e.trenref

para_k.kh_base = [1, 1, 1, 1, 0, 0]*p1(5); %**
proj.para_k = para_k;

cstate.ntind=0;
cstate.ifind=0;
cstate.ntTr=[]; %nucleotide transition rate vector
cstate.cfTr=[]; % conformation transition rate
cstate.time=0;
cstate.dis=0;
cstate.ATPhy=0;
cstate.ntM=zeros(4^proj.N,4^proj.N); %nucleotide (chemical) transiiton counter
cstate.cfM=zeros(2^proj.N,2^proj.N); %conformation transition counter
cstate.ntV=zeros(1,4^proj.N); %residence time in nucleotide space
cstate.cfV=zeros(1,2^proj.N); %residence time in conformation space

proj.cstate=cstate; %current state

proj.map.ind2ntb4=strings(1,4^proj.N);
proj.map.ind2ntb3=strings(1,3^proj.N);
ind2if(1,2^proj.N)=struct('if',[],'apo',[],'hi',[],'at',[],'flag',1);
proj.map.ind2if=ind2if;

ret=proj;

end
