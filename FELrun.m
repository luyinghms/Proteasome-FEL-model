function ret=FELrun(p1,nt1)
%p1=model parameters: the default = [0.3333   -0.2576    0.1984   72.8000    1.0340]
%nucleotide condition (nM) [ATP-gS, ADP, ATP]
if nargin<2
    sz=size(p1);
    nt1=zeros(sz(1),3);
    for i=1:sz(1)
        nt1(i,:)=[0 0 500];
    end
end

N=6; %# of ATPase
Nsteps=100000; %simulate Nsteps
Nrep=40; %repeat Nrep times
sz=size(p1);

cstate.ntind=3408; %index of the initial ATPase conformation
cstate.ifind=52; %index of the initial Nt staus
cstate.ntTr=[]; %nucleotide transition rate vector
cstate.cfTr=0; % conformation transition rate
cstate.time=0;
cstate.dis=0;
cstate.ATPhy=0;
cstate.ntM=zeros(4^N,4^N); %nucleotide (chemical) transiiton counter
cstate.cfM=zeros(2^N,2^N); %conformation transition counter
cstate.ntV=zeros(1,4^N); %residence time in nucleotide space
cstate.cfV=zeros(1,2^N); %residence time in conformation space

ret=cell(1,sz(1));

for i=1:sz(1)
    proj1=setPara1(nt1(i,:),p1(i,:));
    proj1=bdmap(proj1);
    for k=1:Nrep
        proj1.cstate.ntind=cstate.ntind;
        proj1.cstate.ifind=cstate.ifind; %resetting starting point
        proj1=evolN(proj1,proj1.cstate,Nsteps);
    end
    disp('finished test...');
    ret{i}=proj1;
end
end