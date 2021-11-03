function ret=Ntscan(Nt,p1)

if nargin<2
    p1=[1/3, -0.5, 0.75, 150, 5];
    if nargin<1
        NtATP=[0 0 3;0 0 15;0 0 40;0 0 150; 0 0 500; 0 0 2000; 0 0 20000];
        NtADP=[0 0 500; 250 0 500; 500 0 500; 1000 0 500; 2000 0 500; 4000 0 500; 6000 0 500];
        NtgS=[0 0 500; 0 15 500; 0 30 500; 0 45 500; 0 60 500; 0 90 500; 0 120 500];
        Nt=[NtATP;NtADP;NtgS];
    end
end

N=6;
Nsteps=10000; %simulate Nsteps
Nrep=5; %repeat Nrep times

cstate.ntind=3408;
cstate.ifind=52;
cstate.ntTr=[]; %nucleotide transition rate vector
cstate.cfTr=0; % conformation transition rate
cstate.time=0;
cstate.dis=0;
cstate.ATPhy=0;
cstate.ntM=zeros(4^N,4^N); %nucleotide (chemical) transiiton counter
cstate.cfM=zeros(2^N,2^N); %conformation transition counter
cstate.ntV=zeros(1,4^N); %residence time in nucleotide space
cstate.cfV=zeros(1,2^N); %residence time in conformation space

[len,~]=size(Nt);
ret=cell(1,len);
parfor i=1:len    
    proj1=setPara1(Nt(i,:),p1);
    proj1=bdmap(proj1);
    for k=1:Nrep
        proj1.cstate.ntind=cstate.ntind;
        proj1.cstate.ifind=cstate.ifind; %resetting starting point
        proj1=evolN(proj1,proj1.cstate,Nsteps);
    end
    disp('finished test 1');
    ret{i}=proj1;
end

end
