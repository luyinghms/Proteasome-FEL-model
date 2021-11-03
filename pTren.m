function ret=pTren()
p1=[1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.1, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1;
    1/3, -0.2, 0.3, 120, 2.1];
nt1=[0 0 500;
    0 0 40;
    0 0 20000;
    800 0 500;
    4000 0 500;
    0 36 500;
    0 200 500;
    0 0 500;
    0 0 40;
    0 0 20000;
    800 0 500;
    4000 0 500;
    0 36 500;
    0 200 500];

N=6;
Nsteps=10000; %simulate Nsteps
Nrep=5; %repeat Nrep times
sz=size(p1);

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

res=cell(1,sz(1));

parfor i=1:sz(1)
    proj1=setPara2(nt1(i,:),p1(i,:));
    proj1=bdmap(proj1);
    for k=1:Nrep
        proj1.cstate.ntind=cstate.ntind;
        proj1.cstate.ifind=cstate.ifind; %resetting starting point
        proj1=evolN(proj1,proj1.cstate,Nsteps);
    end
    disp('finished test...');
    res{i}=proj1;
end

x=[];
y1=[];
y2=[];
y3=[];
for i=1:sz(1)
    y1=[y1 res{i}.cstate.time];
    y2=[y2 res{i}.cstate.dis];
    y3=[y3 res{i}.cstate.ATPhy];
    x=[x res{i}.para_e.transEn];
end

ret=[x' y1' y2' y3'];