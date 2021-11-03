% Goals
% Calculate energy minimal conformation of each nt binding state
% Return as a cell array
% Each element contains: {
% interface: close 1 or open 0
% whether the subunit is detached: attached 1, detached 0 (apo, no interact with sub)
% free energy, Ei=Eb+Ent+Eif, calculated in MinConf
% }

function ret = MinConf(Ntb, para_e)

[num,N] = size(Ntb);
ret = cell(num,1);
NN = 2^N;
for i=1:num
    interface = dec2base((1:1:NN)-1,2,N);
    [frenergy, attach, apo] = FrEn(Ntb(i,:),interface,para_e);
    [val,~] = min(frenergy);
    notot = find(frenergy==val);  % energy degenercy
    conf.if = [];
    conf.at = [];
    conf.apo = [];
    conf.hi = [];  % find the position of each ATPaes, N+1 = detached 
    conf.energy = val;
    for h=1:length(notot)
        conf.if = [conf.if; str2vec(interface(notot(h),:))];
        conf.at = [conf.at; attach(notot(h),:)];
        conf.apo = [conf.apo; apo(notot(h),:)];
        at1 = attach(notot(h),:);
        lowk=0;
        for k=1:N
            k1=k+1;
            if k1>N 
                k1=1;
            end
            if (at1(k)==1)&&(at1(k1)==0)
                lowk=k;
                break;
            end
        end
        hi1=ones(1,N)*(N+1); 
        for k=1:N
            if at1(lowk)==0
                break;
            end
            hi1(lowk)=k;
            lowk=lowk-1;
            if lowk<1
                lowk=N;
            end
        end
        conf.hi = [conf.hi; hi1];
    end
    ret{i} = conf;
end

end
        
    