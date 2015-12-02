clear all; clc
n = [2,1,0,0,0,0];
T=2700;
P=1;

D = {'H2','O2','H2O','H','O','OH'};

id = zeros(length(D),1);

for i=1:length(D)
    id(i) = hgsid(D(i));
    Data(i)=hgsDB(id(i));
end
    
Out=hgsmix(Data,'s',T,P,n)

n2 = [0.1596, 0.056, 1.7721, 0.0332, 0.0125, 0.1033];
Out2=hgsmix(Data,'s',T,P,n2);

[neq,deltaG] = hgseq(D,n,T,P);

[Cp,Cv,MM,Rg,gamma,a,H,G,S] = hgsprop(D,n,T,P)
