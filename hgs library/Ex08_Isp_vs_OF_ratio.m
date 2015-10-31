%***********************************************************************************************************
%* HGS 1.3 
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 08: Compute the ISP of H2, O2 mixture vs OF ratio

function Ex08_Isp_vs_OF_ratio

close all;

species={'H','H2','H2O','H2O2','HO2','O','O2','OH'};

% Inlet temperature as if the reactives were gas at 300K
Te=300 % K   reactives inlet temperature
Pc=50  % bar chamber pressure 
P2=0.1 % bar nozzle exit    
    
vrof=[]; % vector of OF ratios
visp=[]; % vector of specific impulses



Tcs=1800;
T2s=400;

for rof=2:0.25:8
    
    fprintf('solving for rof=%f \n',rof);
   
    % evaluate mol of each specie at inlet for a given ROF ratio
    nO2=1;
    mO2=nO2*32;
    mH2=mO2/rof;
    nH2=mH2/2;
    
    ni_i=[0;... % H
    nH2;... % H2
    0;...   % H2O
    0;...   % H2O2
    0;...   % HO2
    0;...   % O
    nO2;... % O2
    0];     % OH

    ni_i=ni_i/sum(ni_i); % mole fractions  
    

    % Evaluate inlet properties with HGS assuming gas state
    [~,~,MM,~,~,~,H,~,~]=hgsprop(species,ni_i,Te,Pc);
    n=sum(ni_i); % mixture total number of mols (1)
    m=n*MM*1e-3; % mixture mass kg
    h1G=H/m;     % inlet mixture enthalpy in GAS state kJ/kgK 
                 % we evaluate it just for comparision

    % Inlet enthalpy as if reactives were satured liquid at 10 bar
    % O2 (NIST) hv(404.36 K)-hl(119.62 K)=14.3753 kJ/mol
    % H2 (NIST) hv(413.96K)-hl(31.39K)=10.9495 kJ/mol

    % Enthalpy of O2 liq at Tsat 10 bar (kJ/mol)
    hO2=hgssingle('O2','h',404.36,10)-14.3753; 
    % Enthalpy of H2 liq at Tsat 10 bar (kJ/mol)
    hH2=hgssingle('H2','h',413.96,10)-10.9495; 
    Hin=ni_i(2)*hH2+ni_i(7)*hO2;
    h1=Hin/m; % inlet mixture enthalpy in LIQUID state kJ/kg 

    % We find temperature at nozzle inlet solving for Delta_H=0
    % hgsTp function can't be used as it assumes gas state
    
    Tc=fsolve(@DeltaH,Tcs,optimset('Display','none'));
    fprintf('Chamber outlet temperature Tc=%f \n',Tc);
    Tcs=Tc; % In next interation, we will begin with the solution obtained now
    
    %Tc=myfzero(@DeltaH,300,4500,10,1e-4,1e-5,300);
    ni_calc=hgseq(species,ni_i,Tc,Pc);
    [Cp,Cv,MM,Rg,gamma,a,H,G,S]=hgsprop(species,ni_calc,Tc,Pc);
    m=sum(ni_calc)*MM*1e-3; % mixture mass kg (has to be as before!)
    s=S/m;    
        
    % We use the previous value of T2 to begin the iterations
    % fzero is more robust solver in this case, but decreasing the
    % tolerance we can solve the problem with fsolve, that is faster
    [T2,~,vt ] = hgsisentropic(species,ni_calc,Tc,Pc,P2,'shifting','fsolve',T2s,optimset('Display','iter','TolFun',1e-9,'TolX',1e-3));
    T2s=T2;
    
    fprintf('Nozzle outlet temperature T2=%f \n',T2);
    
    Is=vt/9.81; % Is (optimal expansion, Pe=Pambient)    

    vrof(end+1)=rof;
    visp(end+1)=Is;
    
    plot(vrof,visp); xlabel('ratio O/F'); ylabel('Isp (s)');
    drawnow;
end



    function DeltaH=DeltaH(Tc)
        nc=hgseq(species,ni_i,Tc,Pc);
        [~,~,MMC,~,~,~,HC,~,~]=hgsprop(species,nc,Tc,Pc); 
        nC=sum(nc); % mixture total number of mols (1)
        mc=nC*MMC*1e-3; % mixture mass kg
        hc=HC/mc; % kJ/kgK
        DeltaH=hc-h1;
    end
    
end

