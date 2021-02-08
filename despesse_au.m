% Clear the command window
clc;

% Phyiscal variables
E=([-30,30]);% minumum and maximum for energy scale in eV
PhiH=1;        % hot electrode workfunction in eV
PhiC=2;        % cold electrode workfunction in eV
Vbias=-0.5;     % Appled external voltage in volts
d=5e-10;       % vacuum gap in metres
TH=10000; %Hot electrode temperature (K)
TC=3000; %Cold electrode temperature (K)

% Set the physical constants in SI units
e=1.6e-19;     % Electron charge in Coulombs
coef=0.023e-26;% e^2 / 4 pi epsilon_0 unit of C^2/(F/m) = Joule-metres
hbar=1.05e-34; % reduced Planck's constant in Joule-seconds
m=9.1e-31;     % electron mass in kg
hartree=27.2113834; % conversion factor for atomic units in energy
a0=0.5291772083; %conversion factor for atomic units in length
kb=3.166815422543e-6; %Boltzmann's constant in atomic units of energy per K

% Computational parameters:
NPoints_x=100;
NPoints_E=1000;
NPotentialTerms=30;

% Initialise arrays:
PotentialV=zeros(NPoints_x,1);
Positions=zeros(NPoints_x,1);
Energies=zeros(NPoints_E,1);
xl=zeros(NPoints_E,1);
xr=zeros(NPoints_E,1);
Dh=zeros(NPoints_E,1);
Nh=zeros(NPoints_E,1);
JJ=zeros(NPoints_E,1);

% Computed elements for the numerical analysis:
PhiH=PhiH/hartree;          % hot electrode workfunction in atomic units
PhiC=PhiC/hartree;          % cold electrode workfunction in atomic units
Vbias=Vbias/hartree;
d=d*1e10/a0;           %electrode separation converted to atomic units
E=E/hartree;  % Convert energy range into atomic units
m=1;
e=1;
hbar=1;
coef=1;

dx=d/(NPoints_x-1);   % Interval between points in the discretised potential...
dE=(E(2)-E(1))/(NPoints_E-1); % Interval between calculated energies


% Let's do the first thing... the sum of the terms in the "hot" potential...
for CurrentPointIndex=1:NPoints_x
    
    CurrentValue_x=dx*(CurrentPointIndex-1);
    
    Positions(CurrentPointIndex)=CurrentValue_x;
    
    % First term in the potential formula:
    PotentialV(CurrentPointIndex)=PhiH-...
        (e*Vbias+PhiH-PhiC)*CurrentValue_x/d;
    
    % First term inside the brackets...
    sum=0.25/CurrentValue_x;
    
    % Now for the sum over the terms in the series...
    for CurrentPotentialTerm=1:NPotentialTerms
        sum=sum+0.5*...
            (...
            CurrentPotentialTerm*d/((CurrentPotentialTerm*d)^2-CurrentValue_x^2)...
            -1/(CurrentPotentialTerm*d)...
            );
    end
    PotentialV(CurrentPointIndex)=...
        PotentialV(CurrentPointIndex)-sum*coef;
    
end

disp(['        Vacuum width is = ',num2str(d*a0),' Angstrom']);
disp(['    Work function (hot) = ',num2str(PhiH*hartree),' eV']);
disp(['   Work function (cold) = ',num2str(PhiC*hartree),' eV']);
disp(['                   Bias = ',num2str(Vbias*hartree),' V']);
disp(['Potential maximum value = ',num2str(max(PotentialV)*hartree),' eV']);


% Now for the transmission probability, D_h, as a function of incident
% electron energy...

for CurrentEnergyIndex=1:NPoints_E
    
    Energies(CurrentEnergyIndex)=E(1)+(CurrentEnergyIndex-1)*dE;
    
    % If the energy is less than the maximum, then find the roots...
    if Energies(CurrentEnergyIndex)<max(PotentialV(:))
        for CurrentPointIndex=1:NPoints_x-1
            % Find the values of "x" where the energy goes over the current
            % energy...
            if (PotentialV(CurrentPointIndex)<Energies(CurrentEnergyIndex)) ...
                    &&...
                    (PotentialV(CurrentPointIndex+1)>Energies(CurrentEnergyIndex))
                % Linearly interpolate between these points...
                xl(CurrentEnergyIndex)=Positions(CurrentPointIndex)+...
                    dx*...
                    (Energies(CurrentEnergyIndex)-PotentialV(CurrentPointIndex))/...
                    (PotentialV(CurrentPointIndex+1)-PotentialV(CurrentPointIndex));
                
                % Work out the "trapezium rule" contribution from this point to
                % the right-hand limit...
                
                trapezium=0.5*(Positions(CurrentPointIndex)-xl(CurrentEnergyIndex))...
                    *(PotentialV(CurrentPointIndex+1)-Energies(CurrentEnergyIndex));
                
                
                
                % Find the values of "x" where the energy goes back under the current
                % energy...
            elseif (PotentialV(CurrentPointIndex)>Energies(CurrentEnergyIndex)) ...
                    &&...
                    (PotentialV(CurrentPointIndex+1)<Energies(CurrentEnergyIndex))
                % Linearly interpolate between these points too...
                xr(CurrentEnergyIndex)=Positions(CurrentPointIndex)+...
                    dx*...
                    (Energies(CurrentEnergyIndex)-PotentialV(CurrentPointIndex))/...
                    (PotentialV(CurrentPointIndex+1)-PotentialV(CurrentPointIndex));
                
                trapezium=trapezium+...
                    0.5*(xr(CurrentEnergyIndex)-Positions(CurrentPointIndex))...
                    *(PotentialV(CurrentPointIndex)-Energies(CurrentEnergyIndex));
                
            elseif (PotentialV(CurrentPointIndex)>Energies(CurrentEnergyIndex))
                trapezium=trapezium+...
                    dx*(...
                    0.5*(PotentialV(CurrentPointIndex)+PotentialV(CurrentPointIndex+1))...
                    -Energies(CurrentEnergyIndex)...
                    );
                Dh(CurrentEnergyIndex)=exp(-2*sqrt(2*m)/hbar*trapezium);
            end
        end
        
    else
        % If the energy of the electron is greater than the potential
        % barrier, then the probability of transmission is one!
        Dh(CurrentEnergyIndex)=1;
    end
    
end

tempnan=isnan(Dh);
tempinf=isinf(Dh);
for i=1:NPoints_E
    if tempnan(i)==1 || tempinf(i)==1
        Dh(i)=0;
    end
end

%Let's plot some graphs to help us understand what's been calculated...

% Let's select the energy range to plot based upon where the probability of
% transmissio becomes significant, and then symmetrically about the
% potential maximum...
for i=1:NPoints_E
    if Dh(i)>1e-5
        emin=Energies(i);
        break;
    end
end

emax=2*max(PotentialV)-emin;

root=zeros(NPoints_x,1);
for i=1:NPoints_x
    root(i)=root(i)-0.3;
end

% Now for a plot of the potential barrier:
figure(1);
plot(Positions(:),PotentialV(:),'-o');
hold on
plot(Positions(:),root(:),'r--','LineWidth',1.2)
axis([0,d,emin,0.1]);
%title('Calculated Potential');
xlabel('Position, x (au)');
XTicks([0 1 2 3 4 5 6 7 8 9 9.4486])
%xticklabels({'0','1','2','3','4','5','6','7','8','9','d'})
ylabel('V_h(x) (au)');
grid on;
legend('V_h(x_h)','E=-0.3au')

% and now for the transmission probability...
figure(2);
plot(Dh(:),Energies(:),'-o');
axis([0,1,emin,emax]);
title('Calculated transmission probability');
xlabel('D_h(E)');
ylabel('Electron energy (au)');
grid on;

% So far, so good.  We now need to work out the number density of the
% electrons as a function of the energy on this scale.

kT=kb*TH;
constant=1.1*kT/2/pi^2;

for CurrentEnergyIndex=1:NPoints_E
    
    % Computers struggle with very big and very small numbers!
    if Energies(CurrentEnergyIndex)/kT > -500
    Nh(CurrentEnergyIndex)=2*constant*...
        log(1+exp(-Energies(CurrentEnergyIndex)/kT));
    else
        Nh(CurrentEnergyIndex)=-2*constant*Energies(CurrentEnergyIndex)/kT;
    end
end
%disp(Nh(NPoints_E))
tempnan=isnan(Nh);
tempinf=isinf(Nh);
for i=1:NPoints_E
    if tempnan(i)==1 || tempinf(i)==1
        Nh(i)=0;
    end
end

% and now for a plot of the number density of electrons...
figure(3);
plot(Nh(:),Energies(:),'-o');
axis([0,inf,emin,emax]);
title('Calculated electron number density');
xlabel('N_h(E)');
ylabel('Electron energy (au)');
grid on;

% Finally, to the current...

% The integrand can be generated readily within matlab:
Jtotal=0;
for CurrentEnergyIndex=1:NPoints_E
    
    JJ(CurrentEnergyIndex)=Nh(CurrentEnergyIndex)*Dh(CurrentEnergyIndex);
    
    if CurrentEnergyIndex>1
        
        Jtotal=Jtotal+...
            0.5*(JJ(CurrentEnergyIndex)+JJ(CurrentEnergyIndex-1))*dE;
        
    end
    
    
end

disp(['Total current density = ',num2str(Jtotal),' atomic units']);

% and now for a plot of the number density of electrons...
figure(4);
plot(JJ(:),Energies(:),'-o');
axis([0,inf,emin,emax]);
title('Calculated integrand for current density integral');
xlabel('N_h(E) x D_h(E)');
ylabel('Electron energy (au)');
grid on;

