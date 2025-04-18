syms t kk

global N M NM f MM y0 Numexp Totexp
global loss_values
loss_values = [];  % the loss values array

%% C-rate (Put your C-rate)
Crate = -0.05;

%% Experimental data (Users can input their experimental conditions)
Numexp=68506;                           % this is the number of experimental data points
Totexp=68506;                           % experimental discharging time (or charging time) - for our data, this matches the # data points 
x_change = 1.00;                        % this parameter is used to set the PSO search space, i.e the 'p-value'. 

%% Number of parameters to be identified (Input the number of parameters to be identified)
n_vars = 8;                            % Number of parameters to be identified
kk = sym('kk',[1 n_vars]);             % Arrays for parameters to be identified

%% Node number (Change # of node- N: Cathode, M: Membrance, NM: Cathode) (Put your number of node points)  
N=2;
M=2;
NM=2;

%% Design Parameters
ep=0.335;                 % Porosity at positive 
es=0.47;                  % Porosity at membrane 
en=0.25;                  % Porosity at negative
brugp=2.43;               % Bruggeman coefficient at positive
brugs=2.57;               % Bruggeman coefficient at separator
brugn=2.91;               % Bruggeman coefficient at negative
lp=75.6e-6;               % Thickness at positive (unit:m)
ls=12e-6;                 % Thickness at membrance (unit:m)
ln1=85.2e-6;              % Thickness at negative (unit:m)
Rpp=5.22e-6;              % Radius of solid particle at positive (unit:m)
Rpn=5.86e-6;              % Radius of solid particle at negative (unit:m)
F=96487;                  % Faraday constant (unit: C/mol)
R=8.3143;                 % Ideal gas constant (unit: J/(mol K))
t1=0.363;                 % Transference coefficient                
ap=(3/Rpp)*(1-ep);        % Particle surface area to volume at positive (unit: m^2/m^3)
an=(3/Rpn)*(1-en);        % Particle surface area to volume at negative (unit: m^2/m^3)
T=298.15;                 % Temperature (K)
Acell=0.11;               % Electrode area (m^2)
Capa=5;                   % Nominal capacity (Ah)
iapp = Capa*Crate/Acell;

%% Transport parameters          
c0=1000;                     % Electrolyte concentration (unit:mol/m3) %csol
D1=kk(1)*10^(-9);            % Electrolyte diffusion coefficient (unit:m^2/s)
Kappa=kk(2);                 % Conductivity (unit:S/m)
ctp=51765;                   % Maximum solid phase concentration at positive (unit:mol/m^3)
ctn=29583;                   % Maximum solid phase concentration at negative (unit:mol/m^3)
Dbulk=D1;                    % Electrolyte diffusivity (unit:m^2/s)
sigmap=kk(3);                % Solid phase conductivity at positive (unit:S/m)
sigman=kk(4);               % Solid phase conductivity at negative (unit:S/m)
Dsp=kk(5)*10^(-15);          % Solid particle diffusivity at positive (unit:m^2/s)
Dsn=kk(6)*10^(-14);          % Solid particle diffusivity at negative (unit:m^2/s)

Keffp=Kappa*(ep^brugp);      % Liquid phase conductivity at positive (unit:S/m)
Keffs=Kappa*(es^brugs);      % Liquid phase conductivity at membrane (unit:S/m)
Keffn=Kappa*(en^brugn);      % Liquid phase conductivity at negative (unit:S/m)
D2pos=(ep^brugp)*Dbulk;      % Electrolyte diffusivity at positive (unit:m^2/s)
D2sep=(es^brugs)*Dbulk;      % Electrolyte diffusivity at membrane (unit:m^2/s)
D2neg=(en^brugn)*Dbulk;      % Electrolyte diffusivity at negative (unit:m^2/s)

%% Kinetic parameters
kp=kk(7)*10^(-11);           % Reaction rate constant at positive (unit:m^2.5/(mol^0.5 s))
kn=kk(8)*10^(-12);           % Reaction rate constant at negative (unit:m^2.5/(mol^0.5 s))

%% Step size
h=lp/(N+1);                  % Step size at positive 
h2=ls/(M+1);                 % Step size at separator
h3=ln1/(NM+1);               % Step size at negative

%% Variable declaration
Nt=1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2+1+N+1+M+1+NM+1;
X=cell(Nt,1);

for ii=1:Nt
    X{ii}=symfun(str2sym(sprintf('X_%d(t)',ii)),t);
end

%% Variables 
varsX = [X{1:Nt}];

%% Positive 
% u1: electrolyte concentration
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
u1=sym(zeros(1,1+N+1+M+1+NM+1));
parfor i=1:1+N+1+M+1+NM+1
    u1(i)=X{i};
end

% u2: Surface concentration
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
u2=sym(zeros(1,N+NM));
parfor i=1:N
    u2(i+1)=X{i+1+N+1+M+1+NM+1};
end

for i=1:NM
    u2(i+1+N+1+M+1)=X{i+1+N+1+M+1+NM+1+N};
end

% u3: Average concentration
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
u3=sym(zeros(1,N+NM));
parfor i=1:N
    u3(i+1)=X{i+1+N+1+M+1+NM+1+N+NM};
end

for i=1:NM
    u3(i+1+N+1+M+1)=X{i+1+N+1+M+1+NM+1+N+NM+N};
end

% u4: Solid phase potential
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
u4=sym(zeros(1,N+2+NM+2));
parfor i=1:N+2
    u4(i)=X{i+1+N+1+M+1+NM+1+N+NM+N+NM};
end

for i=1:NM+2
    u4(i+1+N+1+M)=X{i+1+N+1+M+1+NM+1+N+NM+N+NM+N+2};
end

% u5: Liquid potential (V)
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
u5=sym(zeros(1,1+N+1+M+1+NM+1));
parfor i=1:1+N+1+M+1+NM+1
    u5(i)=X{i+1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2};
end

%% Mole flux per area per second
% Positive (jp)
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 

jp=sym(zeros(1,N+1));
parfor i=2:N+1
    
theta=u2(i)*ctp/ctp;
Up=-0.8090*theta+4.4875-0.0428*tanh(18.5138*(theta-0.5542))-17.7326*tanh(15.7890*(theta-0.3117))+17.5842*tanh(15.9308*(theta-0.3120));

jp(i)=2*kp*((u1(i)*c0)^0.5)*((ctp-u2(i)*ctp)^0.5)*((u2(i)*ctp)^0.5)*sinh(0.5*F/R/T*(u4(i)-u5(i)-Up));

end

%% Negative (jn)
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 

jn=sym(zeros(1,1+N+1+M+1+NM));
parfor i=1+N+1+M+1+1:1+N+1+M+1+NM

theta=u2(i)*ctn/ctn;
Un=1.9793*exp(-39.3631*theta)+0.2482-0.0909*tanh(29.8538*(theta-0.1234))-0.04478*tanh(14.9159*(theta-0.2769))-0.0205*tanh(30.4444*(theta-0.6103));

jn(i)=2*kn*((u1(i)*c0)^0.5)*((ctn-u2(i)*ctn)^0.5)*((u2(i)*ctn)^0.5)*sinh(0.5*F/R/T*(u4(i)-u5(i)-Un));

end

% u1: Electrolyte concentration (mol/m3)

% Positive Electrode
dudxf1=1/2/h*(-u1(3)-3*u1(1)+4*u1(2));
dudxb1=1/2/h*(u1(N)+3*u1(N+2)-4*u1(N+1));
dudxf1_2=1/2/h2*(-u1(N+4)-3*u1(N+2)+4*u1(N+3));

bc11=dudxf1;
bc21=D2pos*dudxb1-D2sep*dudxf1_2;

eq1=sym(zeros(1,1+N+1+M+1+NM+1));

eq1(1)=0==bc11;                                          % Boundary equation

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=2:N+1                                           % Internal equation
    
d2udx21=(1/(h^2))*(u1(i-1)-2*u1(i)+u1(i+1));

eq1(i)=diff(u1(i))==(D2pos*d2udx21+ap*(1-t1)*jp(i)/c0)/ep;     

end

eq1(N+2)=0==bc21;                                        % Boundary equation

% Separator
dudxb1_2=1/2/h2*(u1(N+2+M-1)+3*u1(N+2+M+1)-4*u1(N+2+M));
dudxf1_3=1/2/h3*(-u1(N+2+M+1+1+1)-3*u1(N+2+M+1)+4*u1(N+2+M+1+1));

bc31=D2sep*dudxb1_2-D2neg*dudxf1_3;

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=N+2+1:N+2+M                                     % Internal equation

d2udx21=(1/(h2^2))*(u1(i-1)-2*u1(i)+u1(i+1));
eq1(i)=diff(u1(i))==(D2sep*d2udx21)/es;

end

eq1(N+2+M+1)=0==bc31;                                    % Boundary equation
      
% Negative Electrode
dudxb1_3=1/2/h3*(u1(1+N+1+M+1+NM-1)+3*u1(1+N+1+M+1+NM+1)-4*u1(1+N+1+M+1+NM));

bc41=dudxb1_3;

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=N+2+M+1+1:N+2+M+1+NM                            % Internal equation

d2udx21=(1/(h3^2))*(u1(i-1)-2*u1(i)+u1(i+1));
eq1(i)=diff(u1(i))==(D2neg*d2udx21+an*(1-t1)*jn(i)/c0)/en;

end

eq1(1+N+1+M+1+NM+1)=0==bc41;                             % Boundary equation

% u2: Surface concentration
 
% Positive electrode

eq2=sym(zeros(1,N+2+M+1+NM ));

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=2:N+1                                           % Internal equation
    
eq2(i)=0==-u2(i)+u3(i)-jp(i)*Rpp/Dsp/5/ctp;

end

% Negative electrode

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=N+2+M+1+1:N+2+M+1+NM                            % Internal equation
    
eq2(i)=0==-u2(i)+u3(i)-jn(i)*Rpn/Dsn/5/ctn;

end

% u3: Average concentration

% Positive
eq3=sym(zeros(1,1+N+1+M+1+NM));

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=2:N+1                                           % Internal equation
    
eq3(i)=diff(u3(i))==-3*jp(i)/Rpp/ctp;

end

% Negative electrode
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=1+N+1+M+1+1:1+N+1+M+1+NM                        % Internal equation
        
eq3(i)=diff(u3(i))==-3*jn(i)/Rpn/ctn;

end

% u4: Solid Potential 

% positive
dudxf4=1/2/h*(-u4(3)-3*u4(1)+4*u4(2));
dudxb4=1/2/h*(u4(N)+3*u4(N+2)-4*u4(N+1));

bc14=dudxf4+iapp/sigmap; 
bc24=dudxb4;

eq4=sym(zeros(1,1+N+1+M+1+NM+1));
 
eq4(1)=0==bc14;                                          % Boundary equation

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=2:N+1                                           % Internal equation
    
d2udx24=(1/(h^2))*(u4(i-1)-2*u4(i)+u4(i+1));

eq4(i)=0==d2udx24-ap*F*jp(i)/sigmap;

end

eq4(N+2)=0==bc24;  

% Negative 
dudxf4_3=1/2/h3*(-u4(N+2+M+1+1+1)-3*u4(N+2+M+1)+4*u4(N+2+M+1+1));
dudxb4_3=1/2/h3*(u4(1+N+1+M+1+NM-1)+3*u4(1+N+1+M+1+NM+1)-4*u4(1+N+1+M+1+NM));

bc34=dudxf4_3;
bc44=dudxb4_3+iapp/sigman;

eq4(N+2+M+1)=0==bc34;                                     % Boundary equation  

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=N+2+M+1+1:N+2+M+1+NM                             % Internal equation
    
d2udx24=(1/(h3^2))*(u4(i-1)-2*u4(i)+u4(i+1));
eq4(i)= 0==d2udx24-an*F*jn(i)/sigman;

end

eq4(1+N+1+M+1+NM+1)=0==bc44;                              % Boundary equation
 
% u5: Liquid phase potential (V)
 
% Positive Electrode
dudxf5=1/2/h*(-u5(3)-3*u5(1)+4*u5(2));
dudxb5=1/2/h*(u5(N)+3*u5(N+2)-4*u5(N+1));
dudxf5_2=1/2/h2*(-u5(N+2+2)-3*u5(N+2)+4*u5(N+2+1));

bc15=dudxf5;
bc25=Keffp*dudxb5-Keffs*dudxf5_2; 

eq5=sym(zeros(1,1+N+1+M+1+NM+1));
   
eq5(1)=0==bc15;                                           % Boundary equation

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=2:N+1                                            % Internal equation

dudx1=1/2/h*(u1(i+1)-u1(i-1));
dudx4=1/2/h*(u4(i+1)-u4(i-1));
dudx5=1/2/h*(u5(i+1)-u5(i-1));

eq5(i)=0==-sigmap*dudx4-Keffp*dudx5+(2*Keffp*R*T*(1-t1)*dudx1)/F/u1(i)-iapp;

end

eq5(N+2)=0==bc25;                                         % Boundary equation

% Separator
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=N+2+1:N+2+M                                      % Internal equation
    
dudx1=1/2/h2*(u1(i+1)-u1(i-1));
dudx5=1/2/h2*(u5(i+1)-u5(i-1));

eq5(i)=0==-Keffs*dudx5+(2*Keffs*R*T*(1-t1)*dudx1)/F/u1(i)-iapp;

end

dudxb5_2=1/2/h2*(u5(1+N+1+M-1)+3*u5(1+N+1+M+1)-4*u5(1+N+1+M));
dudxf5_3=1/2/h3*(-u5(1+N+1+M+1+1+1)-3*u5(1+N+1+M+1)+4*u5(1+N+1+M+1+1));

bc35=Keffs*dudxb5_2-Keffn*dudxf5_3;

eq5(1+N+1+M+1)=0==bc35;                                   % Boundary equation

% Negative
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor i=1+N+1+M+1+1:1+N+1+M+1+NM                         % Internal equation

dudx1=1/2/h3*(u1(i+1)-u1(i-1));
dudx4=1/2/h3*(u4(i+1)-u4(i-1));
dudx5=1/2/h3*(u5(i+1)-u5(i-1));

eq5(i)=0==-sigman*dudx4-Keffn*dudx5+(2*Keffn*R*T*(1-t1)*dudx1)/F/u1(i)-iapp;

end

bc45=u5(1+N+1+M+1+NM+1);

eq5(1+N+1+M+1+NM+1)=0==bc45;                             % Boundary equation

%% A single-step iteration-free initialization approach (Patented by Venkat Subramanian Group at UT Austin)- Allows for Only Acamedic purpose
mu=10^(-3);
q=1000;
initime=200;
ff=1/2*tanh(q*(t-initime))+1/2;
%ff=1/2*tanh(q*(t))+1/2;

eqn1=sym(zeros(1,1+N+1+M+1+NM+1));

eqn1(1)=-mu*(diff(rhs(eq1(1)))-diff(lhs(eq1(1))))-rhs(eq1(1))+lhs(eq1(1));

% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor ii=2:N+1     
eqn1(ii)=lhs(eq1(ii))-rhs(eq1(ii))*ff;
end
eqn1(N+2)=-mu*(diff(rhs(eq1(N+2)))-diff(lhs(eq1(N+2))))-rhs(eq1(N+2))+lhs(eq1(N+2));
parfor ii=N+2+1:N+2+M     
eqn1(ii)=lhs(eq1(ii))-rhs(eq1(ii))*ff;
end
eqn1(N+2+M+1)=-mu*(diff(rhs(eq1(N+2+M+1)))-diff(lhs(eq1(N+2+M+1))))-rhs(eq1(N+2+M+1))+lhs(eq1(N+2+M+1));
parfor ii=N+2+M+1+1:N+2+M+1+NM     
eqn1(ii)=lhs(eq1(ii))-rhs(eq1(ii))*ff;
end
eqn1(1+N+1+M+1+NM+1)=-mu*(diff(rhs(eq1(1+N+1+M+1+NM+1)))-diff(lhs(eq1(1+N+1+M+1+NM+1))))-rhs(eq1(1+N+1+M+1+NM+1))+lhs(eq1(1+N+1+M+1+NM+1));

eqn2=sym(zeros(1,N+2+M+1+NM));
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor ii=2:N+1     
eqn2(ii)=-mu*(diff(rhs(eq2(ii)))-diff(lhs(eq2(ii))))-rhs(eq2(ii))+lhs(eq2(ii));
end
parfor ii=N+2+M+1+1:N+2+M+1+NM                                 
eqn2(ii)=-mu*(diff(rhs(eq2(ii)))-diff(lhs(eq2(ii))))-rhs(eq2(ii))+lhs(eq2(ii));
end

eqn3=sym(zeros(1,1+N+1+M+1+NM));
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor ii=2:N+1                                                
eqn3(ii)=lhs(eq3(ii))-rhs(eq3(ii))*ff;
end
parfor ii=1+N+1+M+1+1:1+N+1+M+1+NM                                
eqn3(ii)=lhs(eq3(ii))-rhs(eq3(ii))*ff;
end

eqn4=sym(zeros(1,1+N+1+M+1+NM+1));
eqn4(1)=-mu*(diff(rhs(eq4(1)))-diff(lhs(eq4(1))))-rhs(eq4(1))+lhs(eq4(1));
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor ii=2:N+1                                               
eqn4(ii)=-mu*(diff(rhs(eq4(ii)))-diff(lhs(eq4(ii))))-rhs(eq4(ii))+lhs(eq4(ii));
end
eqn4(N+2)=-mu*(diff(rhs(eq4(N+2)))-diff(lhs(eq4(N+2))))-rhs(eq4(N+2))+lhs(eq4(N+2));
eqn4(N+2+M+1)=-mu*(diff(rhs(eq4(N+2+M+1)))-diff(lhs(eq4(N+2+M+1))))-rhs(eq4(N+2+M+1))+lhs(eq4(N+2+M+1));
parfor ii=N+2+M+1+1:N+2+M+1+NM                                 
eqn4(ii)=-mu*(diff(rhs(eq4(ii)))-diff(lhs(eq4(ii))))-rhs(eq4(ii))+lhs(eq4(ii));
end
eqn4(1+N+1+M+1+NM+1)=-mu*(diff(rhs(eq4(1+N+1+M+1+NM+1)))-diff(lhs(eq4(1+N+1+M+1+NM+1))))-rhs(eq4(1+N+1+M+1+NM+1))+lhs(eq4(1+N+1+M+1+NM+1));

eqn5=sym(zeros(1,1+N+1+M+1+NM+1));
eqn5(1)=-mu*(diff(rhs(eq5(1)))-diff(lhs(eq5(1))))-rhs(eq5(1))+lhs(eq5(1));
% 'parfor' can be replaced with 'for' if users do not have the parallel computing toolbox 
parfor ii=2:N+1                                           
eqn5(ii)=-mu*(diff(rhs(eq5(ii)))-diff(lhs(eq5(ii))))-rhs(eq5(ii))+lhs(eq5(ii));
end
eqn5(N+2)=-mu*(diff(rhs(eq5(N+2)))-diff(lhs(eq5(N+2))))-rhs(eq5(N+2))+lhs(eq5(N+2));
parfor ii=N+2+1:N+2+M                                      
eqn5(ii)=-mu*(diff(rhs(eq5(ii)))-diff(lhs(eq5(ii))))-rhs(eq5(ii))+lhs(eq5(ii));
end
eqn5(1+N+1+M+1)=-mu*(diff(rhs(eq5(1+N+1+M+1)))-diff(lhs(eq5(1+N+1+M+1))))-rhs(eq5(1+N+1+M+1))+lhs(eq5(1+N+1+M+1));
parfor ii=1+N+1+M+1+1:1+N+1+M+1+NM                          
eqn5(ii)=-mu*(diff(rhs(eq5(ii)))-diff(lhs(eq5(ii))))-rhs(eq5(ii))+lhs(eq5(ii));
end
eqn5(1+N+1+M+1+NM+1)=-mu*(diff(rhs(eq5(1+N+1+M+1+NM+1)))-diff(lhs(eq5(1+N+1+M+1+NM+1))))-rhs(eq5(1+N+1+M+1+NM+1))+lhs(eq5(1+N+1+M+1+NM+1));

%% Total equations
eqs = [eqn1(1:1+N+1+M+1+NM+1),eqn2(1+1:N+1),eqn2(1+N+2+M+1:NM+N+2+M+1),eqn3(1+1:N+1),eqn3(1+N+2+M+1:NM+N+2+M+1),eqn4(1:N+2),eqn4(1+1+N+1+M:NM+2+1+N+1+M),eqn5(1:1+N+1+M+1+NM+1)];
vars = (varsX);

[MM,f] = massMatrixForm(eqs, vars); % takes all system of eqs and puts into mass matrix. f is the right hand side of the eqns

MM = odeFunction(MM, vars,kk);
f = odeFunction(f, vars,kk);

%% Initial condition guess (Put your initial guess)                                  
U(1:1+N+1+M+1+NM+1)=1;                                                                                        % Electrolyte concentration                                            
U(1+1+N+1+M+1+NM+1:N+1+N+1+M+1+NM+1)=.27;                                                                     % Surface con. at positive                                             
U(1+1+N+1+M+1+NM+1+N:NM+1+N+1+M+1+NM+1+N)=0.9014;                                                             % Surface con. at negative                                              
U(1+1+N+1+M+1+NM+1+N+NM:N+1+N+1+M+1+NM+1+N+NM)=.27;                                                           % Average con. at positive                                       
U(1+1+N+1+M+1+NM+1+N+NM+N:NM+1+N+1+M+1+NM+1+N+NM+N)=0.9014;                                                   % Average con. at negative                                           
U(1+1+N+1+M+1+NM+1+N+NM+N+NM:N+2+1+N+1+M+1+NM+1+N+NM+N+NM)=4.30430037;                                        % Solid phase potential at positive                                      
U(1+1+N+1+M+1+NM+1+N+NM+N+NM+N+2:NM+2+1+N+1+M+1+NM+1+N+NM+N+NM+N+2)=.9202000152e-1;                           % Solid phase potential at negative                                         
U(1+1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2:1+N+1+M+1+NM+1+1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2)=0;                    % Liquid potential

y0 = (U(1:1+N+1+M+1+NM+1+1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2));

x=importdata('voltage/Diag10_V4.txt'); % change this
DiagX_Y = 'Diag10_V4';  % again, set this dynamically

sub_folder_path = fullfile(DiagX_Y);
if ~exist(sub_folder_path, 'dir')
    mkdir(sub_folder_path);
end

%% PSO identificaiton (Input your upper and lower bounds for parameters to be identified. Input your options for PSO)
pp=x_change;               % Deviation for upper and lower bounds (percentage) - original value is 0.3
D10=1;                    % Electrolyte diffusion coefficient  
Kappa0=1.17;              % Conductivity
sigmap0=0.18;             % Solid phase conductivity at positive
sigman0=215;              % Solid phase conductivity at negative
Dsp0=4;                   % Solid particle diffusivity at positive (10-15) 1.48
Dsn0=3.3;                 % Solid particle diffusivity at negative (10-14) 1.74*10-15 or 0.174 * 10-14
kp0=0.7;                  % Reaction rate constant at positive 21.2
kn0=0.7;                  % Reaction rate constant at negative 1112

%Lower_bound = [D10-pp*D10 Kappa0-pp*Kappa0 sigmap0-pp*sigmap0 sigman0-pp*sigman0 Dsp0-pp*Dsp0 Dsn0-pp*Dsn0 kp0-pp*kp0 kn0-pp*kn0];
Lower_bound = [0 0 0 0 0 0 0 0]; % lower boundary when pp => 1
Upper_bound = [D10+pp*D10 Kappa0+pp*Kappa0 sigmap0+pp*sigmap0 sigman0+pp*sigman0 Dsp0+pp*Dsp0 Dsn0+pp*Dsn0 kp0+pp*kp0 kn0+pp*kn0]; 

options = optimoptions('particleswarm', 'SwarmSize', 10, 'Display', 'iter', 'FunctionTolerance', 0.01);

tic
[kk,fval] = particleswarm(@(kk)P2Dobj(kk), n_vars, Lower_bound, Upper_bound, options);
toc

% some code to save data as part of experiments
%loss_filename = sprintf('loss_curve_data_pp%0.2f.mat', x_change); 
%full_loss_path = fullfile(sub_folder_path, loss_filename);
%save(full_loss_path, 'loss_values');

% the loss data does not produce a smooth curve due to how the PSO
% works, so visualising it is not very useful. instead, the final value is useful
% as it represents the loss for the solution that PSO converged on, instead
% of the intermediary guesses for paramteres and the corresponding loss values,
% which are often kind of bad.

D1=kk(1)*10^(-9);          % Electrolyte diffusion coefficient (unit:m2/s)
Kappa=kk(2);               % Conductivity (unit:S/m)
sigmap=kk(3);              % Solid phase conductivity at positive (unit:S/m)
sigman=kk(4);              % Solid phase conductivity at negative (unit:S/m)
Dsp=kk(5)*10^(-15);        % Solid particle diffusivity at positive (unit:m^2/s)
Dsn=kk(6)*10^(-14);        % Solid particle diffusivity at negative (unit:m^2/s)
kp=kk(7)*10^(-11);         % Reaction rate constant at positive (unit:m^2.5/(mol^0.5 s))
kn=kk(8)*10^(-12);         % Reaction rate constant at negative (unit:m^2.5/(mol^0.5 s))

save('identified_params_kinetic','fval','D1','Kappa','sigmap','sigman','Dsp','Dsn','kp','kn');

param_filename = sprintf('identified_params_kinetic_pp%0.2f.txt', x_change); 
% ^ metadata for analysis of how these parameters change over degradation, based on PSO
full_param_path = fullfile(sub_folder_path, param_filename);
fileID = fopen(full_param_path, 'w');

% since it's saved as text, some string manipulation will be rquired via
% python to analyse this data.
fprintf(fileID, 'Optimized Parameters:\n');
fprintf(fileID, 'fval (error): %.4f\n', fval);
fprintf(fileID, 'D1 (Electrolyte diffusion coefficient): %.4e m2/s\n', D1);
fprintf(fileID, 'Kappa (Conductivity): %.4f S/m\n', Kappa);
fprintf(fileID, 'sigmap (Solid phase conductivity at positive): %.4f S/m\n', sigmap);
fprintf(fileID, 'sigman (Solid phase conductivity at negative): %.4f S/m\n', sigman);
fprintf(fileID, 'Dsp (Solid particle diffusivity at positive): %.4e m^2/s\n', Dsp);
fprintf(fileID, 'Dsn (Solid particle diffusivity at negative): %.4e m^2/s\n', Dsn);
fprintf(fileID, 'kp (Reaction rate constant at positive): %.4e m^2.5/(mol^0.5 s)\n', kp);
fprintf(fileID, 'kn (Reaction rate constant at negative): %.4e m^2.5/(mol^0.5 s)\n', kn);
fclose(fileID);
%% Solver

% Adapt solver functions for better conditioning    
M0 = MM(initime,y0(:),kk);
vw = 1./max(abs(M0),[],2);
mw = diag(vw);

F = @(t,y) vw.*f(t,y(:),kk);
M1 = @(t,y) mw*MM(t,y(:),kk);

tsp=100000; % simulation end time 
opt = odeset('Mass', M1,'MStateDependence','weak','RelTol',1e-5,'AbsTol',1e-5,'InitialStep',1e-3,'MaxStep',5,'Events',@stopcondition);
warning('off','all');
tic
[T,Y] = ode15s(F,[0 tsp], y0, opt);
toc

% plot output
figure(1)
p1=plot(T-initime,(Y(:,1+N+1+M+1+NM+1+N+NM+N+NM+1)-Y(:,1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2)),'Linewidth',2);
hold on;

% T and Y are time and voltage data arrays respectively
time_data = T - initime; % adjust time data by subtracting initime
voltage_data = Y(:, 1+N+1+M+1+NM+1+N+NM+N+NM+1) - Y(:, 1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2); % extract voltage data
output_data = [time_data, voltage_data];% single matrix

xlsx_filename = sprintf('tol0.01_pp%0.2f.xlsx', x_change); % data
full_xlsx_path = fullfile(sub_folder_path, xlsx_filename);
writematrix(output_data, full_xlsx_path, 'Sheet', 1); % excel file
disp(['Data written to ', xlsx_filename]);

time1=linspace(0,Totexp,Numexp);
p2=plot(time1,x(:,1),'o','MarkerSize',3,'color','red');
ylim([2.5 4.5]);
set(gca,'FontSize',15);
xlabel('Time(seconds)','FontSize',15);
ylabel('Voltage(V)','FontSize',15);
legend('P2D Model','Experiment');
legend('boxoff');
hold off;
save_fig_1 = fullfile(sub_folder_path, sprintf('tol0.01_pp%0.2f.jpeg', x_change));
saveas(figure(1), save_fig_1);

% plot the data
figure(2);  % Create a new figure window
plot(time_data, voltage_data, 'LineWidth', 2);  % Plot with a line width of 2
title('Normal Voltage Data');
xlabel('Time (s)');
ylabel('Voltage (V)');
grid on;  % turn on the grid
hold off;

% plot the loss curve
figure(3);
plot(loss_values, 'LineWidth', 2);
title('Loss Curve During Optimization');
xlabel('Iteration');
ylabel('Loss Value');
grid on;
hold off;

save_fig_2 = fullfile(sub_folder_path, sprintf('loss_curve_pp%0.2f.jpeg', x_change));
%saveas(figure(3), save_fig_2);

filtered_loss_values = loss_values(loss_values ~= 1000);

% plot the filtered array
figure(4);
plot(filtered_loss_values, 'LineWidth', 2);
ylim([0 0.3]);
title('Loss Curve During Optimization (Filtered)');
xlabel('Iteration');
ylabel('Loss Value');
grid on;

matrix_file_name = fullfile(sub_folder_path, sprintf('filtered_loss_values_pp%0.2f.txt', x_change));
%writematrix(filtered_loss_values, matrix_file_name);

close all;

%% Objective function
function obj=P2Dobj(kk)

global y0 f MM N M NM initime Totexp Numexp loss_values

try
    
% adapt solver functions for better conditioning
M0 = MM(initime,y0(:),kk);
vw = 1./max(abs(M0),[],2); 
mw = diag(vw);    

F = @(t,y) vw.*f(t,y(:),kk);
M1 = @(t,y) mw*MM(t,y(:),kk);

time=linspace(0,Totexp+200,Numexp+3);
opt = odeset('Mass', M1,'MStateDependence','weak','RelTol',1e-5,'AbsTol',1e-5,'InitialStep',1e-3,'MaxStep',5);
warning('off','all');
[T,Y] = ode15s(F,time, y0, opt);

x=importdata('voltage/Diag10_V4.txt'); % this is the loaded experimental data
obj=rms(x(:,1)-(Y(4:Numexp+3,1+N+1+M+1+NM+1+N+NM+N+NM+1)-Y(4:Numexp+3,1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2)));
loss_values = [loss_values; obj];

catch
    obj=1000;    
    loss_values = [loss_values; obj];
end

end

%% Stop condition
function [value,isterminal,direction] = stopcondition(~,Y)

global N M NM

value = (Y(1+N+1+M+1+NM+1+N+NM+N+NM+1)-Y(1+N+1+M+1+NM+1+N+NM+N+NM+N+2+NM+2))-2.5;                             % stop condition 4.2 V
isterminal = 1;                                                                                               % stop the integration
direction = 0;                                                                                                % negative direction

end
