load("thermodynamics_P.dat")
Umean=mean(thermodynamics_P(85:end-2,3))/500;
Pmean=mean(thermodynamics_P(85:end-2,2));
rhovec=500./((thermodynamics_P(85:end-2,4)/3.405).^3);
rhomean=mean(rhovec);
Ustd=std(thermodynamics_P(85:end-2,3))/500;
Pstd=std(thermodynamics_P(85:end-2,2));
rhostd=std(rhovec);
figure(1)
plot(thermodynamics_P(1:end-2,1),thermodynamics_P(1:end-2,end))
figure(2)
plot(thermodynamics_P(1:end-2,1),thermodynamics_P(1:end-2,2))
figure(3)
plot(thermodynamics_P(1:end-2,1),thermodynamics_P(1:end-2,3))
%% p=0.15
close all
load("gr_P_015.dat")
load("gr_rho_01_prove.dat")
load("thermodynamics_P_015.dat")
Umean015=mean(thermodynamics_P_015(85:end-2,3))/500;
Pmean015=mean(thermodynamics_P_015(85:end-2,2));
rhovec015=thermodynamics_P_015(85:end-2,4);
rhomean015=mean(rhovec015);
Ustd015=std(thermodynamics_P_015(85:end-2,3))/500;
Pstd015=std(thermodynamics_P_015(85:end-2,2));
rhostd015=std(rhovec015);
figure(1)
plot(thermodynamics_P_015(1:end-2,1),thermodynamics_P_015(1:end-2,4))
figure(2)
plot(thermodynamics_P_015(1:end-2,1),thermodynamics_P_015(1:end-2,2))
figure(3)
plot(thermodynamics_P_015(1:end-2,1),thermodynamics_P_015(1:end-2,3))
gr_P_015(1:4,1)=zeros(4,1);
figure(4)
plot(gr_P_015(:,2),gr_P_015(:,1))
hold on
plot(gr_rho_01_prove(:,2),gr_rho_01_prove(:,1))

load("gr_P_015_norm.dat")
%gr_P_015_norm(1:2,1)=zeros(2,1);
figure(5)
plot(gr_P_015_norm(:,2),gr_P_015_norm(:,1))
hold on
plot(gr_rho_01_prove(:,2),gr_rho_01_prove(:,1))
%% p=1.5

close all
load("gr_P_15.dat")
load("gr_rho_06.dat")
load("thermodynamics_P_15.dat")
Umean15=mean(thermodynamics_P_15(85:end-2,3))/500;
Pmean15=mean(thermodynamics_P_15(85:end-2,2));
rhovec15=thermodynamics_P_15(85:end-2,4);
rhomean15=mean(rhovec15);
Ustd15=std(thermodynamics_P_15(85:end-2,3))/500;
Pstd15=std(thermodynamics_P_15(85:end-2,2));
rhostd15=std(rhovec15);
figure(1)
plot(thermodynamics_P_15(1:end-2,1),thermodynamics_P_15(1:end-2,4))
figure(2)
plot(thermodynamics_P_15(1:end-2,1),thermodynamics_P_15(1:end-2,2))
figure(3)
plot(thermodynamics_P_15(1:end-2,1),thermodynamics_P_15(1:end-2,3))
%gr_P_015(1,1)=zeros(1,1);
gr_P_15(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_15(:,2),gr_P_15(:,1))
hold on
plot(gr_rho_06(:,2),gr_rho_06(:,1))

%% p=3

close all
load("gr_P_30.dat")
load("gr_rho_07.dat")
load("thermodynamics_P_30.dat")
Umean30=mean(thermodynamics_P_30(85:end-2,3))/500;
Pmean30=mean(thermodynamics_P_30(85:end-2,2));
rhovec30=thermodynamics_P_30(85:end-2,4);
rhomean30=mean(rhovec30);
Ustd30=std(thermodynamics_P_30(85:end-2,3))/500;
Pstd30=std(thermodynamics_P_30(85:end-2,2));
rhostd30=std(rhovec30);
figure(1)
plot(thermodynamics_P_30(1:end-2,1),thermodynamics_P_30(1:end-2,4))
figure(2)
plot(thermodynamics_P_30(1:end-2,1),thermodynamics_P_30(1:end-2,2))
figure(3)
plot(thermodynamics_P_30(1:end-2,1),thermodynamics_P_30(1:end-2,3))
%gr_P_015(1,1)=zeros(1,1);
gr_P_30(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_30(:,2),gr_P_30(:,1))
hold on
plot(gr_rho_07(:,2),gr_rho_07(:,1))

%% p=4.5

close all
load("gr_P_45.dat")
load("thermodynamics_P_45.dat")
Umean45=mean(thermodynamics_P_45(85:end-2,3))/500;
Pmean45=mean(thermodynamics_P_45(85:end-2,2));
rhovec45=thermodynamics_P_45(85:end-2,4);
rhomean45=mean(rhovec45);
Ustd45=std(thermodynamics_P_45(85:end-2,3))/500;
Pstd45=std(thermodynamics_P_45(85:end-2,2));
rhostd45=std(rhovec45);
figure(1)
plot(thermodynamics_P_45(1:end-2,1),thermodynamics_P_45(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 4.5$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_45(1:end-2,1),thermodynamics_P_45(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 4.5$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_45(1:end-2,1),thermodynamics_P_45(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 4.5$','Interpreter','latex')
%gr_P_015(1,1)=zeros(1,1);
gr_P_45(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_45(:,2),gr_P_45(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 4.5$'],'Interpreter','latex')

%% p=6

close all
load("gr_P_60.dat")
load("thermodynamics_P_60.dat")
Umean60=mean(thermodynamics_P_60(85:end-2,3))/500;
Pmean60=mean(thermodynamics_P_60(85:end-2,2));
rhovec60=thermodynamics_P_60(85:end-2,4);
rhomean60=mean(rhovec60);
Ustd60=std(thermodynamics_P_60(85:end-2,3))/500;
Pstd60=std(thermodynamics_P_60(85:end-2,2));
rhostd60=std(rhovec60);
figure(1)
plot(thermodynamics_P_60(1:end-2,1),thermodynamics_P_60(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 6$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_60(1:end-2,1),thermodynamics_P_60(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 6$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_60(1:end-2,1),thermodynamics_P_60(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 6$','Interpreter','latex')
gr_P_60(1,1)=zeros(1,1);
%gr_P_60(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_60(:,2),gr_P_60(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 6$'],'Interpreter','latex')

%% p=7.5

close all
load("gr_P_75.dat")
load("thermodynamics_P_75.dat")
Umean75=mean(thermodynamics_P_75(85:end-2,3))/500;
Pmean75=mean(thermodynamics_P_75(85:end-2,2));
rhovec75=thermodynamics_P_75(85:end-2,4);
rhomean75=mean(rhovec75);
Ustd75=std(thermodynamics_P_75(85:end-2,3))/500;
Pstd75=std(thermodynamics_P_75(85:end-2,2));
rhostd75=std(rhovec75);
figure(1)
plot(thermodynamics_P_75(1:end-2,1),thermodynamics_P_75(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 7.5$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_75(1:end-2,1),thermodynamics_P_75(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 7.5$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_75(1:end-2,1),thermodynamics_P_75(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 7.5$','Interpreter','latex')
gr_P_75(1,1)=zeros(1,1);
%gr_P_75(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_75(:,2),gr_P_75(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 7.5$'],'Interpreter','latex')


%% p=9

close all
load("gr_P_90.dat")
load("thermodynamics_P_90.dat")
Umean90=mean(thermodynamics_P_90(85:end-2,3))/500;
Pmean90=mean(thermodynamics_P_90(85:end-2,2));
rhovec90=thermodynamics_P_90(85:end-2,4);
rhomean90=mean(rhovec90);
Ustd90=std(thermodynamics_P_90(85:end-2,3))/500;
Pstd90=std(thermodynamics_P_90(85:end-2,2));
rhostd90=std(rhovec90);
figure(1)
plot(thermodynamics_P_90(1:end-2,1),thermodynamics_P_90(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 9$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_90(1:end-2,1),thermodynamics_P_90(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 9$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_90(1:end-2,1),thermodynamics_P_90(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 9$','Interpreter','latex')
gr_P_90(1,1)=zeros(1,1);
%gr_P_75(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_90(:,2),gr_P_90(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 9$'],'Interpreter','latex')

%% p=10.5

close all
load("gr_P_105.dat")
load("thermodynamics_P_105.dat")
Umean105=mean(thermodynamics_P_105(85:end-2,3))/500;
Pmean105=mean(thermodynamics_P_105(85:end-2,2));
rhovec105=thermodynamics_P_105(85:end-2,4);
rhomean105=mean(rhovec105);
Ustd105=std(thermodynamics_P_105(85:end-2,3))/500;
Pstd105=std(thermodynamics_P_105(85:end-2,2));
rhostd105=std(rhovec105);
figure(1)
plot(thermodynamics_P_105(1:end-2,1),thermodynamics_P_105(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 10.5$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_105(1:end-2,1),thermodynamics_P_105(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 10.5$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_105(1:end-2,1),thermodynamics_P_105(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 10.5$','Interpreter','latex')
gr_P_90(1,1)=zeros(1,1);
%gr_P_75(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_105(:,2),gr_P_105(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 10.5$'],'Interpreter','latex')

%% p=12

close all
load("gr_P_120.dat")
load("thermodynamics_P_120.dat")
Umean120=mean(thermodynamics_P_120(85:end-2,3))/500;
Pmean120=mean(thermodynamics_P_120(85:end-2,2));
rhovec120=thermodynamics_P_120(85:end-2,4);
rhomean120=mean(rhovec120);
Ustd120=std(thermodynamics_P_120(85:end-2,3))/500;
Pstd120=std(thermodynamics_P_120(85:end-2,2));
rhostd120=std(rhovec120);
figure(1)
plot(thermodynamics_P_120(1:end-2,1),thermodynamics_P_120(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 12$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_120(1:end-2,1),thermodynamics_P_120(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 12$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_120(1:end-2,1),thermodynamics_P_120(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 12$','Interpreter','latex')
gr_P_120(1,1)=zeros(1,1);
%gr_P_75(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_120(:,2),gr_P_120(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 12$'],'Interpreter','latex')

%% p=13.5

close all
load("gr_P_135.dat")
load("thermodynamics_P_135.dat")
Umean135=mean(thermodynamics_P_135(85:end-2,3))/500;
Pmean135=mean(thermodynamics_P_135(85:end-2,2));
rhovec135=thermodynamics_P_135(85:end-2,4);
rhomean135=mean(rhovec135);
Ustd135=std(thermodynamics_P_135(85:end-2,3))/500;
Pstd135=std(thermodynamics_P_135(85:end-2,2));
rhostd135=std(rhovec135);
figure(1)
plot(thermodynamics_P_135(1:end-2,1),thermodynamics_P_135(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 13.5$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_135(1:end-2,1),thermodynamics_P_135(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 13.5$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_135(1:end-2,1),thermodynamics_P_135(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 13.5$','Interpreter','latex')
gr_P_135(1,1)=zeros(1,1);
%gr_P_75(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_135(:,2),gr_P_135(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 13.5$'],'Interpreter','latex')

%% p=15

close all
load("gr_P_150.dat")
load("thermodynamics_P_150.dat")
Umean150=mean(thermodynamics_P_150(85:end-2,3))/500;
Pmean150=mean(thermodynamics_P_150(85:end-2,2));
rhovec150=thermodynamics_P_150(85:end-2,4);
rhomean150=mean(rhovec150);
Ustd150=std(thermodynamics_P_150(85:end-2,3))/500;
Pstd150=std(thermodynamics_P_150(85:end-2,2));
rhostd150=std(rhovec150);
figure(1)
plot(thermodynamics_P_150(1:end-2,1),thermodynamics_P_150(1:end-2,4))
ylabel('$\rho^*$', 'Interpreter','latex')
xlabel('it ($*501$)','Interpreter','latex')
title('Density in LJ MC for $P^* = 15$','Interpreter','latex')
figure(2)
plot(thermodynamics_P_150(1:end-2,1),thermodynamics_P_150(1:end-2,2))
ylabel('Energy (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Energy in LJ MC for $P^* = 15$','Interpreter','latex')
figure(3)
plot(thermodynamics_P_150(1:end-2,1),thermodynamics_P_150(1:end-2,3))
ylabel('Virial Pressure (Reduced Units)', 'Interpreter','latex')
xlabel('it ($*10^{3}$)','Interpreter','latex')
title('Virial Pressure in LJ MC for $P^* = 15$','Interpreter','latex')
gr_P_150(1,1)=zeros(1,1);
%gr_P_75(1:3,1)=zeros(3,1);
figure(4)
plot(gr_P_150(:,2),gr_P_150(:,1))
ylabel('g(r) (Reduced Units)', 'Interpreter','latex')
xlabel('r (reduced units)','Interpreter','latex')
title(['g(r) in LJ MC for $P^* = 150$'],'Interpreter','latex')
