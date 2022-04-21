clear all vars
clc
close all force

%% Sanity check: pressure at focal plane (0); source: spherical focused) transducer (T)
% considering only one strip (x axis)
% need to compare to jinc function
%No phase plate in this code. We are looking for pressure distribution of a focused, circular piston

W_bar = waitbar(0,'Please wait...');

%% Part (a) define a geometry and mesh the geometry to perform Rayleigh intergral...
% Define the plane from the transducer to evaluate the rayleigh intergral.

F = 100e-3;         %Focal length
f = 1.092e6;        %Frequency of Ultrasonic wave
rho0 = 1000;        %Density of water
c0 = 1482;          %sound speed in water
K = 2*pi*f/c0;      %wavenumber
lambda = 2*pi/K;    %wavelength
h = lambda/3;       %square differential element edge length as reported in Terzi et al.
a = 50e-3;          %Radius of the circular transducer
u0 = 1;             %Particle velocity on the surface of the transducer in m/s assume uniform across entire surface of tranducer
size_T=ceil(a/h);          %number of elements to use along edge of mesh
P_T=zeros(size_T,size_T); %Initialize the coordinate mesh on the surface of the transducer->Then fill in the position of each mesh grid
Count=0;
for i=1:length(P_T(:,1))%Index Y
    tic %start the timer
    for j=1:length(P_T(1,:))%Index X
        if ((i-(length(P_T(:,1))+1)/2)^2+(j-(length(P_T(:,1))+1)/2)^2)<=(length(P_T(:,1))/2)^2% Draw a circle on the square matrix
            P_T(i,j)=1;
            Count=Count+1;
        end
    end
    time=toc %stop the timer
    waitbar_show(i,P_T(:,1),time,W_bar)  
end
%imshow(P_T)     %show an image of the transducer
S_T = P_T*pi*a^2/Count; %Assign each mesh square with a differential area dS_T

%introduce focal plane
%D_0 = 6e-3;        %Radius of the focal plane--small window because the beam is quite focused 
D_0 = 10e-3;        %a little wider
z0 = 100e-3;  
size_0 = 111;      %Number of partition points along the entire line
P_0 = zeros(size_0,1); %Initialize the pressure pattern on focal line (2D focal plane --> 1D focal line)

for i = 1:length(P_0(:,1))%Index Y on focal plane
    tic
for j = 1:length(P_0(1,:))%Index X on focal plane
        Sum = 0;
                for k = 1:length(P_T(:,1))
                for l = 1:length(P_T(1,:))
                    if P_T(k,l) ~= 0
                        Y_T = (length(P_T(:,1))/2)-k;%Positon in number of squares not converted to actual distance
                        X_T = -1*((length(P_T(1,:))/2)-l);%Positon in number of squares not converted to actual distance
                        
                        Y_0 = (length(P_0(:,1))/2)-i;%Positon in number of squares not converted to actual distance
                        X_0 = -1*((length(P_0(1,:))/2)-j);%Positon in number of squares not converted to actual distance
          
                        Y_T = Y_T*a/(length(P_T(:,1))/2);%Covert to Real Y coordinate of the point
                        X_T = X_T*a/(length(P_T(:,1))/2);%Covert to Real X coordinate of the point
                        Z_T = 0.5*(2*F - sqrt(4*F^2 - 4*(X_T^2+Y_T^2))); %see supplemental notes for derivation of this geometry

                        Y_0 = Y_0*D_0/(length(P_0(:,1))/2);%Covert to Real Y coordinate of the point
                        X_0 = X_0*D_0/(length(P_0(:,1))/2);%Covert to Real X coordinate of the point
                        R_0 = sqrt((X_0-X_T)^2+(Y_0-Y_T)^2+(z0-Z_T)^2);%Distance between the phase plate and the focal plane

                        Sum = Sum-sqrt(-1)*rho0*K*c0*u0*S_T(k,l)*exp(sqrt(-1)*K*R_0)/(2*pi*R_0);%Perform Rayleigh integral on each mesh point on the transducer
                
                    end
                end
                end
            P_0(i,j)=Sum;%Assign rayleigh intergral result onto each point on the phase plate
            Count=Count+1;%Record how many partitions there are for the phase plate
end    
        time=toc
        waitbar_show(i,P_0(:,1),time,W_bar)
end


set(groot,'DefaultAxesFontSize', 20)
set(groot,'DefaultLineLineWidth', 2)
set(groot,'DefaultAxesLineWidth', 2)

x = linspace(-10e-3,10e-3,length(P_0)); %x-axis in meters; note that this equals sigma since there is no y component
G = K*a^2/(2*F); %focusing gain

%jinc = rho0*c0*u0*G*(2*besselj(1,K*a*x/F)./(K*a*x/F)); %magnitude from Acoustics II lecture, April 11 2022

jinc = -sqrt(-1)*K*rho0*c0 *exp((sqrt(-1))*K*F)/F .* exp(sqrt(-1)*K*x.^2/(2*F)).*(a^2/2*u0).*(2*besselj(1,K*a*x/F)./(K*a*x/F)); %complex pressure solution

subplot(1,2,1)
plot(x*1e3,abs(jinc/(rho0*c0*u0*G)),'color','b') %plot in mm
hold on
%plot(x*1e3,abs(P_0)/abs(max(P_0)),'+','color','k') %normalized by max of P_0 (Terzi et al's scale)
plot(x*1e3,abs(P_0)/(rho0*c0*u0*G),'*','color','r')
title('magnitude','Interpreter','latex')
l = legend('analytical','numerical','Interpreter','latex');
legend('Location','northeast')
xlabel('$x$ [mm]','interpreter','latex');
ylabel('$\rho c_0 u_0 ka^2/2d$','interpreter','latex');  
axis square

subplot(1,2,2)
plot(x*1e3,(angle(jinc/(rho0*c0*u0*G))),'color','b') %plot in mm
hold on
plot(x*1e3,(angle(P_0/(rho0*c0*u0*G))),'*','color','r')
title('phase','Interpreter','latex')
%l = legend('analytical','numerical','Interpreter','latex');
%legend('Location','northeast')
xlabel('$x$ [mm]','interpreter','latex');
ylabel('phase [rad]','interpreter','latex');  
axis square

