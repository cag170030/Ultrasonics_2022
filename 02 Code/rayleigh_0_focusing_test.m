clear all vars
clc
close all force

%% Sanity check: pressure at focal plane (0); source: spherical focused) transducer (T)
%No phase plate in this code. We are looking for pressure distribution of a focused, circular piston
%C.G. is cmoparing to the plots presented in class on Monday.

W_bar = waitbar(0,'Please wait...');

%% Part (a) define a geometry and mesh the geometry to perform Rayleigh intergral...
% Define the plane from the transducer to evaluate the rayleigh intergral.


F = 100e-3;         %Focal length
f = 1.092e6;        %Frequency of Ultrasonic wave
rho0 = 1000;        %Density of water
c0 = 1481;          %sound speed in water
K = 2*pi*f/c0;      %wavenumber
a = 50e-3;          %Radius of the circular transducer
u0 = 1;             %Particle velocity on the surface of the transducer in m/s assume uniform across entire surface of tranducer
%size_T=151;        %Line along the transducer's diameter
size_T=101;          %Line along the transducer's diameter
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
D_0 = 50e-3;        %Radius of the focal plane--change to 5e-3 once accounting for focusing 
z0 = 100e-3;  
%size_0 = 151;      %Number of partition points along the entire line
size_0 = 101;        %use temporarily for quick computation
P_0 = zeros(size_0,size_0); %Initialize the pressure pattern on focal plane

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

                        Sum = Sum-sqrt(-1)*rho0*K*c0*u0*S_T(k,l)*exp(-sqrt(-1)*K*R_0)/(2*pi*R_0);%Perform Rayleigh integral on each mesh point on the transducer
                    end
                end
                end
            P_0(i,j)=Sum;%Assign rayleigh intergral result onto each point on the phase plate
            Count=Count+1;%Record how many partitions there are for the phase plate
end    
        time=toc
        waitbar_show(i,P_0(:,1),time,W_bar)
end

figure(1)
pcolor(angle(P_0))
axis square
figure(2)
pcolor(abs(P_0))
axis square