clear all vars
clc
close all force
%% This script is meant to simulate acoustics wave from a baffled piston incident on a phase plate. To run this code you will need 'waitbar_show.m'
W_bar = waitbar(0,'Please wait...');

%% Part (a) define a geometry and mesh the geometry to perform Rayleigh intergral...
% And also to define the plane from the piston to evaluate the rayleigh
% intergral.
f=2e5;%Frequency of Ultrasonic wave
a=5e-2;%Radius of the circular piston
co=343;%sound speed in air
uo=1;% Particle velocity on the surface of the piston in m/s assume uniform across entire surface of piston
M_size=101;% Meaning a linear line along the piston's diameter is M_size
P_piston=zeros(M_size,M_size);%Initialize the coordinate mesh on the surface of the piston->Then fill in the position of each mesh grid
Count=0;
%P_piston=gpuArray(P_piston);
for i=1:length(P_piston(:,1))%Index Y
    tic
    for j=1:length(P_piston(1,:))%Index X
        if ((i-(length(P_piston(:,1))+1)/2)^2+(j-(length(P_piston(:,1))+1)/2)^2)<=(length(P_piston(:,1))/2)^2% Draw a circle on the square matrix
            P_piston(i,j)=1;
            Count=Count+1;
        end
    end
    time=toc
        waitbar_show(i,P_piston(:,1),time,W_bar)
    
end
imshow(P_piston)% show an imag of the piston
S_piston=P_piston*pi*a^2/Count; %Assigna each mesh square with an area
D_P=10e-2;%Radius of the phase plate
L=20e-2;%Distance between the piston and the phase plate
P_size=101;%Number of partition points along the entire line
Plate_P=zeros(P_size,P_size);%Innitialize the pressure pattern
Count=0;
rho=1.2;%Density in air

for i=1:length(Plate_P(:,1))%Index Y
    tic
    
    for j=1:length(Plate_P(1,:))%Index X
        Sum=0;
        if ((i-(length(Plate_P(:,1))+1)/2)^2+(j-(length(Plate_P(:,1))+1)/2)^2)<=(length(Plate_P(:,1))/2)^2%Draw a circle on the square matrix
            for k=1:length(P_piston(:,1))
                for l=1:length(P_piston(1,:))
                    if P_piston(k,l)~=0
                        Position_pistonY=(length(P_piston(:,1))/2)-k;%Positon in number of squares not converted to actual distance
                        Position_pistonX=-1*((length(P_piston(1,:))/2)-l);%Positon in number of squares not converted to actual distance
                        Position_PlateY=(length(Plate_P(:,1))/2)-i;%Positon in number of squares not converted to actual distance
                        Position_PlateX=-1*((length(Plate_P(1,:))/2)-j);%Positon in number of squares not converted to actual distance
                        Position_pistonY=Position_pistonY*a/(length(P_piston(:,1))/2);%Covert to Real Y coordinate of the point
                        Position_pistonX=Position_pistonX*a/(length(P_piston(:,1))/2);%Covert to Real X coordinate of the point
                        Position_PlateY=Position_PlateY*D_P/(length(Plate_P(:,1))/2);%Covert to Real Y coordinate of the point
                        Position_PlateX=Position_PlateX*D_P/(length(Plate_P(:,1))/2);%Covert to Real X coordinate of the point
                        dist=sqrt((Position_PlateX-Position_pistonX)^2+(Position_PlateY-Position_pistonY)^2+L^2);%Calculate the distance between the phase plate...
                        %and the circular piston
                        
                        Sum=Sum+rho*2*pi*f*uo*S_piston(k,l)*exp(-sqrt(-1)*2*pi*f*dist/co)/(2*pi*dist);%Perform rayleigh integral on each mesh point on the piston
                    end
                end
            end
            Plate_P(i,j)=Sum;%Assign rayleigh intergral result onto each point on the phase plate
            Count=Count+1;%Record how many partitions there are for the phase plate
        end    
    end
        time=toc
        waitbar_show(i,Plate_P(:,1),time,W_bar)
end

figure(1)
pcolor(angle(Plate_P))
figure(2)
pcolor(abs(Plate_P))
%% Part(b) Now that the rayleigh intergral has been performed on the phase plate  
%The next step should be to define the necessary phase on each of the point on the
%phase plate to yield a acoustics vortex beam First better to sanity check
%the current result
