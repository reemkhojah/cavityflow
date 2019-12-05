clear all;close all;fclose all;clc;load('Entry_ypoints.mat');load('Entry_zpoints.mat')
%% Channel variables
um = 0.169;    %mean velocity
Ly = 40e-6;         % domain height
Lz = 150e-6;        % domain width 
ny =400; nz = 1500; % Number of cells
%% Channel flow calculation
a = Lz/2; b = Ly/2; y = -b:Ly/ny:b;z = -a:Lz/nz:a; % Create grid in y direction % Create grid in z direction
alphaa = b/a;    % Aspect ratio of the domain
m = 1.5 + 0.5*alphaa^(-1.4); %  constants
if alphaa <= 1/3
    n = 2;
else
    n = 2 + 0.3*(alphaa - 1/3);
end
for i=1:nz+1
    for j=1:ny+1
        u(i,j) = um * ((m+1)/m) * ((n+1)/n) * (1-(abs(y(j)/b))^n) * (1-(abs(z(i)/a))^m);
    end
end
%% Process and plot streamline entry points
subplot(1,3,1);contourf(y,z,u); c.LineWidth = 0;colorbar; axis equal;xlim([-20e-6 20e-6]); ylim([-75e-6 75e-6]);hold on; % contourf(y,z,u,5,'ShowText','on')
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',10); xlabel( 'Y^* (\mum)', 'fontsize',10); ylabel( 'Z^* (\mum)', 'fontsize',10);box on; grid off;
% ------------------- streamline points from simulation
ypoints = (Entry_ypoints-20)./1000000; zpoints = (Entry_zpoints-75)./1000000; %adjust channel height and width
%  Check the velocity of entry points
subplot(1,3,2);scatter(ypoints,zpoints,'.'); axis equal;xlim([-20e-6 20e-6]); ylim([-75e-6 75e-6]); %check the right place in the channel
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',10); xlabel( 'Y^* (\mum)', 'fontsize',10); ylabel( 'Z^* (\mum)', 'fontsize',10); legend('in','Location','southwest'); box on; grid off;
%--------------------- unifrom points in polygon
ypoints = 10000000.*(ypoints); zpoints = 10000000.*(zpoints);
%------approximate polygon
%conv_vertices = convhull(ypoints, zpoints);
%------Exact boundary
conv_vertices = boundary(ypoints, zpoints);
xq = 1:1:400;yq = 1:1:1500;[Xq, Yq] = meshgrid(xq,yq);Xq = Xq(:);Yq = Yq(:);
calculated_area = polyarea(ypoints(conv_vertices), zpoints(conv_vertices));
inside_points = inpolygon(Xq,Yq,ypoints(conv_vertices), zpoints(conv_vertices));
new_y = Xq(inside_points); new_z = Yq(inside_points);
%%--------------------- Get u values
%------- set the number of y and z entry points to the mesh size
new_y=round(new_y,0); new_z=round(new_z,0); new_y= new_y+200;new_z= new_z+750;
subplot(1,3,3);scatter(new_y,new_z,'.');
axis equal; xlim([0 400]); ylim([0 1500]);
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',10); xlabel( 'Y^* (\mum)', 'fontsize',10); ylabel( 'Z^* (\mum)', 'fontsize',10); legend('in','Location','southwest'); box on; grid off;
savefig('entry in channel')
%-------(polygon) Get U from the matrix accrording to y and z values
pointend = numel(new_y);uuu = zeros(pointend,1);
for i=1:pointend
    yy = new_y(i); zz = new_z(i); %zz = round(zz,0); yy = round(yy,0);
    uu = u(zz,yy);
    uuu(i)= uu;
end
average_u = sum(uuu)/numel(uuu);
calculated_area;
Q_cavity = 2*average_u*(calculated_area/1000000);
%% Normalized flux with channel flux
Q_channel = um*((40*150)/1000000);
%Q_channel = um*Ly*Lz;
Q_norm = Q_cavity/Q_channel;
% Create a table with the data and variable names
T = table(average_u, calculated_area, Q_cavity,Q_channel, Q_norm, 'VariableNames', { 'average_u', 'calculated_area', 'Q_cavity','Q_channel','Q_norm'} );
% Write data to text file
writetable(T, 'Flux_info.txt')
return
