clc; clear all; close all; tic;
%% 1 get streamlines from comsol text file (after deleting information lines from streamline file and retain only numbers)
%put file name in all 5 places next line
load seg_streams_example.txt; x = seg_streams_example(:,1); y = seg_streams_example(:,2); z = seg_streams_example(:,3);st = seg_streams_example(:,4); x = round(x,4); y = round(y,4); z = round(z,4);st = round(st,4); toc;% name columns array
numberofstreamlines = st(end);
%% set streamlines in separate vecotrs x y z 
tic;
i = 1; cont = true; cummIndex = 1;% parameters for while and for loop
while cont
    full=find(st==(i-1)); streamline_length = numel(full); % streamline length % for streamline start with 0 (i-1)
    for j = 1:streamline_length
        savedx(i,j) = x(j+cummIndex-1); savedy(i,j) = y(j+cummIndex-1);savedz(i,j) = z(j+cummIndex-1); %save the x y z points in matrix
    end
    cummIndex = cummIndex + (streamline_length);%cummIndex %
    i=i+1;
    if i == numberofstreamlines; %for full streamline
    %if i == 1000 % for initial testing
        break
    end
end
toc;
%% 2- set conditions to get circulating streamlines and zero non-circulating
tic;
transparancy = 0.2;
[m,n] = size(savedx); 
%Segragate non-circulating streamlines form circulating streamline in the cavity
figure('Name','non-circulating and circulating streamlines')
for ii=1:m;
    A = savedx(ii,:); B = savedy(ii,:);C = savedz(ii,:);  A(A==0) = [];  numA=numel(A); B(B==0) = []; C(C==0) = []; %numA=numel(A);
            if numA <1900 
             savedx(ii,:) = 0; savedy(ii,:) = 0; savedz(ii,:) = 0; 
                    plot3(A,B(1:numA),C(1:numA),'Color', [204/255 237/255 247/255],'LineWidth',1); hold on;%  xlim([100 450]); ylim([0 300]); alpha(0.001); hold on;
            else
                    plot3(A,B(1:numA),C(1:numA),'Color', 'r','LineWidth',0.5); hold on; view([-100 -50 150]);
            end
        if ii == numberofstreamlines %for full streamline
        %if ii == 50 % for initial testing
               break
        end
end
view([-180 -150 350]); daspect([0.8 0.5 0.9]);
savefig('non-circulatingstreamlines')
toc;
%% 3- Get circulating streamline plot
tic; x_cir = savedx(any(savedx,2),:);y_cir = savedy(any(savedy,2),:);z_cir = savedz(any(savedz,2),:); %delete rows of zeros of non circulating flow
[q,r] = size(x_cir);
%% get streamline value of x and y at enrtry and exit
tic;
cavityin = 199; cavityout = 401 ; 
figure('Name','entry and exit');
sz = 80;
Entry = find(x_cir>(cavityin) & x_cir<(cavityin+1)); Entry_ypoints = y_cir(Entry); Entry_zpoints = z_cir(Entry); Exit = find(x_cir>cavityout & x_cir<(cavityout+1)); Exit_ypoints = y_cir(Exit); Exit_zpoints = z_cir(Exit); % take streamline entering 1-2 um before the cavity 
numel(Entry_ypoints)
scatter(Entry_ypoints,Entry_zpoints,sz,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); alpha(transparancy);hold on; scatter(Exit_ypoints,Exit_zpoints,sz,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);alpha(0.1);
hold on
shps = alphaShape(Entry_ypoints,Entry_zpoints,10);shps = plot(shps,'lineStyle','none'); shps.FaceColor = 'r'; alpha 0.9; hold on;shp = alphaShape(Exit_ypoints,Exit_zpoints,10); shp = plot(shp,'lineStyle','none');shp.FaceColor = 'k'; alpha 0.9; 
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color','w'); set(gca,'color',[204/255 237/255 247/255]); set(gca,'fontsize',40); xlabel( 'Y^* (\mum)', 'fontsize',40); ylabel( 'Z^* (\mum)', 'fontsize',40); legend('in','out','Location','eastoutside'); box on; grid off;
savefig('entry_exit')
figure;
shp = alphaShape(Exit_ypoints,Exit_zpoints,10); shp = plot(shp,'lineStyle','none'); shp.FaceColor = 'k'; alpha 0.9; hold on;shps = alphaShape(Entry_ypoints,Entry_zpoints,10);shps = plot(shps,'lineStyle','none'); shps.FaceColor = 'r';alpha 0.7; ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',35); xlabel( 'Y (\mum)', 'fontsize',35); ylabel( 'Z (\mum)', 'fontsize',35); legend('in','out','Location','eastoutside'); box on; grid off;
axis equal;xlim([20 40]); ylim([75 150]);set(gca,'color',[204/255 237/255 247/255]);set(gca, 'XAxisLocation', 'top');set(gca, 'YDir','reverse');
title('Entrance and Exit of Recirculating Cavity Flow');ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',17); xlabel( 'Y (\mum)', 'fontsize',20 ); ylabel( 'Z (\mum)', 'fontsize',20 ); legend('in','out','Location','eastoutside'); box on; set(gca,'color',[204/255 237/255 247/255]); grid off;toc;
save('Entry_ypoints'); save('Entry_zpoints');
save('Exit_ypoints'); save('Exit_zpoints');
savefig('Entry_exit')
