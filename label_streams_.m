clc; clear all; close all; tic;
%-----------------
cavityheight = 150 ; transparancy = 0.1; cavityin = 99; cavityout = 310; %set cavity inlet on x axis %set cavity outlet on x axis %
%-----------------
%set a parula color map
halfcavityheight = cavityheight/2;zrange=(0:halfcavityheight);ran=range(zrange); min_val=min(zrange);max_val=max(zrange);yfloor=floor(((zrange-min_val)/ran)*63)+1; co=zeros(20,3);p=colormap; %data to be plotted %finding range of data %finding maximum value of data %finding minimum value of data
for i=1:halfcavityheight
  a=yfloor(i);co(i,:)=p(a,:);stem3(i,i,zrange(i),'Color',co(i,:)); hold on
end
coll=flip(co);col=[co;coll];
%% 1- get streamlines from comsol text file (after deleting info lines and retain only numbers
% file name
load u001.txt; x = u001(:,1); y = u001(:,2); z = u001(:,3);st = u001(:,4); x = round(x,4); y = round(y,4); z = round(z,4);st = round(st,4); numberofstreamlines = st(end); toc;% name columns array  
% set streamlines in separate vecotrs x y z %fullstreamline = zeros(0,1); % set array to count number of streamlines
tic; i = 1; cont = true; cummIndex = 1; % parameters for while and for loop
while cont
    full=find(st==(i-1)); streamline_length = numel(full); % streamline length % for sstreamline start with 0 (i-1)
    for j = 1:streamline_length
        savedx(i,j) = x(j+cummIndex-1); savedy(i,j) = y(j+cummIndex-1);savedz(i,j) = z(j+cummIndex-1); %save the x y z points in matrix
    end
    cummIndex = cummIndex + (streamline_length); i=i+1;%cummIndex % 
    if i == numberofstreamlines; %for full streamline
        break
    end
end
toc; numberofstreamlines
%% 2- get data one side of the symmetry line (half of the cavity) and eliminate short non-circulating streamlines
tic;[m,n] = size(savedx);
for ii=1:m;
    A = savedx(ii,:); B = savedy(ii,:);C = savedz(ii,:); A(A==0) = []; numA=numel(A); B(B==0) = []; C(C==0) = []; 
        if numA <1000 % set zero to short non-circulating strealines 
            %&& C(50) > halfcavityheight %|  check streamlines plot % numstreaminsidecavity < 100  &&
            savedx(ii,:) = 0; savedy(ii,:) = 0; savedz(ii,:) = 0;
        end
        if ii == numberofstreamlines
               break
        end
end
toc;
%delete rows of zeros of the other half
x_half = savedx(any(savedx,2),:); y_half = savedy(any(savedy,2),:); z_half = savedz(any(savedz,2),:); [q,r] = size(x_half);tic; 
size(x_half)
%% 3- set conditions to separate circulating streamlines and non-circulating streamlines
figure('Name','non-circulatingstreamlines');
for ii=1:q;
    Aa = x_half(ii,:); Bb = y_half(ii,:);Cc = z_half(ii,:); Cnodeci = round(Cc(100),0); Aa(Aa==0) = []; numAa=numel(Aa); Bb(Bb==0) = []; Cc(Cc==0) = []; %B = savedy(:,numA);%C = savedz(:,numA);
    circulating_repeats = numel(unique(Aa))/numA; % for circulating streams
    in_cavity = find(Bb>180); numstreaminsidecavity=numel(in_cavity); % for saturated cavity flow - pick streamlines above channel streamline
            if numstreaminsidecavity < 10 ; %&& numAa <1500 ; %|| num_above < 1 || circulating_repeats == 1;%&& 
               x_half(ii,:) = 0; y_half(ii,:) = 0; z_half(ii,:) = 0;
               plot3(Aa,Bb(1:numAa),Cc(1:numAa),'Color', col(Cnodeci,:)); view([25 70]); xlim([50 350]); ylim([0 250]); zlim([75 150]); daspect([3 2.5 2]); alpha(transparancy);  xlabel( 'X (\mum)','fontsize',15); ylabel( 'Y (\mum)', 'fontsize',15 ); zlabel( 'Z (\mum)', 'fontsize',15); ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 1 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',17); box on; grid off ;hold on;
               zlim([0 150]);
            end
        if ii == r 
               break
        end
end
hold off
set(gca,'color',[0.85 0.85 0.85]); set(gcf, 'Color', 'w'); box off; savefig('non-circulatingstreamlines');
toc;
%% 4- Get circulating streamline plot
tic; x_cir = x_half(any(x_half,2),:);y_cir = y_half(any(y_half,2),:);z_cir = z_half(any(z_half,2),:); [q,r] = size(x_cir);%delete rows of zeros of non circulating flow
upregion = find(abs(z_cir(100))<40); 
z_cir(upregion)= []; y_cir(upregion) = []; x_cir(upregion) = [];
size(x_cir)
figure('Name','circulatingstreamlines');
for ii=1:q;
    AA = x_cir(ii,:); BB = y_cir(ii,:); CC = z_cir(ii,:); numbA = numel(AA);
    AA(AA==0) = [] ; numbAA=numel(AA); BB(BB==0) = []; CC(CC==0) = [];
    % To check the graph circulating streamlines
    CCnodeci = round(CC(100),0);
        plot3(AA,BB(1:numbAA),CC(1:numbAA),'Color', col(CCnodeci,:)); view([25 70]); xlim([50 350]); ylim([0 250]); zlim([75 150]); daspect([3 2.5 2]); alpha(transparancy);  xlabel( 'X (\mum)', 'fontsize',15 ); ylabel( 'Y (\mum)', 'fontsize',15 ); zlabel( 'Z (\mum)', 'fontsize',15); ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 1 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',17); box on; grid off;hold on; 
        zlim([0 150]);
    if ii == r 
        break
    end
end
savefig('circulatingstreamlines');toc;
%% 5- get streamline value of x and y at enrtry and exit
tic; 
figure('Name','entry and exit');
sz = 80;
Entry = find(x_cir>(cavityin) & x_cir<(cavityin+1)); Entry_ypoints = y_cir(Entry); Entry_zpoints = z_cir(Entry); Exit = find(x_cir>cavityout & x_cir<(cavityout+1)); Exit_ypoints = y_cir(Exit); Exit_zpoints = z_cir(Exit); %entry = Entry(1); % take streamline entering 1-2 um before the cavity 
numel(Entry_ypoints)
scatter(Entry_ypoints,Entry_zpoints,sz,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r'); alpha(transparancy);hold on; scatter(Exit_ypoints,Exit_zpoints,sz,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);alpha(0.1); %axis equal; toc; %xlim([30 40]); ylim([40 80]);
hold on
shps = alphaShape(Entry_ypoints,Entry_zpoints,10);shps = plot(shps,'lineStyle','none'); shps.FaceColor = 'r'; alpha 0.5; hold on;shp = alphaShape(Exit_ypoints,Exit_zpoints,10); shp = plot(shp,'lineStyle','none');shp.FaceColor = 'k'; alpha 0.3; 
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',40); xlabel( 'Y^* (\mum)', 'fontsize',40); ylabel( 'Z^* (\mum)', 'fontsize',40); legend('in','out','Location','southwest'); box on; grid off;
axis equal;
savefig('entry_exit')
figure;
shps = alphaShape(Entry_ypoints,Entry_zpoints,10);shps = plot(shps,'lineStyle','none'); shps.FaceColor = 'r'; hold on;shp = alphaShape(Exit_ypoints,Exit_zpoints,10); shp = plot(shp,'lineStyle','none');shp.FaceColor = 'k'; %alpha 0.5; 
ax = gca; ax.BoxStyle = 'full'; ax.LineWidth = 2 ; set(gcf, 'Color', 'w'); set(gca,'color','w'); set(gca,'fontsize',35); xlabel( 'Y^* (\mum)', 'fontsize',35); ylabel( 'Z^* (\mum)', 'fontsize',35); legend('in','out','Location','eastoutside'); box on; grid off;
axis equal;xlim([15 40]); ylim([75 150]);
savefig('entry_exit_nodots')
save('savedx');save('savedy');save('savedz');
save('x_half');save('y_half');save('z_half');
save('x_cir');save('y_cir');save('z_cir');
save('Entry_ypoints'); save('Entry_zpoints'); save('Exit_ypoints'); save('Exit_zpoints');