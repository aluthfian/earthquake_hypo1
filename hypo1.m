clc
clear
%% observation data
BBJI=[-7.46,107.65,615,38]; %bujur,lintang,elevasi masl,tp-ts
SMRI=[-7.05,110.44,203,24]; %bujur,lintang,elevasi masl,tp-ts
JAGI=[-8.47,114.15,171,41]; %bujur,lintang,elevasi masl,tp-ts
%% WGS84 constants
WGS84_a=6378137.0; %meters
WGS84_b=6356752.3142; %meters
WGS84_e2=6.694379990141e-3;
WGS84_e=8.1819190842622e-2;
%% Cartesian coordinate of stations
%BBJI
N=WGS84_a./sqrt(1-power(WGS84_e*sind(BBJI(1)),2));
X=(N+BBJI(3))*cosd(BBJI(1))*cosd(BBJI(2));
Y=(N+BBJI(3))*cosd(BBJI(1))*sind(BBJI(2));
Z=(((1-WGS84_e2)*N)+BBJI(3))*sind(BBJI(1));
BBJI_XYZ=[X,Y,Z];
%SMRI
N=WGS84_a./sqrt(1-power(WGS84_e*sind(SMRI(1)),2));
X=(N+SMRI(3))*cosd(SMRI(1))*cosd(SMRI(2));
Y=(N+SMRI(3))*cosd(SMRI(1))*sind(SMRI(2));
Z=(((1-WGS84_e2)*N)+SMRI(3))*sind(SMRI(1));
SMRI_XYZ=[X,Y,Z];
%JAGI
N=WGS84_a./sqrt(1-power(WGS84_e*sind(JAGI(1)),2));
X=(N+JAGI(3))*cosd(JAGI(1))*cosd(JAGI(2));
Y=(N+JAGI(3))*cosd(JAGI(1))*sind(JAGI(2));
Z=(((1-WGS84_e2)*N)+JAGI(3))*sind(JAGI(1));
JAGI_XYZ=[X,Y,Z];
%% Hypocenter search
lat1=-9.2;
lat2=-8.3;
long1=110.085;
long2=111.330;
depth1=-150000;
depth2=-60000;
numdiv=21;
maxiter=5;
rng('shuffle','philox')
fig = figure;
ax1 = axes;
axis tight manual % this ensures that getframe() returns a consistent size
im=cell(maxiter+1,1);
for iter=0:maxiter
    latSrcRad=lat1+(rand(numdiv,1)*(lat2-lat1));
    longSrcRad=long1+(rand(numdiv,1)*(long2-long1));
    zSrcRad=linspace(depth1,depth2,numdiv);
    epiRMS=zeros(length(latSrcRad),length(longSrcRad),length(zSrcRad));
    BBJI_dT=zeros(length(latSrcRad),length(longSrcRad),length(zSrcRad));
    SMRI_dT=zeros(length(latSrcRad),length(longSrcRad),length(zSrcRad));
    JAGI_dT=zeros(length(latSrcRad),length(longSrcRad),length(zSrcRad));
    for idx1=1:length(latSrcRad)
        for idx2=1:length(longSrcRad)
            for idx3=1:length(zSrcRad)
                Ns=WGS84_a./sqrt(1-power(WGS84_e*sind(latSrcRad(idx1)),2));
                Xs=(Ns+zSrcRad(idx3))*cosd(latSrcRad(idx1))*cosd(longSrcRad(idx2));
                Ys=(Ns+zSrcRad(idx3))*cosd(latSrcRad(idx1))*sind(longSrcRad(idx2));
                Zs=(((1-WGS84_e2)*Ns)+zSrcRad(idx3))*sind(latSrcRad(idx1));
                % velocity src: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2006JB004712
                vpMinvs=9.7e+3;
                tptsBBJI=sqrt((BBJI_XYZ(1)-Xs)^2 +(BBJI_XYZ(2)-Ys)^2 +(BBJI_XYZ(3)-Zs)^2)/vpMinvs;
                tptsSMRI=sqrt((SMRI_XYZ(1)-Xs)^2 +(SMRI_XYZ(2)-Ys)^2 +(SMRI_XYZ(3)-Zs)^2)/vpMinvs;
                tptsJAGI=sqrt((JAGI_XYZ(1)-Xs)^2 +(JAGI_XYZ(2)-Ys)^2 +(JAGI_XYZ(3)-Zs)^2)/vpMinvs;
                BBJI_dT(idx1,idx2,idx3)=abs(tptsBBJI-BBJI(4));
                SMRI_dT(idx1,idx2,idx3)=abs(tptsSMRI-SMRI(4));
                JAGI_dT(idx1,idx2,idx3)=abs(tptsJAGI-JAGI(4));
                epiRMS(idx1,idx2,idx3)=sqrt((BBJI_dT(idx1,idx2,idx3)^2+SMRI_dT(idx1,idx2,idx3)^2+JAGI_dT(idx1,idx2,idx3)^2)/3);
            end
        end
    end
    % plotting
    lat4plot=zeros(numdiv^3,1);
    long4plot=zeros(numdiv^3,1);
    depth4plot=zeros(numdiv^3,1);
    for idxtot=1:numdiv^3
        [idxLatPlot,idxLongPlot,idxDepthPlot]=ind2sub([numdiv numdiv numdiv],idxtot);
        long4plot(idxtot)=longSrcRad(idxLongPlot);
        lat4plot(idxtot)=latSrcRad(idxLatPlot);
        depth4plot(idxtot)=zSrcRad(idxDepthPlot);
    end
    if iter<2
        scatter3(long4plot,lat4plot,depth4plot/1000,'.')
    else
        scatter3(long4plot,lat4plot,depth4plot/1000,180,'p','filled')
    end
    box on
    ax = gca;
    ax.BoxStyle = 'full';
    xlim([109,112])
    ylim([-10,-8])
    zlim([-100,0])
    title('Proses penentuan hiposenter gempa','FontSize',14)
    hipoLong=mean(longSrcRad(idx2));
    hipoLat=mean(latSrcRad(idx1));
    hipoDepth=mean(zSrcRad(idx3));
    string = {'Pusat gempa di',strcat('(',num2str(hipoLat,4),',',num2str(hipoLong,6),')'),...
        ['kedalaman ',num2str(abs(hipoDepth/1e+3),5),' km'],['ralat $\Delta t$ ',num2str(min(min(min(epiRMS))),3),' detik']};
    txt=text(109,-8.2,-15.0,string,'Color','black','FontSize',12,'Interpreter','latex','BackgroundColor','white');
    % hidden axes to put texts
    axHidden=axes('Visible','off','hittest','off'); % Invisible axes
    linkprop([ax1 axHidden],{'CameraPosition' 'XLim' 'YLim' 'ZLim' 'Position'}); % The axes should stay aligned
    set(txt,'Parent',axHidden); % Put the text in the invisible Axes
    xlabel('bujur (derajat)')
    ylabel('lintang (derajat)')
    zlabel('kedalaman (km)')
    % writing image for GIF
    drawnow
    frame = getframe(fig);
    im{iter+1} = frame2im(frame);
    pause(1)
    % creating new bound
    idxMinRMS=find(epiRMS==min(min(min(epiRMS))));
    [idx1,idx2,idx3]=ind2sub(size(epiRMS),idxMinRMS);
    idx1=unique(idx1);
    idx2=unique(idx2);
    idx3=unique(idx3);
    lat_old=[lat1,lat2];
    long_old=[long1,long2];
    depth_old=[depth1,depth2];
    lat1=latSrcRad(min(idx1))-((lat_old(2)-lat_old(1))/(numdiv-1));
    lat2=latSrcRad(max(idx1))+((lat_old(2)-lat_old(1))/(numdiv-1));
    long1=longSrcRad(min(idx2))-((long_old(2)-long_old(1))/(numdiv-1));
    long2=longSrcRad(max(idx2))+((long_old(2)-long_old(1))/(numdiv-1));
    depth1=zSrcRad(min(idx3))-((depth_old(2)-depth_old(1))/(numdiv-1));
    depth2=zSrcRad(max(idx3))+((depth_old(2)-depth_old(1))/(numdiv-1));
end
close
filename = 'hypocenter.gif'; % Specify the output file name
for idx = 1:maxiter+1
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end