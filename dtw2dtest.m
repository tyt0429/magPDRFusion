clc;clear all;close all;
%% Parameter
% PDR height
height = 1.62;
% weight
mag_w = 0.5;
gyro_w = 0.5;
% 取樣區間
time_interval = 0.1;
% 讀取地磁基準網格圖
magnetmap = csvread('0121Line_GRID_mean.csv');

if ~exist('dirName','var')
    dirName='C:\Users\a5609\Desktop';
elseif dirName==0
    dirName='C:\Users\a5609\Desktop';
end
addpath functions\;
indicator = [];
indicator2 = [];

%% 設定網格座標
%長寬
X=7;
Y=38;
% 網格大小
Gridsize=1;
for i=1:Y
    coordinatex(i)=Gridsize/2+(i-1)*Gridsize;
end
for j=1:X
    coordinatey(j)=Gridsize/2+(j-1)*Gridsize;
end
for j=1:X
    for i=1:Y
        coordinate(1,i+(j-1)*Y)=coordinatex(i);
        coordinate(2,i+(j-1)*Y)=coordinatey(j);
    end
end
%% Read File
tmp=dirName;
[FileName0, dirName0] = uigetfile({'*.csv;*.txt','csv,txt'},'Select AndroSensor DATA file (csv) ',dirName);
FileName1=char(FileName0);
dirName1=char(dirName0);
for k=1:1
    FileName=FileName1(k,:);
    dirName=dirName1(1,:);
    if FileName == 0
        dirName=tmp;
        return;
    end
    ds=datastore([dirName FileName]);
    disp(['Read file: ' FileName]);disp(['From:      ' dirName]);
    ds.Delimiter={';' ','};
    ds=CheckFieldName(ds);
    ds.SelectedVariableNames=...
        {'Acc_x','Acc_y','Acc_z','Gra_x','Gra_y','Gra_z','Gyr_x','Gyr_y','Gyr_z','Mag_x','Mag_y','Mag_z','Ori_z','YYYY_MO_DDHH_MI_SS_SSS','TimeSinceStartInMs'};
    data=readall(ds);
    %% gravity projection
    G=[data.Gra_x data.Gra_y data.Gra_z];%手機座標系 重力向量
    M=[data.Mag_x data.Mag_y data.Mag_z];%手機座標系 地磁向量
    i=size(G,1);
    for j=1:i
        a(j)=dot(G(j,:),M(j,:))/norm(G(j,:));%地磁投影至重力方向 純量
        V(j,:)=a(j)* (G(j,:)/norm(G(j,:)));%投影後向量----垂直分量
        H(j,:)=M(j,:)-V(j,:);%水平分量
        b(j)=norm(H(j,:));%水平分量強度
        c(j)=norm(M(j,:));%地磁總強度
        a(j)=-a(j);
    end
    D=[a;b;c];
    EG2=smoothdata(D','movmedian',0.1/(time_interval)*5);%移動平均
    F=EG2';
    %     figure('name','MAG');
    %     plot(D(1,:));
    %     hold on;
    %     plot(F(1,:));
    %     xlabel('Time series')
    %     ylabel('Magnitude(μT)')
    %     title('The Magnitude of MAGNET','FontWeight','bold','FontSize',10)
    dlmwrite("magseries.csv",D);
end
%% PDR
acc = [data.Acc_x data.Acc_y data.Acc_z];
gra = [data.Gra_x data.Gra_y data.Gra_z];
time=data.TimeSinceStartInMs;
%% Step 1 : step detection
acc(:,1)=smoothdata(acc(:,1),'movmedian',0.1/(time_interval)*2);%移動平均
acc(:,2)=smoothdata(acc(:,2),'movmedian',0.1/(time_interval)*2);%移動平均
acc(:,3)=smoothdata(acc(:,3),'movmedian',0.1/(time_interval)*2);%移動平均
for i = 1:size(acc,1)
    mag(i,1) = sqrt(acc(i,1)^2+acc(i,2)^2+acc(i,3)^2)-sqrt(gra(i,1)^2+gra(i,2)^2+gra(i,3)^2);
end
[pks,locs] = findpeaks(mag,'MINPEAKHEIGHT',0.4,'MINPEAKDISTANCE',0.1/(time_interval)*3);
[pks2,locs2] = findpeaks(-mag,'MINPEAKHEIGHT',0.4,'MINPEAKDISTANCE',0.1/(time_interval)*3);
% pks 對應峰值
% locs 對應峰值位數
% MINPEAKHEIGHT 峰值最小高度
% MINPEAKDISTANCE 兩峰值間最小間隔數
locs2O=[0 ;locs2];
%波峰波谷判別
for i = 2:size(locs2O)
    temp=[];
    for j = 1:size(locs)
        if(locs(j)<=locs2O(i)&&locs(j)>=locs2O(i-1))
            temp=[temp;locs(j)];
        end
    end
    if(~isempty(temp))
        dis=locs2O(i)-temp;
        [sortdis,index]=sort(dis);
        wave(i-1,:)=[locs2O(i) temp(index(1))];
    end
end
wave(all(wave==0,2),:)=[];
[locs3] = Zeros_finding_Rowdata_final(mag,data.TimeSinceStartInMs,100,100);%找零點
locs3=locs3';
for k=1:size(wave,1)
    temp2=[];
    for l=1:size(locs3)
        if(locs3(l)>=wave(k,1))
            temp2=[temp2;locs3(l)];
        end
    end
    if(~isempty(temp2))
        [sorttemp2,index2]=sort(temp2);
        locsfin(k)=sorttemp2(1);
    end
end
locsfin=locsfin';
% figure('name','ACC');
% plot(mag);
% hold on;
% plot(locs,pks,'r *');
% plot(locs2,-pks2,'b *');
% plot(locs3,mag(locs3),'g *');
% plot(locsfin,mag(locsfin),'k o');
% legend('Acceleration','peak','valley','zero','result');
% xlabel('Time series')
% ylabel('Magnitude (m/s^2)')
% title('The Magnitude of Acceleration ','FontWeight','bold','FontSize',10);
% set(gca,'FontWeight','bold','fontsize',26);
fprintf('The total step are %d\n',size(locsfin,1))
%% Step 2 : step length
for i=2:size(locsfin,1)
    period(i-1)=(locsfin(i)-locsfin(i-1))*time_interval;
end
period=[(wave(2,2)-wave(1,2))*time_interval period];

SF =1./period;

step = 0.7 + 0.371 * (height - 1.75) + 0.227 * ((SF - 1.79) * height / 1.75);
cdistance(1)=step(1);
for i=2:size(locsfin,1)
    cdistance(i)=cdistance(i-1)+step(i);
end
distance =cdistance(size(locsfin,1));
fprintf('The total distance are %f\n',distance);
%% Step 4 : heading estimation ( gyro & Magnetometer )
gyro = [data.Gyr_x data.Gyr_y data.Gyr_z];
orientationmap=magnetmap(6,:);
north_mag=atan2d(smoothdata(M(:,2),'movmedian',0.1/(time_interval)*10),smoothdata(M(:,1),'movmedian',0.1/(time_interval)*10));
for i=1:size(north_mag,1)
    if( north_mag(i)<0)
        north_mag(i)=north_mag(i)+360;
    end
    if( north_mag(i)>360)
        north_mag(i)=north_mag(i)-360;
    end
    north_mag(i)=north_mag(i)-90;
end
%從第一筆資料開始計算
start=0;
[g_d_az]=gyroheading(gyro,time_interval,start);

GDAZ=[];
MAGAZ=[];
%起始點
NG = 0.5;
EG = 0.5;
updatepoint=[EG NG];

% NM = 0.5;
% EM = 0.5;
j=0;
update=0;
lastnnindex=1;
updatelocsfin=1;
clear coord_g
for i = 2:1:size(mag,1)+1
    
    
    if i == locsfin(find(locsfin==i),1)
        j=j+1;
        NG(i,1) = NG(i-1,1) + step(j) * cosd(g_d_az(i-1,1));
        EG(i,1) = EG(i-1,1) + step(j) * sind(g_d_az(i-1,1));
        GDAZ=[GDAZ;g_d_az(i-1,1)];
        
    else
        NG(i,1) = NG(i-1,1);
        EG(i,1) = EG(i-1,1);
        %         NM(i,1) = NM(i-1,1);
        %         EM(i,1) = EM(i-1,1);
    end
    if EG(i,1) > 38
        EG(i,1) = 38;
    elseif EG(i,1) < 0
        EG(i,1) = 0;
    end
    
    if NG(i,1) > 7
        NG(i,1) = 7;
    elseif NG(i,1) < 0
        NG(i,1) = 0;
    end
    
    if((i-update)==110)
        nnindex=find(i-locsfin(:)<=5);%最接近一步
        sublocsfin=locsfin(lastnnindex:nnindex(1));
        lastnnindex = nnindex(1);
        lastsublocsfin = sublocsfin;
        coord_g = [];
        for k=1:size(sublocsfin,1)-1
            coord_g(k,:) = [EG(sublocsfin(k)) ,NG(sublocsfin(k))];            
        end

        %軌跡網格化
        %        [gridrelationship,sgtitle(str,'FontWeight','bold','fontsize',26);]=Rasterprocessing(F,X,Y,Gridsize,EG(locsfin(updatelocsfin):locsfin(size(sublocsfin,1))),NG(locsfin(updatelocsfin):locsfin(size(sublocsfin,1))),coord_g(updatelocsfin:size(sublocsfin,1)-1,:),sublocsfin(updatelocsfin:size(sublocsfin,1)));
        [gridrelationship,magnetseries]=Rasterprocessing(F,X,Y,Gridsize,EG,NG,coord_g,sublocsfin,lastsublocsfin,updatepoint);
        %序列匹配
        [location,location2,location3,location4,location5,location6,indic, indic2] =Matching(magnetmap,magnetseries,gridrelationship,X,Y,Gridsize,coordinate);
        %         [mag_az]=magnetheading(north_mag,sublocsfin,nnindex,location4,coordinate,orientationmap);
        %位置更新(最接近i的一步的時刻到i)
        updatepoint=[updatepoint;location(1,size(location,2)) location(2,size(location,2))];
        if sqrt((updatepoint(end,1)-EG(end,1))^2+(updatepoint(end,2)-NG(end,1))^2) < 3
            EG(locsfin(size(sublocsfin,1)):i) = location(1,size(location,2));
            NG(locsfin(size(sublocsfin,1)):i) = location(2,size(location,2));
            %地磁扭曲補償與航向計算
            [mag_az]=magnetheading(north_mag,sublocsfin,nnindex,location,coordinate,orientationmap);
        else
            location(1,size(location,2)) = EG(end,1);
            location(2,size(location,2)) = NG(end,1);
            [mag_az]=magnetheading(north_mag,sublocsfin,nnindex,location,coordinate,orientationmap);
  
        end
         %         indicator = [indicator, indic];
        %         indicator2 = [indicator2, indic2];
        %         EG(locsfin(size(sublocsfin,1)):i) = location4(1,size(location4,2));
        %         NG(locsfin(size(sublocsfin,1)):i) = location4(2,size(location4,2));
        %航向更新
        %從i==100(最接近的一步)開始
        start=locsfin(size(sublocsfin,1));
        g_d_az(start)=(mag_az*mag_w+g_d_az(start)*gyro_w);
        [g_d_az]=gyroheading(gyro,time_interval,start,g_d_az);
        %更新時間點
        update=i;
        updatelocsfin=nnindex(1);
        
        
    end
end


for k=1:size(locsfin,1)-1
    coord_g(k,:) = [EG(locsfin(k)) ,NG(locsfin(k))];
    %             coord_m(j,:) = [EM(sublocsfin(j)) ,NM(sublocsfin(j))];
end

%% original pdr
north_mag2 =90
for i = 1:size(gyro,1)
    if i == 1
        g_d2(i,1) = -gyro(i,3) * time_interval * 180 / pi;
        g_d_az2(i,1) = g_d2(i,1) + north_mag2;
        continue
    end
    
    g_d2(i,1) = -gyro(i,3) * time_interval * 180/pi;
    
    if g_d2(i,1) + g_d_az2(i-1,1) > 360
        g_d_az2(i,1) = g_d2(i,1) + g_d_az2(i-1,1) - 360;
        continue
    end
    
    if g_d2(i,1) + g_d_az2(i-1,1) < 0
        g_d_az2(i,1) = g_d2(i,1) + g_d_az2(i-1,1) + 360;
        continue
    end
    
    g_d_az2(i,1) = g_d2(i,1) + g_d_az2(i-1,1);
end
NG2 = 0.5;
EG2 = 0.5;
j=0;
for i = 2:1:size(mag,1)+1
    
    
    if i == locsfin(find(locsfin==i),1)
        j=j+1;
        NG2(i,1) = NG2(i-1,1) + step(j) * cosd(g_d_az2(i-1,1));
        EG2(i,1) = EG2(i-1,1) + step(j) * sind(g_d_az2(i-1,1));
    else
        NG2(i,1) = NG2(i-1,1);
        EG2(i,1) = EG2(i-1,1);
    end
    
    if EG2(i,1) > 38
        EG2(i,1) = 38;
    elseif EG2(i,1) < 0
        EG2(i,1) = 0;
    end
    
    if NG2(i,1) > 7
        NG2(i,1) = 7;
    elseif NG2(i,1) < 0
        NG2(i,1) = 0;
    end
    
end
for i=1:size(locsfin,1)
    coord_g2(i,:) = [EG2(locsfin(i)) ,NG2(locsfin(i))];
end
%% fin

[gridrelationship,magnetseries]=Rasterprocessing(F,X,Y,Gridsize,EG,NG,coord_g,locsfin,updatepoint);

mapshow(coord_g2(:,1),coord_g2(:,2),'Marker','O','Color', 'g');
xlabel('X(m)');ylabel('Y(m)');
axis equal; xlim([-3 42]); ylim([-3.5 8]);
legend('result','position update','original PDR');
set(gca,'FontWeight','bold','fontsize',14)