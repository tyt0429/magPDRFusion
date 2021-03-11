function [gridrelationship,magnetseries]=Rasterprocessing(F,X,Y,gridsize,EG,NG,coord_g,locsfin,updatepoint)
%去找軌跡跟網格交點,用時間去內插在網格的那一段得到網格內的平均磁場
%% to raster
coordinate_x = gridsize/2:gridsize:Y-gridsize/2;
coordinate_y = gridsize/2:gridsize:X-gridsize/2;
for j=1:X
    for i=1:Y
        coordinate(1,i+(j-1)*Y) = coordinate_x(i);
        coordinate(2,i+(j-1)*Y) = coordinate_y(j);
    end
end
coordinateX = (reshape(coordinate(1,:),Y,X))';
coordinateY = (reshape(coordinate(2,:),Y,X))';
% PDR軌跡
PDRpointG = [EG(1),NG(1);coord_g];
time = [1;locsfin];%每一步時間
for i=1:length(PDRpointG)-1
    vector(i,:) = PDRpointG(i+1,:) - PDRpointG(i,:);%每一步向量
end

%展示軌跡
figure('name','trajectory');
%path1
% rectangle('Position',[0 0 39 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
%path2
% rectangle('Position',[0 0 21 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[20 1 1 3],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[21 3 9 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[29 4 1 2],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[30 5 8 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
%path3
% rectangle('Position',[0 0 4 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[4 1 6 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[10 2 6 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[16 3 6 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[22 4 6 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[28 5 5 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);
% rectangle('Position',[33 6 5 1],'FaceColor',[0.0 0.0 0.0 0.3],'EdgeColor',[0 0.0 0.0 0.3]);

mapshow(PDRpointG(:,1),PDRpointG(:,2),'Marker','.','Color', 'b');
hold on;
plot(updatepoint(2:end,1),updatepoint(2:end,2),'*','Color', 'r');

xlabel('x coordinate');
ylabel('y coordinate');
title('The trajectory of the pedestrian','FontWeight','bold','FontSize',10);
axis equal;
grid on;
grid minor;
% 計算與網格交點
clear griddata
count=0;
gridtrajectory=[];
gridtime=[];

%搜尋局部基準圖
for j = 1:X
    for i = 1:Y
        %網格中心點座標
        cx = coordinateX(j,i);
        cy = coordinateY(j,i);
        %網格四角座標形成 polygon
        x1 = [cx-gridsize/2 cx-gridsize/2 cx+gridsize/2 cx+gridsize/2 cx-gridsize/2];
        y1 = [cy-gridsize/2 cy+gridsize/2 cy+gridsize/2 cy-gridsize/2 cy-gridsize/2];
        
        %展示網格
        %mapshow(x1,y1,'Color', 'black');
        
        %計算軌跡與網格交點
        [xi,yi,ii] = polyxpoly(PDRpointG(:,1),PDRpointG(:,2),x1,y1);%ii為交於 polygon 哪一段
        %展示交點
%         mapshow(xi,yi,'Displaytype','point','Marker','o');
        %計算交點與該步起點之向量
        ivector=[xi-PDRpointG(ii(:,1),1) yi-PDRpointG(ii(:,1),2)];
        %內插交點時間
        temp1 = vecnorm(ivector')'; %該步起始向量~交點
        temp2 = vecnorm(vector(ii(:,1),:)')'; %每一步向量
        % (t2*sa+t1*(sb-sa))/sb 108-2meeting/1104ppt
        ti=((temp1)./temp2).*(time(ii(:,1)+1)-time(ii(:,1)))+time(ii(:,1));
        %以網格為單位儲存結果
        count=count+1;
        gdata(count) = struct('gridindex',{[j i]},'gridcentercoordinate',{[cx cy]},'Intersection',{[xi yi ti]});
        %輸出通過網格資料 有交點才記錄到gridtrajectory
        if (isempty(gdata(count).Intersection)~=1)
            gridtrajectory = [gridtrajectory;gdata(count).gridindex];%通過網格之索引
            passtime = sort(gdata(count).Intersection(:,3))';%每個網格通過交點時間
            passtime = [passtime(1) passtime(length(passtime))];
            %通過網格之時間
            if(length(passtime)==2)
                gridtime = [gridtime; passtime];
            else %不只通過兩個網格的問題
                if(ii(:,1)==1)%起終點
                    gridtime = [gridtime;time(1) passtime];
                else
                    gridtime = [gridtime;passtime time(length(time))];
                end
            end
        end
    end
end
% 網格關係&地磁序列萃取
gridrelationship = gridtrajectory - gridtrajectory(1,:);
%對網格內地磁序列做平均
for i=1:size(gridtime,1)
    magnetseries(:,i) = mean(F(:,round(gridtime(i,1)):1:round(gridtime(i,2))),2);
end
% figure('name','LinePDRgridseries');
%
% display=[gridrelationship+1 magnetseries(:,:)'];
% subsize=boundary;
% subimage=nan(subsize(2),subsize(1));
% subimagev=subimage;
% subimageh=subimage;
% subimaget=subimage;
%
% for i=1:size(gridrelationship,1)
%         subimagev(gridrelationship(i,1)+1,gridrelationship(i,2)+1)=display(i,3);
%         subimageh(gridrelationship(i,1)+1,gridrelationship(i,2)+1)=display(i,4);
%         subimaget(gridrelationship(i,1)+1,gridrelationship(i,2)+1)=display(i,5);
% end
%
% subplot(3,1,3);
% h=imagesc([0+gridsize/2 gridsize*subsize(1)-gridsize/2],[0+gridsize/2 gridsize*(subsize(2))-gridsize/2],subimaget);
% set(h,'alphadata',~isnan(subimaget));
% axis xy ;
% axis equal;
%
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% title('總量','FontWeight','bold','FontSize',14);
% h = colorbar;
%
% h.Label.String = 'Magnitude(μT)';
%
% subplot(3,1,2);
% h2=imagesc([0+gridsize/2 gridsize*subsize(1)-gridsize/2],[0+gridsize/2 gridsize*(subsize(2))-gridsize/2],subimageh);
% set(h2,'alphadata',~isnan(subimageh));
% axis xy ;
% axis equal ;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% title('水平分量','FontWeight','bold','FontSize',14);
% h = colorbar;
%
% h.Label.String = 'Magnitude(μT)';
%
% subplot(3,1,1);
% h3=imagesc([0+gridsize/2 gridsize*subsize(1)-gridsize/2],[0+gridsize/2 gridsize*(subsize(2))-gridsize/2],subimagev);
% set(h3,'alphadata',~isnan(subimagev));
% axis xy ;
% axis equal ;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% title('垂直分量','FontWeight','bold','FontSize',14);
% h = colorbar;
%
% h.Label.String = 'Magnitude(μT)';
