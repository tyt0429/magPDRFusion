function [location,location2,location3,location4,location5,location6] = Matching(magnetmap,magnetseries,gridrelationship,X,Y,Gridsize,coordinate)
%%展示地磁基準網格圖
% Displaymap(magnetmap,X,Y,Gridsize);
%讀取測試資料

location2=[];
location3=[];
location5=[];
location6=[];
count=0;
for point=size(magnetseries,2):size(magnetseries,2)
    if(gridrelationship(point,1)<X&&gridrelationship(point,2)<Y)
        windows=point; 
    else
        count=count+1%計算超出邊界幾次
        windows=point-count;
    end
    submagnetseries=magnetseries(:,point-(windows-1):point);%仿真序列
    magnetmap2=[];
    slide=0;
    for x=1:X
        for y=1:Y
            mapindex= gridrelationship(1:point,:)+[x,y];
            mapindex2=(mapindex(:,1)-1)*Y+mapindex(:,2);
            if((mapindex(:,1)<=X)&(mapindex(:,2)<=Y))
                slide=slide+1;                       %可滑動格數
                magnetmap2=[magnetmap(:,mapindex2)];
                %dtw
                distv(slide)= dtw(submagnetseries(1,:),magnetmap2(1,:));
%                 disth(slide)= dtw(submagnetseries(2,:),magnetmap2(2,:));
%                 distt(slide)= dtw(submagnetseries(3,:),magnetmap2(3,:));
                templocation(:,slide)=coordinate(:,(x-1)*Y+y);
                %cor
                v=corrcoef(submagnetseries(1,:),magnetmap2(1,:));
                cov(slide)=v(1,2);
%                 h=corrcoef(submagnetseries(2,:),magnetmap2(2,:));
%                 coh(slide)=h(1,2);
%                 t=corrcoef(submagnetseries(3,:),magnetmap2(3,:));
%                 cot(slide)=t(1,2);
            end
        end
    end
     %dtw
    [sortdistv,index1]=sort(distv);
%     [sortdisth,index2]=sort(disth);
%     [sortdistt,index3]=sort(distt);
%     [sortdismix,index4]=sort(distv*2+disth);
   location(:,point)=(templocation(:,index1(1)))+(fliplr(gridrelationship(point,:))*Gridsize)';
%     location2(:,point)=(templocation(:,index2(1)))+(fliplr(gridrelationship(point,:))*Gridsize)';
%     location3(:,point)=(templocation(:,index3(1)))+(fliplr(gridrelationship(point,:))*Gridsize)';
    %cor
    [sortcov,index5]=sort(cov);
%     [sortcoh,index6]=sort(coh);
%     [sortcot,index7]=sort(cot);
%     [sortcomix,index8]=sort(cov*2+coh);
   location4(:,point)=(templocation(:,index5(slide)))+(fliplr(gridrelationship(point,:))*Gridsize)';
%     location5(:,point)=(templocation(:,index6(slide)))+(fliplr(gridrelationship(point,:))*Gridsize)';
%     location6(:,point)=(templocation(:,index7(slide)))+(fliplr(gridrelationship(point,:))*Gridsize)';
    clear index1 index2 index3 index4 index5 index6 index7 index8 templocation distv disth distt cov coh cot mapindex mapindex2

end

% %% 次序展示
% figure('name','matchingV');
% %繪製序列匹配圖cor
% plot(location4(1,2:size(magnetseries,2)),location4(2,2:size(magnetseries,2)),'-o','Color','g');
% axis equal
% grid on;
% hold on;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% 
% 
% %繪製序列匹配圖dtw
% 
% plot(location(1,2:size(magnetseries,2)),location(2,2:size(magnetseries,2)),'-s','Color','red');
% axis equal
% grid on;
% hold on;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% 
% 
% magnetmap3=reshape(coordinate,2,Y,X);
% for ROW=1:X
% plot(magnetmap3(1,1:Y,ROW),magnetmap3(2,1:Y,ROW),'*','Color','blue');
% axis equal
% grid on;
% hold on;
% end
% legend('corrcoef','dtw','datum');



%% 次序展示
% figure('name','matchingM');
% %繪製序列匹配圖cor
% plot(location6(1,2:size(magnetseries,2)),location6(2,2:size(magnetseries,2)),'-o','Color','g');
% axis equal
% grid on;
% hold on;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% 
% 
% %繪製序列匹配圖dtw
% 
% plot(location3(1,2:size(magnetseries,2)),location3(2,2:size(magnetseries,2)),'-s','Color','red');
% axis equal
% grid on;
% hold on;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% 
% 
% magnetmap3=reshape(coordinate,2,Y,X);
% for ROW=1:X
% plot(magnetmap3(1,1:Y,ROW),magnetmap3(2,1:Y,ROW),'*','Color','blue');
% axis equal
% grid on;
% hold on;
% end
% legend('corrcoef','dtw','datum');
% 
% 
% 
% %% 次序展示
% figure('name','matchingH');
% %繪製序列匹配圖cor
% plot(location5(1,2:size(magnetseries,2)),location5(2,2:size(magnetseries,2)),'-o','Color','g');
% axis equal
% grid on;
% hold on;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% 
% 
% %繪製序列匹配圖dtw
% 
% plot(location2(1,2:size(magnetseries,2)),location2(2,2:size(magnetseries,2)),'-s','Color','red');
% axis equal
% grid on;
% hold on;
% xlabel('X(m)','FontWeight','bold','FontSize',14);
% ylabel('Y(m)','FontWeight','bold','FontSize',14);
% 
% 
% magnetmap3=reshape(coordinate,2,Y,X);
% for ROW=1:X
% plot(magnetmap3(1,1:Y,ROW),magnetmap3(2,1:Y,ROW),'*','Color','blue');
% axis equal
% grid on;
% hold on;
% end
% legend('corrcoef','dtw','datum');