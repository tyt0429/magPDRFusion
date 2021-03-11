function [] = Displaymap(magnetmap,X,Y,Gridsize)
figure('name','Gridmap');
subplot(3,1,3);
sum=reshape(magnetmap(3,:),Y,X);
imagesc([0+Gridsize/2 Gridsize*Y-Gridsize/2],[0+Gridsize/2 Gridsize*X-Gridsize/2],sum');
axis xy ;
axis equal ;
xlabel('X(m)','FontWeight','bold','FontSize',14);
ylabel('Y(m)','FontWeight','bold','FontSize',14);
title('總量','FontWeight','bold','FontSize',14);
colorbar;


subplot(3,1,2);
sum=reshape(magnetmap(2,:),Y,X);
imagesc([0+Gridsize/2 Gridsize*Y-Gridsize/2],[0+Gridsize/2 Gridsize*X-Gridsize/2],sum');
axis xy ;
axis equal ;
xlabel('X(m)','FontWeight','bold','FontSize',14);
ylabel('Y(m)','FontWeight','bold','FontSize',14);
title('水平分量','FontWeight','bold','FontSize',14);
colorbar;


subplot(3,1,1);
sum=reshape(magnetmap(1,:),Y,X);
imagesc([0+Gridsize/2 Gridsize*Y-Gridsize/2],[0+Gridsize/2 Gridsize*X-Gridsize/2],sum');
axis xy ;
axis equal ;
xlabel('X(m)','FontWeight','bold','FontSize',14);
ylabel('Y(m)','FontWeight','bold','FontSize',14);
title('垂直分量','FontWeight','bold','FontSize',14);
colorbar;

figure('name','Orientation map');
sum=reshape(magnetmap(6,:),Y,X);
imagesc([0+Gridsize/2 Gridsize*Y-Gridsize/2],[0+Gridsize/2 Gridsize*X-Gridsize/2],sum');
axis xy ;
axis equal ;
xlabel('X(m)','FontWeight','bold','FontSize',14);
ylabel('Y(m)','FontWeight','bold','FontSize',14);
title('航向','FontWeight','bold','FontSize',14);
colorbar;