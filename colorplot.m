%Loading output file

dat = load('outputs/outse');

[height, length] = size(dat);
%Finding the intervals
dx=dat(((height)^.5+1),1)-dat(1,1);
dy=dat(2,2)-dat(1,2);

%creating a vector x with all the needed divisions 
x=[dat(1,1):dx:dat((height),1)];
%Fixing the floating-point arithmetic issue of 1x399
x=[x,dat(height,1)];

%same with y
y=[dat(1,2):dy:dat(sqrt(height),2)];
y=[y,dat(sqrt(height),2)];

%create a new matrix to have values for potential
potl=zeros(sqrt(height),sqrt(height));
for ix=1:sqrt(height);
        for iy=1:sqrt(height);
            potl(ix,iy)=dat(sqrt(height)*(ix-1)+iy, 3);
        end
end

%Mesh & Contour plot
figure
mesh (x, y, potl','FaceAlpha','0.5')
xlabel('x in AU')
ylabel('y in AU')
zlabel('Effective Potential')
title('Potential in a gravitational system')
axis tight
colorbar
zlim([-3 max(max(potl))])
colormap(jet)    % change color map

maximum = max(max(potl));
[px,py]=find(potl==maximum)

figure
hold on
imagesc(x,y',potl',[-1.6 max(max(potl))])
xlabel('x in AU')
ylabel('y in AU')
grid on
axis('equal')
set(gca,'YDir','normal')
colormap('Jet')

%Plotting maximum point/points
for i=1:size(px)
    scatter(x(px(i)), y(py(i)),'w','filled')
end

% potlex=potl;
% potlex(:,1:200)=-100;
% maximum = max(max(potlex));
% [px,py]=find(potlex==maximum)
% x(px)
% y(py)

%scatter(x(px), y(py),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1)

%Orbit as needed
viscircles([0 0],3,'LineStyle',':', 'Color', 'k', 'LineWidth', 0.1)

%PLotting the position of the planets/sun
%Jupiter
scatter(7.719208023413598E-01,-5.159679042599656, 40, 'MarkerEdgeColor',[0.7 0.2 0.3],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',1)
%Saturn
scatter(-8.487285193590072E+00,3.763807729501108, 40, 'MarkerEdgeColor',[0.9 0.3 0.09],'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',1)
%Sun/M1
scatter(0, 0, 50,'MarkerEdgeColor',[1 0.8 0.2],'MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',1)
%"M2"
scatter(1, 0, 20,'MarkerEdgeColor',[0 0.5 0.8],'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1)
%"M3"
scatter(0, 3,'MarkerEdgeColor',[0.7 0.2 0.3],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',1)

hold off;




