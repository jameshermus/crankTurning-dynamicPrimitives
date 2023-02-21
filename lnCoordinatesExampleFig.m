
function [] = lnCoordinatesExampleFig()

% [Circle, 0 deg, 90 deg, 180 deg, 270] on plot
Phi = [0, 0, (1/4)*pi,0, (7/4)*pi ]; %0.1:0.1:2*pi-0.1;
A1 =  0.25.*[sqrt(0.5),   1, sqrt(0.5), 0.5, sqrt(0.5)];
A2 =  0.25.*[sqrt(0.5), 0.5, sqrt(0.5),   1, sqrt(0.5)];
figure;
for ii = 1:length(Phi)
    
    theta = 0:0.01:2*pi;
    phi = Phi(ii);%pi/8;
    a1 = A1(ii);
    a2 = A2(ii);
    x = a1*cos(theta)';% + 0.01*randn(size(theta));
    y = a2*sin(theta + phi)';% + 0.01*2*pi*randn(size(theta));
%     figure;plot(x,y); axis equal;

    % Area normalize
    [V,D]= eig((1/length(theta))*[x;y]*[x;y]');
    
    % Sort eigenvaulue dex
%     [tmp,dexSort] = sort(diag(D),'descend');
%     V = V(:,dexSort);
%     th1 = atan2(V(2,1),V(1,1));
%     th2 = atan2(V(2,2),V(1,2));
%     
%     a = get_l(th1,x,y);
%     b = get_l(th2,x,y);
%     
%     A = 2*a*b;
%     x = x/A;
%     y = y/A;
    
    for i = 1:4
        theta_d = (i-1)*(pi/4);
        l(i) = get_l(theta_d,x,y);
    end
    
    lx = l(1);
    ly = l(3);
    
    lx_prime = l(2);
    ly_prime = l(4);
    
    ln0(ii) = log(lx/ly);
    ln45(ii) = log(lx_prime/ly_prime);
    
    [V,D]= eig((1/length(theta))*[x;y]*[x;y]');
    A2(ii) = 2*D(1,1)*D(2,2);
    
%     figure(1); plot(x,y,'linewidth',2.5); hold on;
%     xlim([-1 1]); ylim([-1 1]); axis equal;
    
    plot(ln0(ii) + x,ln45(ii) + y,'-','linewidth',2.5); hold on;
    axis equal;  xlim([-1.1 1.1]); ylim([-1.25 1.25]); ylabel('ln(r_{45})'); xlabel('ln(r_{0})');
    xticks([-1 1]); yticks([-1 1]);
    set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin','fontsize',16,'linewidth',2,'box', 'off');
    
end

% figure(2); plot(ln0,ln45,'.','markersize',20); hold on;
% axis equal; % xlim([-1 1]); ylim([-1 1]);
% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');

% figure(3); histogram(ln0);
% figure(4); histogram(ln45);

disp('test');
end

function [l] = get_l(theta_d,x,y)

% Assume mean zero
if((size(x,2) ~= 1) || (size(y,2) ~= 1))
    error('transpose problem with get_l input');
end

% Computer theta
theta_hat = atan2(y, x);
indxs = find(theta_hat(:,1)<0);
theta_hat(indxs,1) = theta_hat(indxs) + 2*pi;

delta_theta = (2*pi)/200;
dex0 = find( abs(wrapToPi(theta_hat-theta_d     )) < delta_theta );
dex1 = find( abs(wrapToPi(theta_hat-(theta_d+pi))) < delta_theta );

l = sqrt( ( mean(x(dex1)) - mean(x(dex0)) )^2 + ( mean(y(dex1)) - mean(y(dex0)) )^2 );

% Check angle for atan2 and theta_d
%             figure;plot(theta_hat);
%
%             figure;
%             plot(x(dex0),y(dex0),'.r','markersize',25); hold on;
%             plot(x(dex1),y(dex1),'.b','markersize',25); hold on;
%             plot(x,y); hold on;

end
