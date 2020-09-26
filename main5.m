close all;clear;clc;
warning off;

%% plane quiver for z_angle = 0

%% Parameters
N=3; k=-0.75; eta=0.1;
nn = 1;nx = 21; ny = 21;nz = 5;
z_alist = linspace(0,2*pi,nz);
%% Calculation
for zz = 1:nz
    w_r = linspace(-1,1,nx); % real of w
    w_i = linspace(-1,1,nx); %  imag of w
    z_a = z_alist(zz); % angle of zeta
    z_l = 1; % length of zeta
    rec = zeros(nx,ny,4,nn); % record of w,zeta,Dw,Dzeta
    for n = 1:nn
        odeFunc=@(y)func(y,N,n-1,k,eta);
        for ii =1:nx
            for jj =1:ny
                y = [w_r(ii) + 1i*w_i(jj); z_l * exp(1i*z_a)];
                if abs(y(1)) >= 1
                    continue;
                end
                rec(ii,jj,:,n) = [y;odeFunc(y)];
            end
        end
    end
    %% plot figures
    for n = 1:nn
        w = rec(:,:,1,n); z = rec(:,:,2,n);
        Dw = rec(:,:,3,n); Dz = rec(:,:,4,n);
        fig = figure;
        set(fig, 'position', get(0,'ScreenSize')); % Fullscreen
        X= real(w); Y=imag(w);
        U = real(Dw); V = imag(Dw);
        % U = real(Dw)./abs(Dw); V = imag(Dw)./abs(Dw);
        quiver(X,Y,U,V);
        hold on;
        set(gca,'FontSize',30);
        xlabel('Re[\omega]');
        ylabel('Im[\omega]');
        axis equal;
        grid on;
        xlim([-1 1]);
        ylim([-1 1]);
        title(['n = ',num2str(n-1),' ¡Ï\zeta = ',num2str(round(z_a/pi,2)),'\pi'],'FontSize',36);
        plotcircle();
    end
end

%% main ode function
function dy = func(y,N,n,k,eta)
persistent an tmpn;
if isempty(an) || isempty(tmpn) || tmpn ~= n
    an = 2*pi / integral(@(t)((1-cos(t)).^n),0,2*pi);
    tmpn = n;
end
I = an/N*sum((1-cos(y)).^n);
w = y(1); z = y(2);
dy = zeros(size(y));
dy(1) = -0.5*1i*(1-abs(w)^2)*conj(z)*(eta+k*I-1);
dy(2) = 1i*(1+eta+k*I)*z-0.5*1i*(conj(w)+w*z^2)*(eta+k*I-1);
end

function plotcircle()
t = linspace(0,2*pi,360);
x = cos(t);
y = sin(t);
plot(x,y,'color','k','linewidth',5);
end
