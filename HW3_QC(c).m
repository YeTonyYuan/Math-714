%QC(c) Long running time!
clear;clc;
n1=8;
N=2^n1; %n=2^9 for the finest mesh
u_approx=Spectral2(N);
u_approx=u_approx(1:N/2+1,1:N/2+1);
Errors=zeros(n1-2,1);
for i=2:n1-1
    n=2^i;
    u=Spectral2(n);
    u=u(1:n/2+1,1:n/2+1);
    e=zeros(n+1);
    for j=1:n/2+1
        for k=1:n/2+1
        e(j,k)=u(j,k)-u_approx((j-1)*2^(n1-i)+1,(k-1)*2^(n1-i)+1);
        end
    end
    Errors(i-1)=norm(e(:),inf);
end
h=2.^[-2:-1:-(n1-1)];
figure
loglog(h,Errors,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^4, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error at t=1', 'FontSize', 24);
xlabel('$1/N$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(N^{-4})$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';



%%
clear;clc;
% p20.m - 2nd-order wave eq. in 2D via FFT (compare p19.m)
% Grid and initial data:
function u=Spectral2(N)
%N = 24; 
x = cos(pi*(0:N)/(N)); y = x';
dt = 10/N^2;
[xx,yy] = meshgrid(x,y);
plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
vv = 0;%Initial condition
vvold = vv;
vv=vvold+dt*exp(-400*((xx-.5).^2 + (yy-.5).^2))...
+dt^3/6*exp(-400*((xx-.5).^2 + (yy-.5).^2))*(800^2*(xx-0.5)^2-800+800^2*(yy-0.5)^2-800); %at t=dt

% Time-stepping by leap frog formula:
[ay,ax] = meshgrid([.56 .06],[.1 .55]); clf
for n = 0:3*plotgap
t = n*dt;
if rem(n+.5,plotgap)<1 % plots at multiples of t=1/3
i = n/plotgap+1;
subplot('position',[ax(i) ay(i) .36 .36])
[xxx,yyy] = meshgrid(-1:1/16:1,-1:1/16:1);
vvv = interp2(xx,yy,vv,xxx,yyy,'cubic');
mesh(xxx,yyy,vvv), axis([-1 1 -1 1 -0.15 1])
colormap([0 0 0]), title(['t = ' num2str(t)]), drawnow
end
uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
D2u= zeros(N+1,N+1); D2uxx=zeros(N+1,N+1);D2uyy=zeros(N+1,N+1);
D4u= zeros(N+1,N+1);
ii = 2:N;
for i = 2:N % 2nd derivs wrt x in each row
v = vv(i,:); V = [v fliplr(v(ii))];
U = real(fft(V));
W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U)); % diff wrt theta
W2 = real(ifft(-[0:N 1-N:-1].^2.*U)); % diff^2 wrt theta
uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).* ...
W1(ii)./(1-x(ii).^2).^(3/2);
end
for j = 2:N % 2nd derivs wrt y in each column
v = vv(:,j); V = [v; flipud(v(ii))];
U = real(fft(V));
W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));% diff wrt theta
W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U)); % diff^2 wrt theta
uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).* ...
W1(ii)./(1-y(ii).^2).^(3/2);
end
D2u=uxx+uyy;
for i = 2:N % 2nd derivs wrt x in each row
v = D2u(i,:); V = [v fliplr(v(ii))];
U = real(fft(V));
W1 = real(ifft(1i*[0:N-1 0 1-N:-1].*U)); % diff wrt theta
W2 = real(ifft(-[0:N 1-N:-1].^2.*U)); % diff^2 wrt theta
D2uxx(i,ii) = W2(ii)./(1-x(ii).^2) - x(ii).* ...
W1(ii)./(1-x(ii).^2).^(3/2);
end
for j = 2:N % 2nd derivs wrt y in each column
v = D2u(:,j); V = [v; flipud(v(ii))];
U = real(fft(V));
W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));% diff wrt theta
W2 = real(ifft(-[0:N 1-N:-1]'.^2.*U)); % diff^2 wrt theta
D2uyy(ii,j) = W2(ii)./(1-y(ii).^2) - y(ii).* ...
W1(ii)./(1-y(ii).^2).^(3/2);
end
D4u=D2uxx+D2uyy;
vvnew = 2*vv - vvold + dt^2*(D2u)+dt^4/12*D4u;
vvold = vv; vv = vvnew;
end
u=vvnew;
end