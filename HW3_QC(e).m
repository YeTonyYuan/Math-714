%QC(e) Long running time!
clear;clc;
Errors=zeros(2,8);
for B=1:8
   n=32;
   u=FDM(B,32);
   u_approx=FDM(B,64);
   e=zeros(33,33);
    for j=1:n+1
        for k=1:n+1
        e(j,k)=u(j,k)-u_approx((j-1)*2+1,(k-1)*2+1);
        end
    end
    Errors(1,B)=norm(e(:),inf);
    
    n=64;
    u=Spectral3(B,64);
   u_approx=Spectral3(B,128);
   e=zeros(33,33);
  for j=1:n/2+1
        for k=1:n/2+1
        e(j,k)=u(j,k)-u_approx((j-1)*2+1,(k-1)*2+1);
        end
    end  
    Errors(2,B)=norm(e(:),inf);
    



end


h=1:1:8;
figure
plot(h,Errors(1,:),'o-', 'LineWidth', 2)
hold on; 
plot(h,Errors(2,:),'o-', 'LineWidth', 2)

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error at t=1', 'FontSize', 24);
xlabel('$B$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("Finite Difference", "Spectral",'FontSize', 24 );
lgd.Location = 'northwest';




       
%%
function [u]=FDM(B,N)

dx = 1/N;
dy = dx;
sigma = 1/sqrt(2); gamma = 1/sqrt(2); %Courant-Friedrich Stability Condition
dt = sigma*(dx);
t = 0:dt:1; x = 0:dx:1; y = 0:dy:1; 
u = zeros(length(x),length(y),length(t));

u(:,:,1)=0;
u(:,:,2) = transpose(sin(B*pi*x))*sin(B*pi*y)*dt;

%u(x,y,dt)


for n=2:length(t)-1
    for i=2:length(x)-1
        for j=2:length(y)-1
            u(i,j,n+1)= (sigma^2)*(u(i+1,j,n)-2*u(i,j,n)+u(i-1,j,n))...
                +(gamma^2)*(u(i,j+1,n)-2*u(i,j,n)+u(i,j-1,n)) + 2*u(i,j,n) - u(i,j,n-1);
        end
    end
    
    for j=2:length(y)-1
        u(1,j,n+1)= (sigma^2)*(u(3,j,n)-2*u(2,j,n)+u(1,j,n))...
                +(gamma^2)*(u(1,j+1,n)-2*u(1,j,n)+u(1,j-1,n)) + 2*u(1,j,n) - u(1,j,n-1);
        u(end,j,n+1)= (sigma^2)*(u(end-2,j,n)-2*u(end-1,j,n)+u(end,j,n))...
                +(gamma^2)*(u(end,j+1,n)-2*u(end,j,n)+u(end,j-1,n)) + 2*u(end,j,n) - u(end,j,n-1);
    end
    
    for i=2:length(x)-1
         u(i,1,n+1)= (sigma^2)*(u(i+1,1,n)-2*u(i,1,n)+u(i-1,1,n))...
                +(gamma^2)*(u(i,3,n)-2*u(i,2,n)+u(i,1,n)) + 2*u(i,1,n) - u(i,1,n-1);
         u(i,end,n+1)= (sigma^2)*(u(i+1,end,n)-2*u(i,end,n)+u(i-1,end,n))...
                +(gamma^2)*(u(i,end-2,n)-2*u(i,end-1,n)+u(i,end,n)) + 2*u(i,end,n) - u(i,end,n-1);
    end
    
    u(1,1,n+1)= (sigma^2)*(u(3,1,n)-2*u(2,1,n)+u(1,1,n))...
                +(gamma^2)*(u(1,3,n)-2*u(1,2,n)+u(1,1,n)) + 2*u(1,1,n) - u(1,1,n-1);
            
    u(end,1,n+1)= (sigma^2)*(u(end-2,1,n)-2*u(end-1,1,n)+u(end,1,n))...
                +(gamma^2)*(u(end,3,n)-2*u(end,2,n)+u(end,1,n)) + 2*u(end,1,n) - u(end,1,n-1);
            
    u(1,end,n+1)= (sigma^2)*(u(3,end,n)-2*u(2,end,n)+u(1,end,n))...
                +(gamma^2)*(u(1,end-2,n)-2*u(1,end-1,n)+u(1,end,n)) + 2*u(1,end,n) - u(1,end,n-1);        
            
    u(end,end,n+1)= (sigma^2)*(u(end-2,end,n)-2*u(end-1,end,n)+u(end,end,n))...
                +(gamma^2)*(u(end,end-2,n)-2*u(end,end-1,n)+u(end,end,n)) + 2*u(end,end,n) - u(end,end,n-1);        
            
end
u=u(:,:,length(t));

end

%%
% p20.m - 2nd-order wave eq. in 2D via FFT (compare p19.m)
% Grid and initial data:
function u=Spectral3(B,N)
%N = 24; 
N=128;
x = cos(pi*(0:N)/(N)); y = x';
dt = 10/N^2;
[xx,yy] = meshgrid(x,y);
plotgap = round((1/3)/dt); dt = (1/3)/plotgap;
vv = 0;%Initial condition
vvold = vv;
vv=vvold+dt*sin(B*pi*xx)*sin(B*pi*yy)...
+dt^3/6*(-2*B^2*pi^2)*sin(B*pi*xx)*sin(B*pi*yy); %at t=dt

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