clear;
Spectral(64)

%%
%QC(c)
clear;clc;
N=2^5; %n=2^9 for the finest mesh
u_approx=Spectral(N);
Errors=zeros(4,1);
for i=1:4
    n=2^i;
    u=Spectral(n);
    e=zeros(n+1);
    for j=1:n+1
        for k=1:n+1
        e(j,k)=u(j,k)-u_approx((j-1)*2^(5-i)+1,(k-1)*2^(5-i)+1);
        end
    end
    Errors(i)=norm(e(:),inf);
end
h=2.^[-1:-1:-4];
figure
loglog(h,Errors,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^4, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error at t=1', 'FontSize', 24);
xlabel('$\Delta x$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';



%%
function u=Spectral(N) %Initial conditions from HW2
% p20.m - 2nd-order wave eq. in 2D via FFT (compare p19.m)
% Grid and initial data:

x = cos(pi*(0:N)/(2*N)); y = x'; %x between [0,1]!!!!!!!
dt = 3/N^2;
[xx,yy] = meshgrid(x,y);
plotgap = round((1/3)/dt); dt = (1/3)/plotgap;

vv = 0;%Initial condition
vvold = vv;
vv=vvold+dt*exp(-400*((xx-.5).^2 + (yy-.5).^2))...
+dt^3/6*exp(-400*((xx-.5).^2 + (yy-.5).^2))*(800^2*(xx-0.5)^2-800+800^2*(yy-0.5)^2-800); %at t=dt



% Time-stepping by leap frog formula:

for n = 0:3*plotgap
t = n*dt;

uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
D2u= zeros(N+1,N+1); D2uxx=zeros(N+1,N+1);D2uyy=zeros(N+1,N+1);
D4u= zeros(N+1,N+1);
ii = 1:N+1;
for i = 1:N+1 % 2nd derivs wrt x in each row
v = vv(i,:); 
uxx(i,:)=chebfft(chebfft(v));
end
for j = 1:N+1 % 2nd derivs wrt x in each row
v = vv(:,j); 
uyy(:,j)=chebfft(chebfft(v));
end
D2u=uxx+uyy;

for i = 1:N+1 % 2nd derivs wrt x in each row
v = D2u(i,:); 
D2uxx(i,:)=chebfft(chebfft(v));
end
for j = 1:N+1 % 2nd derivs wrt x in each row
v = D2u(:,j); 
D2uyy(:,j)=chebfft(chebfft(v));
end
D4u=D2uxx+D2uyy;

vvnew = 2*vv - vvold + dt^2*D2u+ dt^4/12*D4u;
vvold = vv; vv = vvnew;
end

u=vv;

end












%%
function w = chebfft(v)
N = length(v)-1; if N==0, w=0; return, end
x = cos((0:N)'*pi/(2*N)); %x between [0,1]
ii = 0:N-1;
v = v(:); V = [v; flipud(v(2:N))]; % transform x -> theta
U = real(fft(V));
W = real(ifft(1i*[ii 0 1-N:-1]'.*U));
w = zeros(N+1,1);
w(2:N) = -W(2:N)./sqrt(1-x(2:N).^2); % transform theta -> x
w(1) = sum(ii'.^2.*U(ii+1))/N + .5*N*U(N+1);
w(N+1) = sum((-1).^(ii+1)'.*ii'.^2.*U(ii+1))/N + ...
.5*(-1)^(N+1)*N*U(N+1);
end
