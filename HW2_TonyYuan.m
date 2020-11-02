%% QB
clear;clc;
error=[];
for  N=10:10:120
x=0:1/N:1;
v=exp(-400*(x-0.5).^2);
xq=0:1/N^2:1;
fq=exp(-400*(xq-0.5).^2);
vq=interp1(x,v,xq);
e=norm(fq-vq,inf);
error=[error,e]
end


%% QC
clear;clc;
N=10; %n=2^11 for the finest mesh
u_approx=FDM(2^N);
Errors=zeros(4,1);
for i=6:N-1
    n=2^i;
    u=FDM(n);
    e=zeros(n+1);
    for j=1:n+1
        for k=1:n+1
        e(j,k)=u(j,k)-u_approx((j-1)*2^(N-i)+1,(k-1)*2^(N-i)+1);
        end
    end
    Errors(i-5)=norm(e(:),inf);
end
h=2.^[-6:-1:-(N-1)];
figure
loglog(h,Errors,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^2, 'LineStyle', '-')

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
function [u]=FDM(N)
dx = 1/N;
dy = dx;
sigma = 1/sqrt(2); gamma = 1/sqrt(2); %Courant-Friedrich Stability Condition
dt = sigma*(dx);
t = 0:dt:1; x = 0:dx:1; y = 0:dy:1; 
u = zeros(length(x),length(y),length(t));

u(:,:,1)=0;
u(:,:,2) = transpose(exp(-400*(x-0.5).^2))*exp(-400*(y-0.5).^2)*dt;

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





