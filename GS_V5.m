%Question C(b)
clear;clc;

N=9;
u_approx=GS(2^N);
Errors=zeros(N-1,1);
for i=1:N-1
    n=2^i;
    u=GS(n);
    e=zeros(n+1);
    for j=1:n+1
        for k=1:n+1
        e(j,k)=u(j,k)-u_approx((j-1)*2^(N-i)+1,(k-1)*2^(N-i)+1);
        end
    end
    Errors(i)=norm(e(:),inf);
end
h=2.^[-1:-1:-(N-1)];
figure
loglog(h,Errors,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error for BC f using G-S', 'FontSize', 24);
xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^2$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';


%% Question C(g)
clear;clc;

N=8;
u_approx=GS_2(2^N);
Errors=zeros(N-1,1);
for i=1:N-1
    n=2^i;
    u=GS_2(n);
    e=zeros(n+1);
    for j=1:n+1
        for k=1:n+1
        e(j,k)=u(j,k)-u_approx((j-1)*2^(N-i)+1,(k-1)*2^(N-i)+1);
        end
    end
    Errors(i)=norm(e(:),inf);
end
h=2.^[-1:-1:-(N-1)];
figure
loglog(h,Errors,h,h )






%%

clear;clc;

N=8;
Errors=zeros(N-1,1);
for i=1:N
    n=2^i;
    h=1/n;
    u=GS_2(n);
    e_2=zeros(n+1,n);
    x = 0:h:1;
    y = 0:h:1;

    [X Y] = meshgrid(x,y) ;
    sol = cos(2*pi.*Y(:,2:end-1)).*(exp(-2*pi.*X(:,2:end-1))*exp(4*pi)...
        - exp(2*pi.*X(:,2:end-1)))/( exp(4*pi)-1 ) ;
    e_2=u(:,2:end-1)-sol;
    Errors(i-1)=norm(e_2(:),inf);
end
h=2.^[-1:-1:-N];
figure
loglog(h,Errors,h,h.^2 )
 
        

%%
function [u] = GS(n) %Gauss-Seidel iterations using BC f.
% ...

h = 1 / n;
u = zeros(n+1, n+1);
for i=1:n+1
    u(i,1) = cos(2*pi*(i-1)*h);
end

for iter=1:n^2
    for i=2:n        
    u(1,i) = (1/4)*(2*u(2,i) +u(1,i-1)...
  +u(1,i+1));
    u(end,i) = (1/4)*(2*u(end-1,i) +u(end,i-1)...
  +u(end,i+1));
    end    
    for x=2:n
        for y=2:n
            u(x,y) = (1/4)*(  u(x-1,y) + u(x+1,y) ...
                                   + u(x,y-1) + u(x,y+1) );
        end
    end
   
end

end

%%

function [u] = GS_2(n) %Gauss-Seidel iterations using BC f^.

h = 1 / n;
u = zeros(n+1, n+1);
for i=1:n+1
    u(i,1) = cos(2*pi*(i-1)*h);
    if u(i,1) >=0
        u(i,1)=1;
    else
        u(i,1)=-1;
    end
end

for iter=1:n^3   %# of iterations n^2 log n
    for i=2:n        
    u(1,i) = (1/4)*(2*u(2,i) +u(1,i-1)...
  +u(1,i+1));
    u(end,i) = (1/4)*(2*u(end-1,i) +u(end,i-1)...
  +u(end,i+1));
    end    
    for x=2:n
        for y=2:n
            u(x,y) = (1/4)*(  u(x-1,y) + u(x+1,y) ...
                                   + u(x,y-1) + u(x,y+1) );
        end
    end
   
end

end






