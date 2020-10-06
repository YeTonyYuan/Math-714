%I used some codes that were posted to Canvas.

%Question C(b)
clear;clc;

N=8; %n=2^8 for the finest mesh
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
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';






%% Question C(g)
%Error compared with true solution for BC f^ using G-S

clear;clc;

N=8; %n=2^8
Errors=zeros(N,1);
for i=1:N
    n=2^i;
    h=1/n;
    u=GS_2(n);
    x = 0:h:1;
    y = 0:h:1;

    [X Y] = meshgrid(x,y) ;
    sol = cos(2*pi.*Y(:,2:end-1)).*(exp(-2*pi.*X(:,2:end-1))*exp(4*pi)...
        - exp(2*pi.*X(:,2:end-1)))/( exp(4*pi)-1 ) ;
    e_2=u(:,2:end-1)-sol;
    Errors(i)=norm(e_2(:),inf);
end
h=2.^[-1:-1:-N];
figure
loglog(h,Errors)
hold on; 
loglog(h, h.^(-1), 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error compared with true solution for BC $\hat{f}$ using G-S', 'Interpreter','latex','FontSize', 24);
xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^{-1})$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';

%% Question D
clear;clc;
N=8;
u_approx=MG(2^N);
Errors=zeros(N-1,1);
for i=1:N-1
    n=2^i;
    u=MG(n);
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

title('Error for BC f using multigrid', 'FontSize', 24);
xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';


%% Question E
% Error compared with finer grids for BC f^ using G-S
clear;clc;

N=8; %n=2^8 as reference solution
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
loglog(h,Errors)
hold on; 
loglog(h, h, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error compared with reference solution for BC $\hat{f}$ using G-S','Interpreter','latex', 'FontSize', 24);
xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';
        

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

for iter=1:n^2*log(n)  
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
function [u] = MG(n)

h = 1 / n;
u = zeros(n+1, n+1);

for q=1:n
for i=1:n+1
    u(i,1) = cos(2*pi*(i-1)*h);
end

for iter=1:10
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
u_h=u(:,2:end-1);

%Forming Delta:
    dx=h; dy=h;
    x = 0:h:1;
    y = 0:h:1;
    
    % Building the 1D Dirichelet minus Laplacian matrix
    e = ones(n-1,1);
    Asp = spdiags([-e 2*e -e], -1:1, n-1, n-1);
  
    % Building the 1D Neumann minus Laplacian matrix
    e = ones(n+1,1);
    Bsp = spdiags([-e 2*e -e], -1:1, n+1, n+1);
    Bsp(1,1:2) = [1 -1 ]*dy;
    Bsp(end,end-1:end) = [ -1 1 ]*dy;
    
    % creating the identities (here carefull with the 
    % boundaries)
    I_A = speye(n+1,n+1);   I_A(1,1) = 0;   I_A(end,end) = 0;
    I_B = speye(n-1,n-1);
    
    % assembling the 2D minus Laplacian
    Delta = kron(Asp/dx^2,I_A) + kron(I_B,Bsp/dy^2);
    
    % writing the source term
    f = cos(2*pi*y);    f(1) = 0;   f(end) = 0;
    e_1 = zeros(n-1,1);   e_1(1) = 1;
    F = kron(e_1, f).';
    f = F(:)/dx^2;
    
    r_h=f- Delta*u_h(:);
    r_h=reshape(r_h,n+1,n-1);
    r_H=zeros(n/2+1,n/2-1);
    for i=1:n/2+1
        for j=1:n/2-1
        r_H(i,j)=r_h(2*i-1,2*j);
        end
    end
    r_H=r_H(:);
    
    % Forming Delta_H
    dx_H=2*h;dy_H=2*h;
    x_H = 0:2*h:1;
    y_H = 0:2*h:1;
    
    % Building the 1D Dirichelet minus Laplacian matrix
    e = ones(n/2-1,1);
    Asp = spdiags([-e 2*e -e], -1:1, n/2-1, n/2-1);
  
    % Building the 1D Neumann minus Laplacian matrix
    e = ones(n/2+1,1);
    Bsp = spdiags([-e 2*e -e], -1:1, n/2+1, n/2+1);
    Bsp(1,1:2) = [1 -1 ]*dy_H;
    Bsp(end,end-1:end) = [ -1 1 ]*dy_H;
    
    % creating the identities (here carefull with the 
    % boundaries)
    I_A = speye(n/2+1,n/2+1);   I_A(1,1) = 0;   I_A(end,end) = 0;
    I_B = speye(n/2-1,n/2-1);
    
    % assembling the 2D minus Laplacian
    Delta_H = kron(Asp/dx_H^2,I_A) + kron(I_B,Bsp/dy_H^2);
    e_H=Delta_H\ r_H;
    e_H=reshape(e_H,n/2+1,n/2-1);
    
    for i=1:n/2+1
        for j=1:n/2-1
            u(2*i-1,2*j+1)= u(2*i-1,2*j+1)+e_H(i,j);
        end
    end
end    

end






