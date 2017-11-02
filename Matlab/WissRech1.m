%% Wissenschatfliches Rechnen 1 - Zumbusch

%% Aufgabenblatt 1

% DGL:   - (d^2u/dx^2 + d^2u/dy^2) =  x(1-x)+y(1-y);



ngrid=100;


x=linspace(0,1,ngrid+2);
x=x(2:end-1);
y=x;

[X,Y]=meshgrid(x,y);

boundary=0;

rho=zeros(ngrid,ngrid);
rho(2:end-1,2:end-1)=-(X(2:end-1,2:end-1).*(1-X(2:end-1,2:end-1))+Y(2:end-1,2:end-1).*(1-Y(2:end-1,2:end-1)));
rho_vec=reshape(rho(2:end-1,2:end-1)',[1,(ngrid-2)^2]);

blockdiag=ones(1,(ngrid-2)^2-1);
blockdiag(mod([1:(ngrid-2)^2-1],ngrid-2)==0)=0;
laplace=diag(-4*ones(1,(ngrid-2)^2),0)+diag(blockdiag,1)+diag(blockdiag,-1)+diag(ones(1,(ngrid-2)^2-(ngrid-2)),ngrid-2)+diag(ones(1,(ngrid-2)^2-(ngrid-2)),-(ngrid-2));

clear('blockdiag');

rho_vec=rho_vec/laplace*(x(2)-x(1));

sol=zeros(ngrid,ngrid);
sol(2:end-1,2:end-1)=reshape(rho_vec,[ngrid-2,ngrid-2])';

figure(2);
imagesc(sol);
figure(1)
imagesc(rho)





























%