clear
clc
% length = 6m       semi major length = 3m
% height  = 4       semi minor length = 2m
nx=41;ny=41;                                            % Number of points in X and Y
minerr=1e-06;                                           % Minimum error

%% Physical domain points of X location
for n=1:nx
    
for k=2:ny-1
X(1,k) = 12 - (k-1)*12/(nx-1);
X(nx,k) = 12 - (k-1)*12/(nx-1);
end

X(1:5,n) = 12 -(n-1)*12/(nx-1);   X(37:41,n) = 12 - (n-1)*12/(nx-1);
X(1:5,1) = 12;    X(37:41,1) = 12;
X(17:25,n) = 18 + (n-1)*27/(nx-1);
X(1:5,ny) = 0;    X(37:41,ny) = 0;
X(17:25,ny) = 45;

for k=6:16
    X(k,1) = X(k-1,1) + 0.5;
    X(k,n) = X(k-1,n) + n*3.75/(ny-1);
    X(k,ny) = X(k-1,ny) + 3.75;
end
for k=36:-1:26
    X(k,1) = X(k-1,1) - 0.5;
    X(k,n) = X(k+1,n) + n*3.75/(ny-1);
    X(k,ny) = X(k-1,ny) - 3.75;
end
%% Physical domain points of Y location
Y(1,n) = 10;      Y(21,n) = 10;     Y(nx,n) = 10;
Y(1,1) = 10;      Y(21,1) = 10;     Y(nx,1) = 10;
Y(1,ny) = 10;     Y(21,ny) = 10;    Y(nx,ny) = 10;

Y(5:17,n) = 12 + (n-1)*(8)/(ny-1);       Y(25:37,n) = 8 - (n-1)*(8)/(ny-1);
Y(5:17,ny) = 20;       Y(25:37,ny) = 0;
for k=1:4
    Y(k,1) = 10 + (k-1)*0.5;
    Y(k,n) = 10 + (k-1)*(Y(5,n) - Y(37,n))/8;
    Y(k,ny) = 10 + (k-1)*2.5;
end
for k=18:24
    Y(k,1) = 12 - (k-17)*0.5;
    Y(k,n) = Y(17,n) - (k-17)*(Y(17,n) - Y(25,n))/8;
    Y(k,ny) = 20 - (k-17)*2.5;
end

for k=38:41
    Y(k,1) = 8 + (k-37)*0.5;
    Y(k,n) = 10 + (k-ny)*(Y(5,n) - Y(37,n))/8;
    Y(k,ny) = 0 + (k-37)*2.5;
end
Y(1,1:ny) = 10;
Y(1,1:ny) = 10;

end

% N+1 iteration level Boundary points

newX = X;
newY = Y;

for t=1:2000
%%  computing the parameters Alpha Beta and Gamma

for j=2:ny-1
    for i=2:nx-1
        
alpha(i,j)=(1/4)*((X(i,j+1)-X(i,j-1))^2 + (Y(i,j+1) - Y(i,j-1))^2);
alpha(1,j)=(1/4)*((X(1,j+1)-X(1,j-1))^2 + (Y(1,j+1) - Y(1,j-1))^2);
 
beta(i,j)=(1/4)*((X(i+1,j)-X(i-1,j))*(X(i,j+1)-X(i,j-1))+...
(Y(i+1,j) - Y(i-1,j))*(Y(i,j+1) - Y(i,j-1)));
beta(1,j)=(1/4)*((X(2,j)-X(nx-1,j))*(X(1,j+1)-X(1,j-1))+...
(Y(2,j) - Y(nx-1,j))*(Y(1,j+1) - Y(1,j-1)));
 
gamma(i,j)=(1/4)*((X(i+1,j)-X(i-1,j))^2 + (Y(i+1,j) - Y(i-1,j))^2);
gamma(1,j)=(1/4)*((X(2,j)-X(nx-1,j))^2 + (Y(2,j) - Y(nx-1,j))^2);
 
%   solving the poisson equation by iteration method
 
newX(i,j)=((-0.5)/(alpha(i,j)+gamma(i,j)))*(0.5*beta(i,j)*(X(i+1,j+1)-X(i-1,j+1)-X(i+1,j-1) + X(i-1,j-1))...
-alpha(i,j)*(X(i+1,j)+X(i-1,j))-gamma(i,j)*(X(i,j+1)+X(i,j-1)));
 
newY(i,j)=((-0.5)/(alpha(i,j)+gamma(i,j)))*(0.5*beta(i,j)*(Y(i+1,j+1)-Y(i-1,j+1)-Y(i+1,j-1) + Y(i-1,j-1))...
-alpha(i,j)*(Y(i+1,j)+Y(i-1,j))-gamma(i,j)*(Y(i,j+1)+Y(i,j-1)));
 
newX(1,j)=((-0.5)/(alpha(1,j)+gamma(1,j)))*(0.5*beta(1,j)*(X(2,j+1)-X(nx-1,j+1)-X(2,j-1) + X(nx-1,j-1))...
-alpha(1,j)*(X(2,j)+X(nx-1,j))-gamma(1,j)*(X(1,j+1)+X(1,j-1)));
 
newY(1,j)=((-0.5)/(alpha(1,j)+gamma(1,j)))*(0.5*beta(1,j)*(Y(2,j+1)-Y(nx-1,j+1)-Y(2,j-1) + Y(nx-1,j-1))...
-alpha(1,j)*(Y(2,j)+Y(nx-1,j))-gamma(1,j)*(Y(1,j+1)+Y(1,j-1)));
    end
end

newX(nx,:) = newX(1,:);
newY(nx,:) = newY(1,:);

Er1 = 0; Er2 = 0;
for j=2:ny-1
    for i=2:nx-1
        Er1 = Er1 + (abs(newX(i,j)-X(i,j)));
        Er2 = Er2 + (abs(newY(i,j)-Y(i,j)));
    end
end
        Err = (Er1 + Er2)/nx;
    X=newX;
    Y=newY;
    if Err<minerr
        break
    end

end

%% To plot Mesh in Physical domain
figure(1)
hold on
axis equal
for m=1:nx
plot(X(m,:),Y(m,:),'b');
end
for m=1:ny
plot(X(:,m),Y(:,m),'Color',[0 0 0]);
end
% print(gcf,'MESH.jpg','-dpng','-r300');
