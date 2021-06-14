clear
clc
%%

numX = 101; numY = numX;

nx = numX ; ny = ceil(numY/2) ;
chord = 1 ;              % chord length of aerofoil
thickness = 0.3 ; 
err = 1e-06;             % error for gauss seidal for solving elliptic pde

X = zeros(nx,ny);
Y = zeros(nx,ny);

alpha = zeros(nx,ny);   beta = zeros(nx,ny);   gamma = zeros(nx,ny);

tmpX = linspace(0.0,3.0*chord,(ny-2*ceil(ny/4))+2) ;
tmpX = tmpX(2:length(tmpX)-1) ;

X(:,ny) = linspace( 2*chord, 3*chord,nx)	;
X(:,1)  = linspace( 1*chord, 0,nx)	;
X(1,:)  = linspace(1,2,ny) ;
X(nx,:) = [ zeros(1,ceil(ny/4)) tmpX (zeros(1,ceil(ny/4)) + 3*chord ) ] 	;

Y(:,ny) = zeros(1,nx)   ;
Y(:,1)  = zeros(1,nx)   ;

for i=1:ny
    Y(1,i) = 0.6*(0.2969*sqrt(X(1,i)-1) - 0.1260*(X(1,i)-1) - 0.3516*(X(1,i)-1)^2 + 0.2843*(X(1,i)-1)^3 - 0.1015*(X(1,i)-1)^4);
end

Y(nx,:) = [linspace(0,4*thickness,ceil(ny/4)) (zeros(1,(ny-2*ceil(ny/4)))+4.5*thickness) linspace(4.5*thickness,0,ceil(ny/4)) ] 	;

%%
% N+1 iteration level Boundary points

newX = X;
newY = Y;

Err = 1 ;   iter = 0 ;

while ( Err > err )
    %%  computing the parameters Alpha Beta and Gamma
    
    iter = iter + 1 ;
    
    for j=2:ny-1
        for i=2:nx-1
            
            alpha(i,j)=(1/4)*((X(i,j+1)-X(i,j-1))^2 + (Y(i,j+1) - Y(i,j-1))^2);
            
            
            beta(i,j)=(1/4)*((X(i+1,j)-X(i-1,j))*(X(i,j+1)-X(i,j-1))+...
                (Y(i+1,j) - Y(i-1,j))*(Y(i,j+1) - Y(i,j-1)));
            
            gamma(i,j)=(1/4)*((X(i+1,j)-X(i-1,j))^2 + (Y(i+1,j) - Y(i-1,j))^2);
            
            %   solving the poisson equation by iteration method
            
            newX(i,j)=((-0.5)/(alpha(i,j)+gamma(i,j)+14^-9))*(0.5*beta(i,j)*(newX(i+1,j+1)-newX(i-1,j+1)-newX(i+1,j-1) + newX(i-1,j-1))...
                -alpha(i,j)*(newX(i+1,j)+newX(i-1,j))-gamma(i,j)*(newX(i,j+1)+newX(i,j-1)));
            
            newY(i,j)=((-0.5)/(alpha(i,j)+gamma(i,j)+14^-9))*(0.5*beta(i,j)*(newY(i+1,j+1)-newY(i-1,j+1)-newY(i+1,j-1) + newY(i-1,j-1))...
                -alpha(i,j)*(newY(i+1,j)+newY(i-1,j))-gamma(i,j)*(newY(i,j+1)+newY(i,j-1)));
            
        end
    end
    
    
    Er1 = 0; Er2 = 0;
    for j=2:ny-1
        for i=2:nx-1
            Er1 = Er1 + (newX(i,j)-X(i,j))^2;
            Er2 = Er2 + (newY(i,j)-Y(i,j))^2;
        end
    end
    Err = (sqrt(Er1/(nx*ny)) + sqrt(Er2/(nx*ny)))/2;
    X=newX;
    Y=newY;
    
end

XX = [ X  flip(X(:,1:ny-1),2) ] ;
YY = [ Y -flip(Y(:,1:ny-1),2) ] ;

fprintf(' No iterations = %d for creating Mesh \n', iter)

%% uncomment below lines if you want to see grid points and mesh
% plot(XX,YY,'k*')  %% uncomment if you want to see mesh grid points


hold on             %% from here it plots mesh grid
axis equal
for m=1:numX
    plot(XX(m,:),YY(m,:),'b');
end
for m=1:numY
    plot(XX(:,m),YY(:,m),'Color',[0 0 0]);
end

% print(gcf,'NACA0012.jpg','-dpng','-r300');