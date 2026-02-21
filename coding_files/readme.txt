function [A,E,R_2] = fit_linear(Z,Y)
n = size(Z,1);
Z = [ones(n,1) Z];
T = transpose(Z);
A = (T * Z)\(T*Y);
Y_hat = Z * A;
E = Y - Y_hat;
Sr = transpose(E) * E;
St = sum((Y - mean(Y)).^2);
R_2 = 1 - (Sr / St);
end

data = readmatrix('yacht_hydrodynamics.csv');

Z = data(:,1:6);   
Y = data(:,7);     

[A,E,R2] = fit_linear(Z,Y);

disp('Fit parameters:')
disp(A)

disp("R^2 = " + R2_B)


x_cube = ones(1,6);

x_cube = [1 x_cube];

y_cube = x_cube * A;

disp("Residuary resistance for cube boat = " + y_cube)

function [A] = fit_nonlinear(x,y,model,seeds)
p = seeds;                 
h = 1e-6;
max = 100;
for i = 1:max

    Y_hat = model(x,p);
    E = y - Y_hat;

    z = zeros(length(x),length(p));

    for j = 1:length(p)
        K = zeros(length(p),1);
        K(j) = h;
        Y_diff = model(x,p+K);
        z(:,j) = (Y_diff - Y_hat)/h;
    end

    T = transpose(z);
    D = ((T*z) + h*eye(size(T*z))) \ (T*E);
    p = p + 0.1*D;
end
A = p;
end


data = readmatrix('rheo_data.csv');

X = data(:,1);     
Y = data(:,2);     

% code for model 1
Z = X;    
[A_1,E,R2] = fit_linear(Z,Y);

T = A_1(1);
N = A_1(2);

%code for model 2

Hmodel = @(x,p) p(1) + p(2)*x.^p(3);

seeds_1 = [1; 1; 1];

A_2 = fit_nonlinear(X,Y,Hmodel,seeds_1);

%code for model 3

Pmodel = @(x,p) p(1) + p(2)*x + p(3)*x.^p(4);


seeds_2 = [100; 1; 1; 1];

A_3 = fit_nonlinear(X,Y,Pmodel,seeds_2);


x = linspace(min(X), max(X), 200);

y_1 = T + N*x;
y_2 = Hmodel(x, A_2);
y_3 = Pmodel(x, A_3);

%graph

semilogx(X,Y,'o','MarkerSize',10,'LineWidth',1.5)
hold on
semilogx(x,y_1,'LineWidth',2)
semilogx(x,y_2,'LineWidth',2)
semilogx(x,y_3,'LineWidth',2)

xlabel('Shear Rate [1/s]','FontSize',14)
ylabel('Shear Stress [Pa]','FontSize',14)
legend('Data','Bingham','Herschel-Bulkley','HB Plus','Location','best')
set(gca,'FontSize',20)
