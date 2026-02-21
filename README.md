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
%then for the script
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
