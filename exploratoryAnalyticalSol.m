function exploratoryAnalyticalSol
global H L cao
H = 6; % Height in cm
L = 4; % Length in cm
cao = 0.1/(8.206*10^(-5)*300.15); % Concentration in mol/m^3


h = H / 25; % Increment in height
l = L / 25; % Increment in length

c = zeros(26,26); % Initializing concentration matrix

% Create concentration matrix
for i = 1:26
   for j = 1:26
       c(i,j) = findsol(l*(i-1),h*(j-1)); % Substitute each x and y into function
   end  
end

% Oops I got the height and length on the opposite sides
c = transpose(c);

% Plot as contour plot
    figure
    [C,h] = contour((0:l:L), (0:h:H), c);
    clabel(C,h)
    xlabel('x (cm)'), ylabel('y (cm)')
end

function s = findsol (x,y)

global H L cao

i = 0;
n = 1;

% Calculating sum part of the function
while n < 99
i = i+((-1)^n-1)/(n*pi*sinh(n*pi*H/L))*sin(n*pi*x/L)*sinh(n*pi*y/L);
n = n + 1;
end

% Multiply by constants before the sum part
s = -2*i*cao;

end