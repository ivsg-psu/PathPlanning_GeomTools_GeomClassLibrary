close all

while(1)

    for ith_frame = 1:48
        drawframe(ith_frame);
        pause(0.01);
    end
end

% https://www.mathworks.com/matlabcentral/communitycontests/contests/6/entries/14362

function drawframe(f)
% This function draws a frame in an animation
% visualising the effect of zooming in on a part of a Julia set,
% 
% The code here is a remix of Julia sets:
%     https://uk.mathworks.com/matlabcentral/communitycontests/
%             contests/6/entries/14057
% For further explanation, see the comments there.
% Frame geometry.
nFrame = 48;
x0 = -0.030006;
y0 = 0.019988;
% Define frame-dependent range.
dxy = 1.5 * 10^(-0.05*(f-1));
nxy = 300;
% Complex constant.
r = 1;
theta = 13 * pi / 24;
c = r * cos(theta) + 1i * r * sin(theta);
% Iiterations and divergence criterion.
itMax = 50;
zDiverged = 1.5;
% Grid of points in the complex plane.
x = linspace(x0-dxy,x0+dxy, nxy);
min(x)
max(x)
y = linspace(y0-dxy, y0+dxy, nxy);
[X,Y] = meshgrid(x,y);
Z = complex(X,Y);
C = ones(size(X))*c;
V = zeros(size(C));
for k = 1:itMax
    Z = Z.^2 + C;
    V = V + (abs(Z)<zDiverged);
end
% Visualise the Julia set.
axes(Position=[0 0 1 1])
imagesc(V);
colormap(parula);
axis off
end