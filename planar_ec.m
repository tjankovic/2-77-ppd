clear
%Script to explore sensitivities and configuration stability for planar
%exact constraint toy
%Shien Yang Lee 2017, GPLv3

%the current version assumes a square planar object with side lengths d
%coordinate system: right handed cartesian centered at center of square

d = 10; %object side length
contact_pts = [4 -5;... %contact points between pegs and object
              -3 -5;...
               5  0];

theta = deg2rad(0); %angle of weight vector in deg, CW from +x

W = 10; %weight of planar object

syms r1x r2x r3x r1y r2y r3y
R = [r1x r1y;... %reaction forces at each contact
     r2x r2y;...
     r3x r3y];
R(abs(contact_pts)~=5) = zeros(1,3);
var_list = reshape(transpose(R),1,6); %eliminate non-normal components of forces (frictionless assumption)
var_list = var_list(var_list~=0);

%TO-DO: check contact points for validity (must lie on perimeter)
for i=1:size(contact_pts,1)
    pt = contact_pts(i,:);
    assert(any(abs(pt)==d./2))
end

%calculate moment contributions of reaction forces
R = [R [0; 0; 0]];
contact_pts = padarray(contact_pts, [0 1],'post'); %pad to 3d for cross product function
M = cross(contact_pts,R,2); %cross product along 2nd dim gives moment contributions in vector form
m = M(:,3); %convert to row vector since problem is 2D

%static equilibrium equations
eqn_m = m(1) + m(2) + m(3) == 0;
eqn_x = R(1,1) + R(2,1) + R(3,1) + W*cos(theta) == 0;
eqn_y = R(1,2) + R(2,2) + R(3,2) - W*sin(theta) == 0;

[A,B] = equationsToMatrix([eqn_m, eqn_x, eqn_y], var_list);
sol = linsolve(A,B);
sol = eval(sol)
