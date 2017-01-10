Point(1) = {0,0,0,0.1};
Point(2) = {1,0,0,0.1};
Point(3) = {1,1,0,0.1};
Point(4) = {0,1,0,0.1};
Point(5) = {0.5, 0.5, 0, 0.01};
Point(6) = {0.5, 0.6, 0, 0.01};

Point(7) = {0.7, 0.6, 0, 0.01};
Point(8) = {0.75, 0.65, 0, 0.01};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Circle(6) = {6, 5, 6};
Circle(7) = {8, 7, 8};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Delete { Surface{1}; }

Line Loop(2) = {6};
Line Loop(3) = {7};

Plane Surface(1) = {1,2,3}; 

