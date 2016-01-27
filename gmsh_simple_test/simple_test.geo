cl_1 = 0.02;

Point(1) = {0, 0, 0, cl_1};
Point(2) = {0.1, 0, 0, cl_1};
Point(3) = {0.1, 0.1, 0, cl_1};
Point(4) = {0, 0.1, 0, cl_1};
Point(5) = {0, 0, 0.1, cl_1};
Point(6) = {0.1, 0, 0.1, cl_1};
Point(7) = {0.1, 0.1, 0.1, cl_1};
Point(8) = {0, 0.1, 0.1, cl_1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Line Loop(3) = {9, -8, -12, 4};
Plane Surface(3) = {3};
Line Loop(4) = {9, 5, -10, -1};
Plane Surface(4) = {4};
Line Loop(5) = {10, 6, -11, -2};
Plane Surface(5) = {5};
Line Loop(6) = {12, -7, -11, 3};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Physical Volume("Odlew") = {1};

Physical Surface("WB_3R_") = {1, 2, 3, 4, 5, 6};
