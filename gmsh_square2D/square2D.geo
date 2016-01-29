lc = 0.025;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.1, 0.0, 0.0, lc};
Point(3) = {0.1, 0.1, 0.0, lc};
Point(4) = {0.0, 0.1, 0.0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};

Plane Surface(6) = {5};

Physical Surface("ODLEW") = {6};

Physical Line("WB_3R_ODLEW") = {1, 2, 3, 4};

