lc_1 = 0.025;
lc_2 = 0.025;

Point(1) = {0.0, 0.0, 0.0, lc_1};
Point(2) = {0.1, 0.0, 0.0, lc_1};
Point(3) = {0.1, 0.1, 0.0, lc_1};
Point(4) = {0.0, 0.1, 0.0, lc_1};

Point(102) = {0.1, 0.0, 0.0, lc_1};
Point(103) = {0.1, 0.1, 0.0, lc_1};
Point(104) = {0.0, 0.1, 0.0, lc_1};
Point(105) = {0.2, 0.0, 0.0, lc_2};
Point(106) = {0.2, 0.2, 0.0, lc_2};
Point(107) = {0.0, 0.2, 0.0, lc_2};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {102, 105};
Line(105) = {105, 106};
Line(106) = {106, 107};
Line(107) = {107, 104};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(101) = {102, 103, -107, -106, -105, -104};

Plane Surface(1) = {1};
Plane Surface(101) = {101};

Physical Surface("ODLEW") = {1};
Physical Surface("FORMA") = {101};

Physical Line("WB_3R_FORMA") = {105, 106};
Physical Line("WB_4R_ODLEW") = {2, 3};

Periodic Line {2} = {102};
Periodic Line {3} = {103};

