# Plane-plane angle restraint on the base of A1 and base of U2
 &rst  restraint="angle(plane(:1@N9,:1@N1,:1@N3,:1@N7), plane(:2@N1,:2@C4,:2@N3,:2@C6 ))", r1=0., r2=180., r3=180., r4=360., rk2 = 30.,  rk3 = 30.,   /

#restraint 2, featuring new, simpler symmetric parabolic well definition
 &rst restraint="distance(com(:1@N9,:1@N1,:1@N3,:1@N7), com(:2@N1, :2@C2, :2@N3,:2@C4,:2@C5, :2@C6))", r0=5.0, k0=3.0, /
