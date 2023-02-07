#  idb -dbx -I ~/amber10/src/sander -gui ~/amber10/exe/sander
 &rst  restraint="angle(plane(:1@N9,:1@N1,:1@N3,:1@N7), com(:2@N1, :2@C2, :2@N3, :2@C4,:2@C5,:2@C6 ))", r1=0., r2=180., r3=180., r4=360., rk2 = 30.,  rk3 = 30.,   /

#restraint 2
 &rst restraint="distance(com(:1@N9,:1@N1,:1@N3,:1@N7), com(:2@N1, :2@C2, :2@N3,:2@C4,:2@C5, :2@C6))", r1=0., r2=5.0, r3=5.0, r4=100.0, rk2=3.0, rk3=3.0, /
