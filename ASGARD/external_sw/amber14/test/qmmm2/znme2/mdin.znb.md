ZnMe2: 10 steps MD - Pure QM with PM3/ZnB, no temperature coupling
 &cntrl      
   imin =0,
   irest=1,
   ntx=5,
   dt=0.0005,
   ntpr=1,   
   ntb=0, 
   cut=999., 
   ntt=0,    
   nstlim=10,
   ifqnt=1   
 /                                                       
 &qmmm                                                   
   qmmask='@*',                                        
   qm_theory='ZnB', qmcharge=0                           
 /                                                       
EOF   
