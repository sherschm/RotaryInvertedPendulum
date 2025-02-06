
##define system parameters
r=0.008 # L-bar radius (in m)
m1=0.01 # L-bar horizontal section mass (in kg)
m2=0.015 # L-bar vertical section length (in kg)
l1=0.125 # L-bar horizontal section length (in m)
l2=0.255 # L-bar vertical section length (in m)
m=m1+m2 #total L-bar mass
g=9.81 #gravity
rc=[0;0.75*l1;-0.75*l2] #position of L-bar centre of mass, measured in body frame (see diagrams)

#Inertia tensor of L-bar horizontal section measured in body frame (see diagrams)
I1=[(1/4)*m1*r^2+(1/3)*m1*l1^2 0 0;0 (1/2)*m1*r^2 0;0 0 (1/4)*m1*r^2+(1/3)m1*l1^2] 
#Inertia tensor of L-bar vertical section measured in body frame (see diagrams), using the parallel axis theorem
I2=[(1/4)*m2*r^2+(1/12)*m2*l2^2 0 0;0 (1/4)*m2*r^2+(1/12)*m2*l2^2 0;0 0 0.5*m2*r^2]+m2*([0;l1;-l2/2]'*[0;l1;-l2/2]*I(3)-[0;l1;-l2/2]*[0;l1;-l2/2]')
#Combined inertia tensor of L-bar, measured in body frame (see diagrams)
Ip=I1+I2