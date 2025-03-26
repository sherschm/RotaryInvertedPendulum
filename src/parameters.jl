
##define system parameters
r=0.004 # L-bar radius (in m)
ri=0.0032
dens_Al=2700
dens_per_m=dens_Al*pi*(r^2-ri^2)


#m1=0.01 # L-bar horizontal section mass (in kg)
#m2=0.015 # L-bar vertical section length (in kg)
l1=0.125 # L-bar horizontal section length (in m)
l2=0.264 # L-bar vertical section length (in m)
6
#m1=0.01753
#m2=0.03577
m1=Float16(dens_per_m*l1)*17.4/19.01
m2=Float16(dens_per_m*l2)*17.4/19.01 # L-bar vertical section length (in kg)

m_tip=0.02#0.002

m=m1+m2+m_tip #total L-bar mass
g=9.81 #gravity
#rc=[0;0.75*l1;-0.75*l2] #position of L-bar centre of mass, measured in body frame (see diagrams)
#rc=[0;0.104;-0.0855] #position of L-bar centre of mass, measured in body frame (see diagrams)

rc=(1/m)*(m1*[0;l1/2;0]+m2*[0;l1;-l2/2]+m_tip*[0;l1;-l2])

#rc=[0;CoM_y;CoM_z] #position of L-bar centre of mass, measured in body frame (see diagrams)

function I_tube(r1,r2,h,mi,offset,axis)
    #Inertia of thick-walled cylinders
    if axis == "x"
        Ix=0.5*mi*(r2^2+r1^2)
        Iy=(1/12)*mi*(3*(r2^2+r1^2)+h^2)
        Iz=Iy
    elseif axis == "y"
        Ix=(1/12)*mi*(3*(r2^2+r1^2)+h^2)
        Iy=0.5*mi*(r2^2+r1^2)
        Iz=Ix
    elseif axis == "z"
        Ix=(1/12)*mi*(3*(r2^2+r1^2)+h^2)
        Iy=Ix
        Iz=0.5*mi*(r2^2+r1^2)
    end

    parallel_axis_terms=mi*(offset'*offset*I(3)-offset*offset')

    I_out=Diagonal([Ix;Iy;Iz])+parallel_axis_terms

    return I_out
end

#Inertia tensor of L-bar horizontal section measured in body frame (see diagrams)
#I1=[(1/4)*m1*r^2+(1/3)*m1*l1^2 0 0;0 (1/2)*m1*r^2 0;0 0 (1/4)*m1*r^2+(1/3)m1*l1^2] 
#Inertia tensor of L-bar vertical section measured in body frame (see diagrams), using the parallel axis theorem
#I2=[(1/4)*m2*r^2+(1/12)*m2*l2^2 0 0;0 (1/4)*m2*r^2+(1/12)*m2*l2^2 0;0 0 0.5*m2*r^2]+m2*([0;l1;-l2/2]'*[0;l1;-l2/2]*I(3)-[0;l1;-l2/2]*[0;l1;-l2/2]')
#Combined inertia tensor of L-bar, measured in body frame (see diagrams)
#Ip=I1+I2
I1=I_tube(r,ri,l1,m1,[0;l1/2;0],"y")
I2=I_tube(r,ri,l2,m2,[0;l1;-l2/2],"z")
#I1=[]

I3=m_tip*Diagonal([0;l1^2;l2^2])

Ip=I1+I2+I3

#Ip=Ip/10
#Ip=zeros(3,3)