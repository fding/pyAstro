from __future__ import division
from visual import *
from numpy import *
#from astrometry import *
from units import *
#from ephempyexample import *
from random import *
#Solar System, ephemeris
class stuff:
    v=array([0.,0.,0.])
    pos=array([0.,0.,0.])
    mass=0
    a=1
    e=0
    acc=array([0.,0.,0.])
    k1=array([0.,0.,0.])
    k2=array([0.,0.,0.])
    k3=array([0.,0.,0.])
    k4=array([0.,0.,0.])
    j1=array([0.,0.,0.])
    j2=array([0.,0.,0.])
    j3=array([0.,0.,0.])
    j4=array([0.,0.,0.])
    def __init__(self):
        a=0


####
###Planets:
sun=stuff()
sun.pos=array([-6.243886115088563e8,1.062083997699179e8,2.400305101962636e6])
sun.mass=1.9891e30
sun.j1=array([0.,0.,0.])
sun.j2=array([0.,0.,0.])
sun.j3=array([0.,0.,0.])
sun.j4=array([0.,0.,0.])

mercury=stuff()
mercury.mass=3.3022e23
mercury.pos=array([-4.868325826752351e10,2.199652906945730e10,6.200633221988793e9])
mercury.v=array([-3.021965047917072e4,-4.227543190653812e4,-6.802678950812648e3])
mercury.j1=array([0.,0.,0.])
mercury.j2=array([0.,0.,0.])
mercury.j3=array([0.,0.,0.])
mercury.j4=array([0.,0.,0.])

venus=stuff()
venus.mass=4.8685e24
#venus.a=108208930000
#venus.e=0.0068
venus.pos=array([-8.112438268561418e10, 7.110584284762563e10,5.620948075856388e09])
venus.v=array([-2.329321731782643e4,-2.645659210508201e4,9.820917050329498e2])
venus.j1=array([0.,0.,0.])
venus.j2=array([0.,0.,0.])
venus.j3=array([0.,0.,0.])
venus.j4=array([0.,0.,0.])
earth=stuff()
earth.mass=5.9736e24
#earth.a=149598261*1000
#earth.e=0.01671123
earth.pos=array([ -2.629469163081263e10,1.449571069968961e11,-1.207457863666117e6])
earth.v=array([-2.982509096017993e4,-5.315222157007565e03,-7.367084768115095E-01])
earth.j1=array([0.,0.,0.])
earth.j2=array([0.,0.,0.])
earth.j3=array([0.,0.,0.])
earth.j4=array([0.,0.,0.])
moon=stuff()
moon.mass=7.349e22
#moon.a=384399000
#moon.e=0.0529
moon.pos=array([-2.648940524592341e10,1.446321279486018e11,-2.055339080217481e7])
moon.v=array([-2.895706787392246e4,-5.878199883958462e3,7.710557110927274e1])
moon.j1=array([0.,0.,0.])
moon.j2=array([0.,0.,0.])
moon.j3=array([0.,0.,0.])
moon.j4=array([0.,0.,0.])
jupiter=stuff()
jupiter.mass=1.8981e27
#jupiter.a=778547200*1000
#jupiter.e=0.048775
jupiter.pos=array([7.333218281948025e11,9.712434449977794e10,-1.682478939648275e10])
jupiter.v=array([-1.873166210649385e3,1.357358248115560e4,-1.455243707655640e1])
jupiter.j1=array([0.,0.,0.])
jupiter.j2=array([0.,0.,0.])
jupiter.j3=array([0.,0.,0.])
jupiter.j4=array([0.,0.,0.])
mars=stuff()
mars.mass=6.4185e23
#mars.a=227939100*1000
#mars.e=0.093315
mars.pos=array([8.498014333528024e10,-1.930225884994938e11,-6.146530377332136e9])
mars.v=array([2.306604064424532e4, 1.188752890774497e4,-3.170516409506519e2])
mars.j1=array([0.,0.,0.])
mars.j2=array([0.,0.,0.])
mars.j3=array([0.,0.,0.])
mars.j4=array([0.,0.,0.])

saturn=stuff()
saturn.mass=5.68319e26
saturn.pos=array([-1.408584239721548e12,-2.646072651210334e11,6.065554088410646e10])
saturn.v=array([1.265778199845038e3,-9.513741849415672e3,1.147823978528284e2])
saturn.j1=array([0.,0.,0.])
saturn.j2=array([0.,0.,0.])
saturn.j3=array([0.,0.,0.])
saturn.j4=array([0.,0.,0.])

objects=[sun,earth,moon,jupiter,venus]

sun.radius=1e10
venus.radius=1e10
earth.radius=1e10
mercury.radius=1e9
moon.radius=1e10
mars.radius=1e10
jupiter.radius=1e10
moon.color=color.red
earth.color=color.blue
sun.v=array([1.425074259611765,-1.068014319018527e1,-1.402385389576676e-2])


sun.acc=array([0.,0.,0.])
mercury.acc=array([0.,0.,0.])
venus.acc=array([0.,0.,0.])
earth.acc=array([0.,0.,0.])
moon.acc=array([0.,0.,0.])
mars.acc=array([0.,0.,0.])
jupiter.acc=array([0.,0.,0.])
saturn.acc=array([0.,0.,0.])



###
k=0.01720209895
def orbital_elements(rvect,rdot,t,mu=1):
    
    r=mag(rvect)
    v=mag(rdot)
    hvect=cross(rvect,rdot)
    N=cross(array([0,0,1]),hvect)
    evect=cross(rdot,hvect)/mu-rvect/r

    a=1/(2/r-v**2/mu)
    e=mag(evect)
    i=angle(acos(vdot(hvect,array([0,0,1]))/mag(hvect)))
    Omega=angle(acos(vdot(N,array([1,0,0]))/mag(N)))
    if N[1]<0:
        Omega=PI*2.0-Omega
    w=angle(acos(vdot(evect,N)/(e*mag(N))))
    if vdot(cross(N,evect),hvect)<0:
        w=PI*2.0-w
        
    cosnu=1./e*(mag(hvect)**2/mu/r-1)
    E=0
    if abs(r*cosnu/a+e)<=1:
        E=acos(r*cosnu/a+e)
    elif r*cosnu/a+e>1:
        E=0
    elif r*cosnu/a+e<-1:
        E=pi
    if vdot(cross(evect,rvect),hvect)<0:
        E=2*pi-E
    M=E-e*sin(E)
    T=t-M*1/k*sqrt(a**3/mu)
    return [a,e,i.shift('-pitopi').deg(),Omega.deg(),w.deg(),T]



def solve_kepler(M,e):
    f=lambda x:x-e*sin(x)-M
    fderiv=lambda x:1-e*cos(x)
    tolerance=1e-11
    E=M
    while abs(f(E))>tolerance:
        E=E-f(E)/fderiv(E)
    return E
def calculate_orbit(orbit_elements,t,mu=1):
    a=orbit_elements[0]
    e=orbit_elements[1]
    i=deg_to_rad(orbit_elements[2])
    Omega=deg_to_rad(orbit_elements[3])
    w=deg_to_rad(orbit_elements[4])
    T=orbit_elements[5]
    M=k*sqrt(mu/a**3)*(t-T)
    M=fmod(M,2*pi)
    E=solve_kepler(M,e)
    x=a*(cos(E)-e)
    y=a*sqrt(1-e**2)*sin(E)
    r_orb=array([x,y,0])
    M=array([[cos(w)*cos(Omega)-sin(w)*sin(Omega)*cos(i),-sin(w)*cos(Omega)-cos(w)*sin(Omega)*cos(i),sin(Omega)*sin(i)],[cos(w)*sin(Omega)+sin(w)*cos(Omega)*cos(i),-sin(w)*sin(Omega)+cos(w)*cos(Omega)*cos(i),-cos(Omega)*sin(i)],[sin(w)*sin(i),cos(w)*sin(i),cos(i)]])

    r_ecl=dot(M,r_orb)
    return r_ecl


def gravity_simulate(objects,t,dt=1000):
    #Euler's method
    G=6.67e-11
    time=0
    while time<t*24.*3600.:
        time+=dt
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].pos-objects[j].pos
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)

                objects[i].acc+=(-force/objects[i].mass)
                objects[j].acc+=(force/objects[j].mass)
            objects[i].pos=objects[i].pos+objects[i].v*dt+0.5*objects[i].acc*(dt**2)
            objects[i].v=objects[i].v+objects[i].acc*dt
        for i in range(len(objects)):
            objects[i].acc=array([0.,0.,0.])
    return objects

def gravity_simulate2(objects,t,dt=1000):
    #RK 4 integrator
    G=6.67e-11
    time=0
    for i in range(len(objects)):
        objects[i].j1=array([0.,0.,0.])
        objects[i].j2=array([0.,0.,0.])
        objects[i].j3=array([0.,0.,0.])
        objects[i].j4=array([0.,0.,0.])
    while time<t*24.*3600.:
        time+=dt
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].pos-objects[j].pos
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)
                objects[j].j1+=force/objects[j].mass
                objects[i].j1+=-force/objects[i].mass
            objects[i].k2=objects[i].v+objects[i].j1*0.5*dt
            objects[i].new_pos=objects[i].pos+objects[i].k1*0.5*dt
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].new_pos-objects[j].new_pos
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)
                objects[j].j2+=force/objects[j].mass
                objects[i].j2+=-force/objects[i].mass
            objects[i].k3=objects[i].v+objects[i].j2*0.5*dt
            objects[i].new_pos=objects[i].pos+objects[i].k2*0.5*dt
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].new_pos-objects[j].new_pos
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)
                objects[j].j3+=force/objects[j].mass
                objects[i].j3+=-force/objects[i].mass
            objects[i].k4=objects[i].v+dt*objects[i].j3
            objects[i].new_pos=objects[i].pos+objects[i].k3*dt
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].new_pos-objects[j].new_pos
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)
                objects[j].j4+=force/objects[j].mass
                objects[i].j4+=-force/objects[i].mass
        for i in range(len(objects)):
            objects[i].pos+=dt/6.*(objects[i].v+2*objects[i].k2+2*objects[i].k3+objects[i].k4)
            objects[i].v+=dt/6.*(objects[i].j1+2*objects[i].j2+2*objects[i].j3+objects[i].j4)
            objects[i].j1=array([0.,0.,0.])
            objects[i].j2=array([0.,0.,0.])
            objects[i].j3=array([0.,0.,0.])
            objects[i].j4=array([0.,0.,0.])
    return objects

def gravity_simulate3(objects,t,dt=1000):
    #Verlet integrator
    #G=6.674e-11
    G=6.6727005490826302397214410070979e-11
    time=0
    dt=dt/10.
    for i in range(len(objects)):
        objects[i].pos1=copy(objects[i].pos)
        objects[i].pos2=copy(objects[i].pos)

    for count in range(10):
        for i in range(len(objects)):
            objects[i].acc=array([0.,0.,0.])
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].pos1-objects[j].pos1
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)
                objects[i].acc+=(-force/objects[i].mass)
                objects[j].acc+=(force/objects[j].mass)
            objects[i].v=copy(objects[i].v+objects[i].acc*dt)  
            objects[i].pos1=copy(objects[i].pos1+objects[i].v*dt)     #-0.5*objects[i].acc*dt*dt
    dt=dt*10.
    time+=dt
    while abs(time)<abs(t*24.*3600.):
        time+=dt
        for i in range(len(objects)):
            objects[i].acc=array([0.,0.,0.])
        for i in range(len(objects)):
            for j in range(i+1,len(objects)):
                dist=objects[i].pos1-objects[j].pos1
                magd=mag(dist)
                force= G*objects[i].mass*objects[j].mass*dist/(magd*magd*magd)
                objects[i].acc+=(-force/objects[i].mass)
                objects[j].acc+=(force/objects[j].mass)
            objects[i].pos=copy(2*objects[i].pos1-objects[i].pos2+objects[i].acc*dt*dt)
            objects[i].pos2=copy(objects[i].pos1)
            objects[i].pos1=copy(objects[i].pos)
    for i in range(len(objects)):
        objects[i].v=copy((objects[i].pos1-objects[i].pos2)/dt)
        objects[i].pos=copy(objects[i].pos1)

    return objects

def solar_system_simulate(t,T=2455562.500000,extra_obj=[],stepnumber=1e4):
    #Start of simulation is January 1, 2011
    dt=24*3600.*(t-T)/stepnumber
    print dt

    if T>2455562.500000:
        objects=gravity_simulate3(objects,T-2455562.500000,dt) #steps it up to the time of start
    objects+=extra_obj
    return gravity_simulate3(objects,t-T,dt)
                
'''    
def rk4(xvect,vvect,a,dt):
    k1=dt*a(xvect)
    j1=dt*vvect
    k2=dt*a(xvect+j1/2)
    j2=dt*(vvect+dt*k1/2.)
    k3=dt*a(xvect+j2/2.)
    j3=dt*(vvect+dt*k2/2.)
    k4=dt*a(xvect+j3)
    j4=dt*(vvect+dt*k3)
    vnew=vvect+(k1+2*k2+2*k3+k4)/6
    xnew=xvect+(j1+2*j2+2*j3+j4)/6
    return [xnew,vnew]
'''    

def ephemerides(orbit_elements,t,mu=1,EL=EL_SB):
    c=173.1446327
    ephem=Ephemeris('405')
    T=orbit_elements[5]
    r_asteroid=calculate_orbit(orbit_elements,t,mu)
    r_earth=-ephem.position(t,10,2)
    rho=recliptic_to_req(r_asteroid)-r_earth
    t=t-mag(rho)/c
    r_asteroid=calculate_orbit(orbit_elements,t,mu)
    rho=recliptic_to_req(r_asteroid)-r_earth
    rho=center_to_surface(rho,EL[0],LST(t,mode='JD'))
    ans=req_to_eq(rho)
    return ans

def ephemeris(orb,t):
    ans=ephemerides(orb,t)
    return [ans[0].time(),ans[1].dms()]

def f(rvect,vvect,t,mu=1.0):
    r0=mag(rvect)
    v2=vdot(vvect,vvect)
    rv=vdot(rvect,vvect)
    return 1-1/(2*r0*r0*r0)*t*t+rv/(2*r0**5)*t**3+1/24.*(3/(r0**3)*(v2/r0**2)-15.*rv*rv/(r0**7)-2./(r0**6))*t**4+1/8.*rv/(r0**5)*(3*v2/r0**2-2/r0**3-7*rv*rv/(r0**7))*t**5

def g(rvect,vvect,t,mu=1.0):
    r0=mag(rvect)
    v2=vdot(vvect,vvect)
    rv=vdot(rvect,vvect)
    return t-1/(6*r0*r0*r0)*t*t*t+rv/(4*r0**5)*(t**4)+1/120.*(9/(r0**3)*(v2/r0**2)-45.*rv*rv/(r0**7)-8./(r0**6))*t**5

lorb_list=[]
lrx_list=[]
lry_list=[]
lrz_list=[]
lvx_list=[]
lvy_list=[]
lvz_list=[]
def DetermineOrbit(RADecFilename):
    RADecFile=open(RADecFilename)
    coord_perm_t_list=[]
    ephem=Ephemeris('405')
    uncertainty_list=[]
    for line in RADecFile:
        if line[0]=='%':
            continue
        if line=='' or line==' ':
            continue
        line_parts=line.split(' | ')
        RAh=float(line_parts[1])
        RAm=float(line_parts[2])
        RAs=float(line_parts[3])
        Rau=float(line_parts[4])
        Decd=float(line_parts[5])
        Decm=float(line_parts[6])
        Decs=float(line_parts[7])
        Decu=float(line_parts[8])
        t=float(line_parts[0])
        coord_perm_t_list.append([t,angle([RAh,RAm,RAs],'time'),angle([Decd,Decm,Decs],'dms')])
        uncertainty_list.append([Rau,Decu])
        #uncertainty_list.append([0,0])
    #r2,v2,t2=GaussianMethod(coord_list[0],coord_list[1],coord_list[2],t_list[0],t_list[1],t_list[2])
    N=100
    r2,v2,t2=LaplaceMethod(coord_perm_t_list[0][1:],coord_perm_t_list[1][1:],coord_perm_t_list[2][1:],coord_perm_t_list[3][1:],coord_perm_t_list[4][1:],coord_perm_t_list[0][0],coord_perm_t_list[1][0],coord_perm_t_list[2][0],coord_perm_t_list[3][0],coord_perm_t_list[4][0])
    print orbital_elements(r2,v2,t2) #Inserted on 10/3/2011
    r2,v2,t2=differential_corrections(r2,v2,t2,coord_perm_t_list[5:])
    r2,v2,t2=differential_corrections(r2,v2,t2,coord_perm_t_list[5:])
    r2,v2,t2=differential_corrections(r2,v2,t2,coord_perm_t_list[5:])
    r2,v2,t2=differential_corrections(r2,v2,t2,coord_perm_t_list[5:])
    r2,v2,t2=differential_corrections(r2,v2,t2,coord_perm_t_list[5:])
    r2,v2,t2=differential_corrections(r2,v2,t2,coord_perm_t_list[5:])
    real_orb_el=orbital_elements(r2,v2,t2)
    print "Preliminary Results: "
    print real_orb_el  
    print "r: ", r2
    print "v: ", v2
    print "t: ", t2
    for i in range(N):
        coord_t_list=[]
        print i+1
        for j in range(len(coord_perm_t_list)):
            ra=angle(gauss(coord_perm_t_list[j][1].hours(),uncertainty_list[j][0]/3600.),'hours')
            dec=angle(gauss(coord_perm_t_list[j][2].deg(),uncertainty_list[j][1]/3600.),'deg')
            coord_t_list.append([coord_perm_t_list[j][0],ra,dec])
        r2,v2,t2=LaplaceMethod(coord_t_list[0][1:],coord_t_list[1][1:],coord_t_list[2][1:],coord_t_list[3][1:],coord_t_list[4][1:],coord_t_list[0][0],coord_t_list[1][0],coord_t_list[2][0],coord_t_list[3][0],coord_t_list[4][0])
        r2,v2,t2=differential_corrections(r2,v2,t2,coord_t_list[5:])
        r2,v2,t2=differential_corrections(r2,v2,t2,coord_t_list[5:])
        r2,v2,t2=differential_corrections(r2,v2,t2,coord_t_list[5:])
        r2,v2,t2=differential_corrections(r2,v2,t2,coord_t_list[5:])
        r2,v2,t2=differential_corrections(r2,v2,t2,coord_t_list[5:])
        r2,v2,t2=differential_corrections(r2,v2,t2,coord_t_list[5:])
        lorb_list.append(orbital_elements(r2,v2,t2))
        lrx_list.append(r2[0])
        lry_list.append(r2[1])
        lrz_list.append(r2[2])
        lvx_list.append(v2[0])
        lvy_list.append(v2[1])
        lvz_list.append(v2[2])
    
    rxdev=std(lrx_list)
    rydev=std(lry_list)
    rzdev=std(lrz_list)
    vxdev=std(lvx_list)
    vydev=std(lvy_list)
    vzdev=std(lvz_list)
    orb_list=array(lorb_list)
    amean=mean(orb_list[:,0])
    astddev=std(orb_list[:,0])
    emean=mean(orb_list[:,1])
    estddev=std(orb_list[:,1])
    imean=mean(orb_list[:,2])
    istddev=std(orb_list[:,2])
    Omean=mean(orb_list[:,3])
    Ostddev=std(orb_list[:,3])
    wmean=mean(orb_list[:,4])
    wstddev=std(orb_list[:,4])
    Tmean=mean(orb_list[:,5])
    Tstddev=std(orb_list[:,5])
    print "Uncertainty: ", rxdev,rydev,rzdev,vxdev,vydev,vzdev
    return real_orb_el, [astddev,estddev,istddev,Ostddev,wstddev,Tstddev],[amean,emean,imean,Omean,wmean,Tmean]

#Test Arguments:
#([[-7,20,24.3],[17,06,26.13]],[[-11,28,26.5],[17,20,20.02]],[[-18,13,4.4],[17,49,47.27]],array([-0.35795,0.87291,0.37843]),array([-0.46601,0.82851,0.35918]),array([-0.63547,0.72609,0.31478]),2455025.503,2455032.498,2455044.499)
#Should get 0.4196,-1.093,-0.4137

#([angle([17,06,26.13],'time'),angle([-7,20,24.3],'dms')],[angle([17,20,20.02],'time'),angle([-11,28,26.5],'dms')],[angle([17,49,47.27],'time'),angle([-18,13,4.4],'dms')],array([-0.35795,0.87291,0.37843]),array([-0.46601,0.82851,0.35918]),array([-0.63547,0.72609,0.31478]),2455025.503,2455032.498,2455044.499)

def divPoly(poly,root):
    result=zeros(len(poly)-1)
    for i in range(len(poly)-1):
        poly[i+1]=poly[i+1]+poly[i]*root
        result[i]=poly[i]
    return lambda r: sum(result[i]*(r**(len(poly)-2-i)) for i in range(len(poly)-1)),result
def fact(n):
    ans=1
    if n==0:
        return ans
    else:
        return n*fact(n-1)

def LaplaceMethod(coord1,coord2,coord3,coord4,coord5,t1,t2,t3,t4,t5,corrected=False):
#t4-t2 should be small, but t5-t1 should be larger
    speedc=173.1446327 #speed of light in AU per day
    AUm=1.495978707e11
    
    latitude=EL_SB[0]
    ephem=Ephemeris('405')
    tau1=k*(t2-t3)
    tau3=k*(t4-t3)
    rho1hat=eq_to_req(coord1,1.00)
    rho2hat=eq_to_req(coord2,1.00)
    rho3hat=eq_to_req(coord3,1.00)
    rho4hat=eq_to_req(coord4,1.00)
    rho5hat=eq_to_req(coord5,1.00)
    drho3hat=((rho2hat-rho3hat)*tau3*tau3-(rho4hat-rho3hat)*tau1*tau1)/(tau1*tau3*(tau3-tau1))
    tau1=k*(t1-t3)
    tau3=k*(t5-t3)
    ddrho3hat=-2*((rho1hat-rho3hat)*tau3-(rho5hat-rho3hat)*tau1)/(tau1*tau3*(tau3-tau1))
    '''
    b=array([rho1hat,rho2hat,rho4hat,rho5hat])
    A=zeros([4,4])
    for i in range(4):
        for j in range(4):
            if i==0:
                A[i,j]=1/fact(j)*(k*(t1-t3))**j
            if i==1:
                A[i,j]=1/fact(j)*(k*(t2-t3))**j
            if i==2:
                A[i,j]=1/fact(j)*(k*(t4-t3))**j
            if i==3:
                A[i,j]=1/fact(j)*(k*(t5-t3))**j
    ddrho3hat=linalg.solve(A,b)[2]
    print ddrho3hat
    '''
    R3center=ephem.position(t3,10,2)
    R3=ephem.position(t3,10,2)-obs_vector(latitude,LST(t3,mode='JD'))/AUm
    dt=1e-6
    R3_2=ephem.position(t3+dt,10,2)-obs_vector(latitude,LST(t3+dt,mode='JD'))/AUm
    R3_1=ephem.position(t3-dt,10,2)-obs_vector(latitude,LST(t3-dt,mode='JD'))/AUm
    dR3=1/k*(R3_2-R3_1)/(2*dt)
    #ddR3=-1/(k*k)*(R3_2+R3_1-2*R3)/(dt**2)
    ddR3=-1.000003/mag(R3center)**3*R3center
    
    A=TP(drho3hat,ddR3,rho3hat)/TP(drho3hat,ddrho3hat,rho3hat)
    B=TP(drho3hat,R3,rho3hat)/TP(drho3hat,ddrho3hat,rho3hat)
    a=-(A**2-2*A*dot(rho3hat,R3)+dot(R3,R3))
    b=-(2*A*B-2*B*dot(rho3hat,R3))
    c=-B*B
    r=solve(lambda r:r**8+a*r**6+b*r**3+c,1.5)
    coef=[1.0,0,a,0,0,b,0,0,c]
    while r<0 or abs(r)<1.1:
        newpoly,coef=divPoly(coef,r)
        r=solve(newpoly,1.7)
    rho3=A+B/r**3
    C=TP(ddR3,ddrho3hat,rho3hat)/TP(drho3hat,ddrho3hat,rho3hat)
    D=TP(R3,ddrho3hat,rho3hat)/TP(drho3hat,ddrho3hat,rho3hat)
    drho3=C+D/r**3
    rvect=(rho3*rho3hat-R3)
    vvect=(drho3*rho3hat+rho3*drho3hat-dR3)
    '''
    print 'geovects:',obs_vector(latitude,LST(t3-dt,mode='JD'))/AUm,obs_vector(latitude,LST(t3,mode='JD'))/AUm,obs_vector(latitude,LST(t3+dt,mode='JD'))/AUm
    print 'rho3hat,drho3hat','ddrho3hat',rho3hat,drho3hat,ddrho3hat
    print 'R3,dR3,ddR3', R3, dR3, ddR3
    print 'A', A
    print 'B', B
    print 'a', a
    print 'b', b
    print 'c', c
    print rvect,vvect,t3-rho3/c
    '''
    return array(req_to_recliptic(rvect)),array(req_to_recliptic(vvect)),t3-rho3/speedc





  
def GaussianMethod(coord1,coord2,coord3,t1,t2,t3):
    c=173.1446327 #speed of light in AU per day
    AUm=1.495978707e11
    r1guess=array([-0.81477363, -0.87351652,  0.50165778])
    v1guess=array([0,0,0])
    r2guess=array([-0.81477363, -0.87351652 , 0.50165778])
    v2guess=array([0,0,0])
    r3guess=array([-0.81477363, -0.87351652,  0.50165778])
    v3guess=array([0,0,0])
    latitude=EL_SB[0]
    rho1hat=eq_to_req(coord1,1.0)
    rho2hat=eq_to_req(coord2,1.0)
    rho3hat=eq_to_req(coord3,1.0)
    g1=obs_vector(latitude,LST(t1,mode='JD'))
    g2=obs_vector(latitude,LST(t2,mode='JD'))
    g3=obs_vector(latitude,LST(t3,mode='JD'))
    tau1=k*(t1-t2)
    tau3=k*(t3-t2)
    count=0
    count2=0

    ephem=Ephemeris('405')
    R1=ephem.position(t1,10,2)-g1/AUm
    R2=ephem.position(t2,10,2)-g2/AUm
    R3=ephem.position(t3,10,2)-g3/AUm
    
    while True:
        oldr2=copy(r2guess)
        f1=f(r2guess,v2guess,tau1)
        f3=f(r2guess,v2guess,tau3)
        g1=g(r2guess,v2guess,tau1)
        g3=g(r2guess,v2guess,tau3)
        if count==0:
            T3=-tau1
            T1=tau3
            T2=tau3-tau1
            B1=1/12.*(T1*T3+T2*(T3-T1))
            B2=1/12.*(T1*T3+T2**2)
            B3=1/12.*(T1*T3-T2*(T3-T1))
            C1=T1*(1+B1/mag(r1guess)**3)
            C2=T3*(1+B3/mag(r3guess)**3)
            C3=T2*(1-B2/mag(r2guess)**3)
            a1=C1/C3
            a3=C2/C3
        else:
            a1=g3/(f1*g3-f3*g1)
            a3=-g1/(f1*g3-f3*g1)
        U2=Tan(coord2[0])
        V2=Tan(coord2[1])/Cos(coord2[0])
        U1=Tan(coord1[0])
        V1=Tan(coord1[1])/Cos(coord1[0])
        U3=Tan(coord3[0])
        V3=Tan(coord3[1])/Cos(coord3[0])
        D=(U1-U2)*(V3-V2)-(U3-U2)*(V1-V2)

        P1=R1[1]-U1*R1[0]
        Q1=R1[2]-V1*R1[0]
        P2=R2[1]-U2*R2[0]
        Q2=R2[2]-V2*R2[0]
        P3=R3[1]-U3*R3[0]
        Q3=R3[2]-V3*R3[0]
        P=a1*P1-P2+a3*P3
        Q=a1*Q1-Q2+a3*Q3
        x1=(P*(V3-V2)-Q*(U3-U2))/(a1*D)
        x3=(Q*(U1-U2)-P*(V1-V2))/(a3*D)
        r1guess=array([x1,U1*x1-P1,V1*x1-Q1])
        r3guess=array([x3,U3*x3-P3,V3*x3-Q3])
        y2=a1*(U1*x1-P1)+a3*(U3*x3-P3)
        z2=a1*(V1*x1-Q1)+a3*(V3*x3-Q3)
        r2guess=array([(P2+y2)/U2,y2,z2])
        #rho1=(a1*TP(R1,rho2hat,rho3hat)-TP(R2,rho2hat,rho3hat)+a3*TP(R3,rho2hat,rho3hat))/(a1*TP(rho1hat,rho2hat,rho3hat))
        #rho2=(a1*TP(rho1hat,R1,rho3hat)-TP(rho1hat,R2,rho3hat)+a3*TP(rho1hat,R3,rho3hat))/(-TP(rho1hat,rho2hat,rho3hat))
        #rho3=(a1*TP(rho2hat,R1,rho1hat)-TP(rho2hat,R2,rho1hat)+a3*TP(rho2hat,R3,rho1hat))/(a3*TP(rho2hat,rho3hat,rho1hat))
        rho1=mag(r1guess+R1)
        rho2=mag(r2guess+R2)
        rho3=mag(r3guess+R3)
        #r2guess=rho2*rho2hat-R2
        #r3guess+R3=rho3*rho3hat
        if count2==0:
            count2+=1
            v2guess=((tau3**2)*(r1guess-r2guess)-tau1**2*(r3guess-r2guess))/(tau1*tau3*(tau3-tau1))
        else:
            v2guess=f3/(g1*f3-g3*f1)*r1guess-f1/(g1*f3-g3*f1)*r3guess
        tau1=k*(t1-rho1/c-t2+rho2/c)
        tau3=k*(t3-rho3/c-t2+rho2/c)
        v2guess=1/g1*(r1guess-f1*r2guess)
        if (mag(r2guess-oldr2)<1e-10):
            break
    return array(req_to_recliptic(r2guess)),array(req_to_recliptic(v2guess)),t2-rho2/c


def TP(v1,v2,v3):
    return dot(cross(v1,v2),v3)
'''
dec1=angle([-7,20,24.3],'dms')
dec2=angle([-11,28,26.5],'dms')
dec3=angle([-18,13,4.4],'dms')
ra1=angle([17,06,26.13],'time')
ra2=angle([17,20,20.02],'time')
ra3=angle([17,49,47.27],'time')
GaussianMethod([ra1,dec1],[ra2,dec2],[ra3,dec3],2455025.503,2455032.498,2455044.499)
'''
def differential_corrections(r0,v0,t0,coord_t_list):
    n=len(coord_t_list)
    
    A=zeros([2*n,7])
    b=zeros([2*n])
    diff=1e-1
    dx=diff*array([1,0,0])
    dy=diff*array([0,1,0])
    dz=diff*array([0,0,1])
    
    for i in range(n):
        t=coord_t_list[i][0]
        ephem11=ephemerides(orbital_elements(r0+dx,v0, t0), t)
        ephem12=ephemerides(orbital_elements(r0-dx,v0, t0), t)
        ephem21=ephemerides(orbital_elements(r0+dy,v0, t0), t)
        ephem22=ephemerides(orbital_elements(r0-dy,v0, t0), t)
        ephem31=ephemerides(orbital_elements(r0+dz,v0, t0), t)
        ephem32=ephemerides(orbital_elements(r0-dz,v0, t0), t)
        ephem41=ephemerides(orbital_elements(r0,v0+dx, t0), t)
        ephem42=ephemerides(orbital_elements(r0,v0-dx, t0), t)
        ephem51=ephemerides(orbital_elements(r0,v0+dy, t0), t)
        ephem52=ephemerides(orbital_elements(r0,v0-dy, t0), t)
        ephem61=ephemerides(orbital_elements(r0,v0+dz, t0), t)
        ephem62=ephemerides(orbital_elements(r0,v0-dz, t0), t)
        b[i]=angle.shift(coord_t_list[i][1]-ephemerides(orbital_elements(r0,v0, t0), t)[0],'-pitopi').rad()
        A[i,0]=angle.shift(ephem11[0]-ephem12[0],'-pitopi').rad()/(2*diff)
        A[i,1]=angle.shift(ephem21[0]-ephem22[0],'-pitopi').rad()/(2*diff)
        A[i,2]=angle.shift(ephem31[0]-ephem32[0],'-pitopi').rad()/(2*diff)
        A[i,3]=angle.shift(ephem41[0]-ephem42[0],'-pitopi').rad()/(2*diff)
        A[i,4]=angle.shift(ephem51[0]-ephem52[0],'-pitopi').rad()/(2*diff)
        A[i,5]=angle.shift(ephem61[0]-ephem62[0],'-pitopi').rad()/(2*diff)

        b[n+i]=angle.shift(coord_t_list[i][2]-ephemerides(orbital_elements(r0,v0, t0), t)[1],'-pitopi').rad()
        A[n+i,0]=angle.shift(ephem11[1]-ephem12[1],'-pitopi').rad()/(2*diff)
        A[n+i,1]=angle.shift(ephem21[1]-ephem22[1],'-pitopi').rad()/(2*diff)
        A[n+i,2]=angle.shift(ephem31[1]-ephem32[1],'-pitopi').rad()/(2*diff)
        A[n+i,3]=angle.shift(ephem41[1]-ephem42[1],'-pitopi').rad()/(2*diff)
        A[n+i,4]=angle.shift(ephem51[1]-ephem52[1],'-pitopi').rad()/(2*diff)
        A[n+i,5]=angle.shift(ephem61[1]-ephem62[1],'-pitopi').rad()/(2*diff)
    sol=linalg.lstsq(A,b)
    r0=r0+array(sol[0][:3])
    v0=v0+array(sol[0][3:6])
    return r0,v0,t0
    
    
