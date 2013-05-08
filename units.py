#Conversions, unit addition
from __future__ import division
from visual import *
from numpy import *
from math import *
from types import *
#meters in AU


class angle(object):
    radval=0
    mode='0to2pi'
    dispmode='rad'
    __class__='angle'
    def __init__(self,num,dispmode='rad',mode='0to2pi'):
        self.mode=mode
        self.dispmode=dispmode
        if (str(type(num))=="<class 'units.angle'>"):
            self.radval=num.radval
        elif dispmode=='rad':
            self.radval=num
        elif dispmode=='deg':
            self.radval=num*pi/180.
        elif dispmode=='hours':
            self.radval=num*pi/12.
        elif dispmode=='dms':
            self.radval=dms_to_rad(num)
        elif dispmode=='time':
            self.radval=time_to_rad(num)
    def deg(self):
        return (180.)*self.radval/pi
    def rad(self):
        return self.radval
    def dms(self):
        if self.mode=='0to2pi':
            return rad_to_dms(self.radval)
        if self.mode=='-pitopi':
            degval=self.deg()
            sign=1
            if degval<0:
                sign=-1
                degval=abs(degval)
            ans=deg_to_dms(degval)
            return [sign*ans[0],ans[1],ans[2]]
    def hours(self):
        return 12.*self.radval/pi
    def time(self):
        if self.mode=='0to2pi':
            return rad_to_time(self.radval)
        if self.mode=='-pitopi':
            degval=self.deg()
            sign=1
            if degval<0:
                sign=-1
                degval=abs(degval)
            ans=deg_to_time(degval)
            return [sign*ans[0],ans[1],ans[2]]
    def shift(self,mode):
        if mode=='0to2pi':
            while self.radval>2.0*pi:
                self.radval=self.radval-2.0*pi
            while self.radval<0:
                self.radval=self.radval+2.0*pi
        if mode=='-pitopi':
            while self.radval>pi:
                self.radval=self.radval-2.0*pi
            while self.radval<(-pi):
                self.radval=self.radval+2.0*pi
        self.mode=mode
        return angle(self.radval,mode=self.mode)
    def displayMode(self,a):
        self.dispmode=a
        return self
    def __add__(self,a):
        if type(a)==type(angle(0)):
            return angle(self.radval+a.radval).shift(self.mode)
        return angle(self.radval+a).shift(self.mode)
    def __sub__(self,a):
        if type(a)==type(angle(0)):
            return angle(self.radval-a.radval).shift(self.mode)
        return angle(self.radval-a).shift(self.mode)
    def __mul__(self,a):
        return angle(self.radval*a).shift(self.mode)
    def __div__(self,a):
        if type(a)==type(angle(0)):
            return self.radval/a.radval
        return angle(self.radval/a)
    def __truediv__(self,a):
        if type(a)==type(angle(0)):
            return self.radval/a.radval
        return angle(self.radval/a)
    def __str__(self):
        if self.dispmode=='rad':
            return str(self.radval)
        if self.dispmode=='deg':
            return str(self.radval*180./pi)
        if self.dispmode=='hours':
            return str(self.radval*12./pi)
        if self.dispmode=='dms':
            a=self.dms()
            return str(a[0])+" d "+str(a[1])+" m "+str(a[2])+" s"
        if self.dispmode=='time':
            a=self.time()
            return str(a[0])+" h "+str(a[1])+" m "+str(a[2])+" s"
    
    def __repr__(self):
        if self.dispmode=='rad':
            return "angle("+str(self.radval)+",'rad')"
        if self.dispmode=='deg':
            return "angle("+str(self.radval*180./pi)+",'deg')"
        if self.dispmode=='hours':
            return "angle("+str(self.radval*12./pi)+",'hours')"
        if self.dispmode=='dms':
            a=self.dms()
            return "angle(["+str(a[0])+", "+str(a[1])+", "+str(a[2])+"],dispmode='dms')"
        if self.dispmode=='time':
            a=self.time()
            return "angle(["+str(a[0])+", "+str(a[1])+", "+str(a[2])+"],dispmode='time')"



     
def Sin(x):
    """
Sin(x,mode='deg'): returns the sine of angle x measured in units
given in mode. By default, mode is 'deg'. Sin also takes dms, time, and radian
units.
    """
    return sin(x.rad())
def Cos(x,mode='deg'):
    """
Cos(x,mode='deg'): returns the cosine of angle x measured in units
given in mode. By default, mode is 'deg'. Cos also takes dms, time, and radian
units.
    """
    return cos(x.rad())
def Tan(x,mode='deg'):
    """
Tan(x,mode='deg'): returns the tangent of angle x measured in units
given in mode. By default, mode is 'deg'. Tan also takes dms, time, and radian
units.
    """
    return tan(x.rad())
def Asin(x):
    return angle(asin(x))
def Acos(x):
    return angle(acos(x))
def Atan(x):
    return angle(atan(x))

##class uv:
##    u=0
##    n=0
##    def __init__(self,num,uncertainty):
##        self.n=num
##        self.u=uncertainty
##    def __add__(self,other):
##        return uv(self.n+other.n,self.u+other.u)
##    def __radd__(self,other):
##        return uv(self.n+other,self.u)
##    def __sub__(self,other):
##        return uv(self.n-other.n,self.u+other.u)
##    def __rsub__(self,other):
##        return uv(self.n-other,self.u)
##    def __mul__(self,other):
##        if type(other)==type(4) or type(other)==type(4.0):
##            return uv(self.n*other,other*self.u)
##        return uv(self.n*other.n,other.n*self.u+other.u*self.n)
##    def __rmul__(self,other):
##        return uv(self.n*other, self.u*other)
##    def __div__(self,other):
##        return uv(self.n/other.n,other.n*self.u+other.u*self.n)
##    def __rdiv__(self,other):
##        return uv(self.n/other, self.u/other)
##    def __truediv__(self,other):
##        if type(other)==type(4) or type(other)==type(4.0):
##            return uv(self.n/other,self.u/other)
##        return uv(self.n/other.n,other.n*self.u+other.u*self.n)
##    def __pow__(self,other):
##        return uv(self.n**other,(self.n**other)*other*self.u/self.n)
##    def __neg__(self):
##        return uv(-self.n,self.u)
##    def __abs__(self):
##        return uv(abs(self.n),self.u)
##    def __str__(self):
##        return str(self.n)+" pm "+str(self.u)
##    def __float__(self):
##        return float(self.n)
#Conversions:
def dms_to_deg(deg):
    """
dms_to_deg(deg) changes angle deg measured in dms to decimal degrees
    """
    deg[0]=fmod(deg[0],360.0)
    sign=1
    if (deg[0]<0):
        deg[0]=-deg[0]
        return -(deg[0]+deg[1]/60.+deg[2]/3600.) 
    return (deg[0]+deg[1]/60.+deg[2]/3600.)
def deg_to_dms(degree):
    """
deg_to_dms(deg) changes angle deg measured in decimal degrees to dms
    """
    degree=fmod(degree,360.0)
    if degree<0:
        degree+=360.
    deg=[0,0,0]
    deg[0]=floor(degree)
    deg[1]=floor(60*(degree-deg[0]))
    deg[2]=3600.*(degree-deg[0]-deg[1]/60.)
    return deg
def deg_to_rad(degree):
    """
deg_to_rad(deg) changes angle deg measured in decimal degrees to radians
    """
    return degree/360.*2.*pi
def rad_to_deg(rad):
    """
rad_to_deg(rad) changes angle rad measured in radians to decimal degrees
    """
    return rad/2./pi*360.
def deg_to_time(deg):
    """
deg_to_time(deg) changes decimal degrees to time units
    """
    deg=fmod(deg,360.0)
    if deg<0:
        deg+=360
    time=[0,0,0]
    time[0]=floor(deg/360*24)
    time[1]=floor(60*24*(deg-time[0]/24*360)/360)
    time[2]=3600*24.*(deg-time[0]/24.*360.-time[1]/24./60.*360.)/360.
    return time
def time_to_deg(time):
    """
time_to_deg(time) changes time units to decimal degrees
    """
    if time[0]<0:
        return -(abs(time[0])/24*360+time[1]/24/60*360+time[2]/24/3600*360)
    return (time[0]/24.*360.+time[1]/24./60.*360.+time[2]/24./3600.*360.)
def dms_to_time(dms):
    """
dms_to_time(dms) changes dms units to time units
    """
    return deg_to_time(dms_to_deg(dms))
def time_to_dms(time):
    """
time_to_dms(time) changes time units to dms
    """
    return deg_to_dms(time_to_deg(time))
def time_to_rad(t):
    return deg_to_rad(time_to_deg(t))
def rad_to_time(r):
    return deg_to_time(rad_to_deg(r))
def dms_to_rad(dms):
    return deg_to_rad(dms_to_deg(dms))
def rad_to_dms(r):
    return deg_to_dms(rad_to_deg(r))
def n_time(time):
    time[1]+=(time[2]-fmod(time[2],60.0))/60
    time[2]=fmod(time[2],60.0)
    time[0]+=(time[1]-fmod(time[1],60))/60
    time[1]=fmod(time[1],60.0)
    return time





AUm=1.495978707e11
#longitude of Santa Barbara
EL_SB=[angle([34,26,44.0],'dms'),angle([240,20,17.0],'dms')]

PI=angle(pi)

def HourAngle(LST,RA):
    """
HA(LST,RA): gives the hour angle for an object with right ascension RA to LST
    """
    return LST-RA
'''Superfluous now
def add_time(time1,time2):
    return deg_to_time(time_to_deg(time1)+time_to_deg(time2))
def add_dms(deg1,deg2):
    return deg_to_dms(dms_to_deg(deg1)+dms_to_deg(deg2))
'''
def n_dms(deg):
    deg[1]+=(deg[2]-fmod(deg[2],60.0))/60
    deg[2]=fmod(deg[2],60.0)
    deg[0]+=(deg[1]-fmod(deg[2],60))/60
    deg[1]=fmod(deg[1],60.0)
    deg[0]=fmod(deg[0],360.0)
    return deg

def req_to_recliptic(coord):
    e=deg_to_rad(23.43799)
    return [coord[0],cos(e)*coord[1]+coord[2]*sin(e),-coord[1]*sin(e)+coord[2]*cos(e)]
def recliptic_to_req(coord):
    e=deg_to_rad(23.43799)
    return [coord[0],cos(e)*coord[1]-coord[2]*sin(e),coord[1]*sin(e)+coord[2]*cos(e)]
def req_to_eq(vect):
    rmag=mag(vect)
    coord=array(vect)/rmag
    dec=asin(coord[2])
    ra=acos(coord[0]/cos(dec))
    if abs(sin(ra)*cos(dec)-coord[1])>1e-7:
        ra=2*pi-ra
    return [angle(ra),angle(dec)]
def eq_to_req(coord,R):
    ra=coord[0]
    dec=coord[1]
    return array([R*Cos(ra)*Cos(dec),R*Cos(dec)*Sin(ra),R*Sin(dec)])
def eq_to_ecliptic(coord):
    tolerance=0.000001
    a=coord[0]
    d=coord[1]
    e=angle(23.43799,'deg')
    beta=Asin(-Sin(a)*Cos(d)*Sin(e)+Sin(d)*Cos(e))
    l=Asin(1/Cos(beta)*(Sin(a)*Cos(d)*Cos(e)+Sin(d)*Sin(e)))
    if abs(Cos(a)*Cos(d)-Cos(beta)*Cos(l))<tolerance:
        return [beta,l]
    return [beta,angle(pi)-l]

def eq_to_local(coord, LST, EL=EL_SB):
    phi=EL[0]
    d=coord[1]
    HA=LST-coord[0]
    h=Asin(Sin(phi)*Sin(d)+Cos(phi)*Cos(d)*Cos(HA))
    A=Asin(-Cos(d)*Sin(HA)/Cos(h))
    tolerance=1e-6
    if abs(Cos(A)*Cos(h)-Cos(phi)*Sin(d)+Sin(phi)*Cos(d)*Cos(HA))>tolerance:
        A=PI-A
    return [h,A]

def local_to_eq(local,LST,EL=EL_SB):
    phi=EL[0]
    z=PI/2-local[0]
    A=local[1]
    d=Asin(Sin(phi)*Cos(z)+Cos(phi)*Sin(z)*Cos(A))
    HA=Asin(-Cos(PI/2-z)*Sin(A)/Cos(d))
    tolerance=1e-6
    if abs(Cos(d)*Cos(HA)-Cos(phi)*Cos(z)+Sin(phi)*Sin(z)*Cos(A))>tolerance:
        HA=PI-HA
    return [angle(LST-HA),angle(d)]
    
def LAT_to_LST(RA_Sun,time):
    partial=RA_Sun+time
    return partial+angle([12,0,0],'time')
def LST_to_LAT(RA_Sun,LST):
    return LST-PI-RA_Sun
def transit_data_LAT(RA_Sun,coord,EL=EL_SB):
    lat_transit=LST_to_LAT(RA_Sun,coord[0])
    p=EL[0]
    d=coord[1]
    HA=Acos(-Sin(p)*Sin(d)/Cos(p)/Cos(d))
    return [lat_transit-HA,lat_transit,lat_transit+HA]
def LAT_to_UT(LAT,EL=EL_SB):
    UT=LAT-EL[1]
    return UT
def UT_to_LAT(UT,EL=EL_SB):
    LAT=UT+EL[1]
    return LAT    
def angular_dist(coord1,coord2):
    """
angular_dist(coord1,coord2): gives angular distance in dms
between coord1 and coord2
    """
    a1=coord1[0]
    a2=coord2[0]
    d1=coord1[1]
    d2=coord2[1]
    v1=[Cos(a1)*Cos(d1),Sin(a1)*Cos(d1),Sin(d1)]
    v2=[Cos(a2)*Cos(d2),Sin(a2)*Cos(d2),Sin(d2)]
    ans=Acos(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
    return ans
''' Extranuous now
def minus_dms(deg1,deg2):
    """
minus_dms(deg1,deg2): subtracts two angles given in dms
    """
    return deg_to_dms(dms_to_deg(deg1)-dms_to_deg(deg2))
def minus_time(t1,t2):
    """
minus_time(time1,time2): subtracts two angles given in time
    """
    return deg_to_time(time_to_deg(t1)-time_to_deg(t2))
'''
def LST(solar_time,EL=EL_SB,mode='calender'):
    start=JulianDate([2011,7,1,7,0,0])
    if mode=='JD':
        days=solar_time-start
    else:
        days=JulianDate(solar_time)-start
    H=1.0027379038908448*24.*days
    LST=fmod(17.622006749999997+H+EL[1].hours()-EL_SB[1].hours(),24.)
    LST=angle(LST,'hours')
    return LST
    
def JulianDate(cal):
    JD2000=2451544.500000
    y=cal[0]
    month=cal[1]
    day=cal[2]
    hour=cal[3]
    minute=cal[4]
    second=cal[5]
    l=floor((y-2001)/4+1)
    mlist=[31,28,31,30,31,30,31,31,30,31,30,31]
    JDN=366.*l+365.*(y-2000-l)-floor((y-2001)/100)+floor((y-2001)/400)+sum(mlist[0:month-1])+day-1
    if (y%4==0) and ((y%100!=0) or (y%400==0)):
        if month>=3:
            JDN+=1
    return JD2000+JDN+hour/24.+minute/(24.*60.)+second/(24.*3600.)


def eq_to_standard(ad,AD):
    '''
    Converts between equatorial coordinates and local standard coordinates
    '''
    d=ad[1]
    a=ad[0]
    D=AD[1]
    A=AD[0]
    xi=Sin(a-A)/(Sin(D)*Tan(d)+Cos(D)*Cos(a-A))
    eta=(Tan(d)-Tan(D)*Cos(a-A))/(Tan(D)*Tan(d)+Cos(a-A))
    return [xi,eta]
def standard_to_eq(XE,AD):
    '''
    Converts standard flat coordinates to equatorial coordinates
    '''
    xi=XE[0]
    eta=XE[1]
    D=AD[1]
    A=AD[0]
    a=A+Atan(xi/(Cos(D)-eta*Sin(D)))
    d=Atan((eta*Cos(D)+Sin(D))*Sin(a-A)/xi)
    return [a,d]

###converted everything before
def apparentRADec(realRA,realDec,LST):
    '''
    Gives the RA and Dec an observer would see, after correcting for atmospheric refraction
    '''
    n=angle([0,0,58.2],'dms').rad()
    Local=eq_to_local([realRA,realDec],LST)
    z=PI/2.0-Local[0]
    apparentz=solve(lambda x:angle.rad(angle.shift(z-x-n*Tan(x),'-pitopi')),z)
    return local_to_eq([PI/2.-apparentz,Local[1]],LST)
def realRADec(apparentRA,apparentDec,LST):
    '''
    Gives astrometric RA and Dec given apparent RA and Dec, taking into account atmospheric refraction
    '''
    n=angle([0,0,58.2],'dms').rad()
    Local=eq_to_local([apparentRA,apparentDec],LST)
    apparentz=PI/2.0-Local[0]
    z=apparentz+n*Tan(apparentz)
    return local_to_eq([PI/2.-z,Local[1]],LST)


def transit_data_UT(RA_Sun,coord,EL=EL_SB):
    """
transit_data_UT(RA_Sun,coord,EL): return [rise, transit, set] times
of an object with coordinates [Dec, RA]. EL is the [latitude,longitude] location
of the observer.    """
    ans=transit_data_LAT(RA_Sun,coord,EL)
    return [LAT_to_UT(ans[0],EL),LAT_to_UT(ans[1],EL),LAT_to_UT(ans[2],EL)]

def solve(f,guess=0):
    tolerance=1e-10
    oldguess=guess
    while abs(f(guess))>tolerance:
        oldguess=guess
        guess=guess-f(guess)/(deriv(f)(guess))
    return guess

def deriv(f,dx=1e-9):
    return lambda x:(f(x+dx)-f(x-dx))/(2*dx)


