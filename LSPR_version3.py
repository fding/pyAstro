from __future__ import division
from numpy import *
import pyfits
from numpy.random import *
import numpy.linalg as linalg
from Astro_functions import *
from random import *


def centroid(a):
    y_len=len(a)
    x_len=len(a[0])
    x_start=-(x_len-1)/2
    y_start=(y_len-1)/2
    x=0
    y=0
    s=0
    u_x=0
    u_y=0
    u_s=0
    for i in range(x_len):
        for j in range(y_len):
            s+=a[j,i]
            u_s+=sqrt(abs(a[j,i]))
            y+=a[j,i]*(y_start-j)
            x+=a[j,i]*(x_start+i)
            u_x+=sqrt(abs(a[j,i]))*abs(y_start-j)
            u_y+=sqrt(abs(a[j,i]))*abs(x_start+i)
    u_x=abs(x/s*(u_x/x+u_s/s))
    u_y=abs(y/s*(u_y/y+u_s/s))
    x=x/s
    y=y/s
    return [[x,u_x],[y,u_y]]
def bilinear_interpolation(pt_list):
    xsum=sum(pt_list[:,0])
    x2sum=sum(pt_list[:,0]**2)
    ysum=sum(pt_list[:,1])
    y2sum=sum(pt_list[:,1]**2)
    zsum=sum(pt_list[:,2])
    xysum=sum(pt_list[:,0]*pt_list[:,1])
    xzsum=sum(pt_list[:,0]*pt_list[:,2])
    yzsum=sum(pt_list[:,1]*pt_list[:,2])
    A=array([[x2sum,xysum,xsum],[xysum,y2sum,ysum],[xsum,ysum,len(pt_list)]])
    b=array([xzsum,yzsum,zsum])
    coef=linalg.solve(A,b)
    return [coef[0],coef[1],coef[2]]

def eq_to_standard(da,DA):
    d=da[0]
    a=da[1]
    D=DA[0]
    A=DA[1]
    xi=Sin(a-A)/(Sin(D)*Tan(d)+Cos(D)*Cos(a-A))
    eta=(Tan(d)-Tan(D)*Cos(a-A))/(Tan(D)*Tan(d)+Cos(a-A))
    return [xi,eta]
def standard_to_eq(XE,DA):
    xi=XE[0]
    eta=XE[1]
    D=DA[0]
    A=DA[1]
    a=A+rad_to_deg(atan(xi/(Cos(D)-eta*Sin(D))))
    d=rad_to_deg(atan((eta*Cos(D)+Sin(D))*Sin(a-A)/xi))
    return [d,a]
def apparentRADec(realRA,realDec,LST):
    n=dms_to_rad([0,0,58.2])
    Local=eq_to_local([realDec,realRA],LST)
    z=pi/2-dms_to_rad(Local[0])
    apparentz=solve(lambda x:z-x-n*tan(x),z)
    return local_to_eq([rad_to_dms(pi/2-apparentz),Local[1]],LST)
def realRADec(apparentRA,apparentDec,LST):
    n=dms_to_rad([0,0,58.2])
    Local=eq_to_local([apparentDec,apparentRA],LST)
    apparentz=pi/2-dms_to_rad(Local[0])
    z=apparentz+n*tan(apparentz)
    return local_to_eq([rad_to_dms(pi/2-z),Local[1]],LST)
def deriv(f):
    dx=1e-9
    return lambda x:(f(x+dx)-f(x))/dx

def solve(f,guess=0):
    tolerance=1e-9
    while abs(f(guess))>tolerance:
        guess=guess-f(guess)/deriv(f)(guess)
    return guess
    
def LSPR(filename,imagefile):
    """
LSPR uses the method of least squares to find the transform
taking x,y coordinates in an image to the corresponding RA, Dec coordinates.

In the file named filename, a nine column table should be given,
with columns separated by ' | ':
Star name | x coordinate | y coordinate | RA hours | RA min | RA seconds
| Dec degree | Dec mintes | Dec seconds
The first row should contain useful information about the file.

The function will compute a least square best fit linear equation
to the specified (x,y) coordinates and (RA,Dec) coordinates and
outputs the equation and residuals of the interpolated RA Dec coordinates
from the actual RA and Dec given in the files, in the format:
Star name | RA residual ([h,m,s]) | Dec residual ([d,m,s])

Finally, LSPR will return a pair of functions RA_eq, Dec_eq.
Each function could take two coordinates x and y, and output the
RA and Dec respectively.
    """
    image=pyfits.getdata(imagefile)
    center=[(image.shape[0]-1)/2,(image.shape[1]-1)/2]
    hdr=pyfits.getheader(imagefile)
    observation_date=hdr['Date-Obs'].split('-')
    year=int(observation_date[0])
    month=int(observation_date[1])
    day=int(observation_date[2].split('T')[0])
    observation_time=hdr['Time-OBS'].split(':')
    hour=int(observation_time[0])
    minute=int(observation_time[1])
    second=int(float(observation_time[2]))
    time_diff=48.5
    observation_LST=LST([year,month,day,hour,minute,second+time_diff])
    
    f=open(filename)
    data1=[]
    data2=[]
    names=[]
    i=0
    for line in f:
        if i==0:
            i+=1
            continue
        str_list=line.split(" | ")
        if str_list[0][0]=='%':
            continue
        names.append(str_list[0])
        x=float(str_list[1])
        y=float(str_list[2])
        RaH=float(str_list[3])
        RaM=float(str_list[4])
        RaS=float(str_list[5])
        DecD=float(str_list[6])
        DecM=float(str_list[7])
        DecS=float(str_list[8])
        [Dec,RA]=apparentRADec([RaH,RaM,RaS],[DecD,DecM,DecS],observation_LST)
        data1.append([x,y,time_to_deg(RA)])
        data2.append([x,y,dms_to_deg(Dec)])
    
    coef1=bilinear_interpolation(array(data1))
    coef2=bilinear_interpolation(array(data2))
    centerRA0=coef1[0]*center[0]+coef1[1]*center[1]+coef1[2]
    centerDec0=coef2[0]*center[0]+coef2[1]*center[1]+coef2[2]
    centerrange=1.0
    centerRA_best=centerRA0
    centerDec_best=centerDec0
    bestchi2=0
    RA_eq_best=lambda x,y: 0
    Dec_eq_best=lambda x,y:0
    coef1best=[]
    coef2best=[]
    for count in range(300):
        if count==0:
            centerRA=copy(centerRA0)
            centerDec=copy(centerDec0)
        else:
            centerRA=centerRA_best+gauss(0,.30)*1/3600.
            centerDec=centerDec_best+gauss(0,.30)*1/3600.
        
        xi_list=[]
        eta_list=[]
        for i in range(len(data1)):
            standard=eq_to_standard([data2[i][2],data1[i][2]],[centerDec,centerRA])
            xi_list.append([data1[i][0],data1[i][1],standard[0]])
            eta_list.append([data1[i][0],data1[i][1],standard[1]])
        coef1=bilinear_interpolation(array(xi_list))
        coef2=bilinear_interpolation(array(eta_list))

        xy_to_xi=lambda x,y: coef1[0]*x+coef1[1]*y+coef1[2]
        xy_to_eta=lambda x,y: coef2[0]*x+coef2[1]*y+coef2[2]
    
        eq1_RA=lambda x,y: standard_to_eq([xy_to_xi(x,y),xy_to_eta(x,y)],[centerDec,centerRA])[1]    
        eq2_Dec=lambda x,y: standard_to_eq([xy_to_xi(x,y),xy_to_eta(x,y)],[centerDec,centerRA])[0] 
        RA_eq=lambda x,y: realRADec(deg_to_time(eq1_RA(x,y)),deg_to_dms(eq2_Dec(x,y)),observation_LST)[1]
        Dec_eq=lambda x,y: realRADec(deg_to_time(eq1_RA(x,y)),deg_to_dms(eq2_Dec(x,y)),observation_LST)[0]
        ra_chi2=0
        dec_chi2=0
        for i in range(len(data1)):
            temp=realRADec(deg_to_time(data1[i][2]),deg_to_dms(data2[i][2]),observation_LST)
            RA_P=RA_eq(data1[i][0],data1[i][1])
            RA_A=temp[1]
            Dec_P=Dec_eq(data2[i][0],data2[i][1])
            Dec_A=temp[0]
            dec_res=dms_to_deg(Dec_A)-dms_to_deg(Dec_P)
            dec_chi2+=dec_res**2
            ra_res=time_to_deg(RA_A)-time_to_deg(RA_P)
            ra_chi2+=ra_res**2
        tot=min(ra_chi2,dec_chi2)
        if count==0:
            print sqrt(1/(len(data1)-3)*tot)*3600.
            bestchi2=tot
            centerRA_best=centerRA
            centerDec_best=centerDec
            RA_eq_best=RA_eq
            Dec_eq_best=Dec_eq
            coef1best=[coef1[0],coef1[1],coef1[2]]
            coef2best=[coef2[0],coef2[1],coef2[2]]
        elif tot<bestchi2:
            bestchi2=tot
            centerRA_best=centerRA
            centerDec_best=centerDec
            RA_eq_best=RA_eq
            Dec_eq_best=Dec_eq
            coef1best=[coef1[0],coef1[1],coef1[2]]
            coef2best=[coef2[0],coef2[1],coef2[2]]
    coef1=[coef1best[0],coef1best[1],coef1best[2]]
    coef2=[coef2best[0],coef2best[1],coef2best[2]]
    centerRA=centerRA_best
    centerDec=centerDec_best
    RA_eq=RA_eq_best
    Dec_eq=Dec_eq_best
    
    eq1_str= "xi="+str(coef1[0])+"x + "+str(coef1[1])+"y + "+str(coef1[2])
    eq2_str= "eta="+str(coef2[0])+"x + "+str(coef2[1])+"y + "+str(coef2[2])
    print eq1_str
    print eq2_str
    fout=open(filename[:-4]+"_residuals.txt",'w')
    fout.write("Star Residuals: \n")
    fout.write("Julian Date: "+str(JulianDate([year,month,day,hour,minute,second])+time_diff/(24.*3600.))+"\n")
    fout.write("Local Sidereal Time: "+str(LST([year,month,day,hour,minute,second+time_diff]))+"\n")
    fout.write(eq1_str+"\n")
    fout.write(eq2_str+"\n")
    fout.write("Center RA: "+str(realRADec(deg_to_time(centerRA),deg_to_dms(centerDec),observation_LST)[1])+"\n")
    fout.write("Center Dec: "+str(realRADec(deg_to_time(centerRA),deg_to_dms(centerDec),observation_LST)[0])+"\n")
    fout.write("\n")
    fout.write("Angular Residual \n")
    ra_chi2=0
    dec_chi2=0
    dec_res_list=[]
    ra_res_list=[]
    for i in range(len(data1)):
        temp=realRADec(deg_to_time(data1[i][2]),deg_to_dms(data2[i][2]),observation_LST)
        RA_P=RA_eq(data1[i][0],data1[i][1])
        RA_A=temp[1]
        Dec_P=Dec_eq(data2[i][0],data2[i][1])
        Dec_A=temp[0]
        dec_res=abs(dms_to_deg(Dec_A)-dms_to_deg(Dec_P))
        dec_res_list.append(dec_res)
        dec_chi2+=abs(dec_res)**2
        ra_res=abs(time_to_deg(RA_A)-time_to_deg(RA_P))
        ra_res_list.append(ra_res)
        ra_chi2+=abs(ra_res)**2
        fout.write(names[i]+" | "+str(deg_to_time(ra_res))+" | "+str(deg_to_dms(dec_res))+"\n")
    fout.write("\n")
    fout.write("Maximal RA Residual: "+str(deg_to_time(max(ra_res_list))))
    fout.write("\n")
    fout.write("Maximal Dec Residual: "+str(deg_to_dms(max(dec_res_list))))
    fout.write("\n")
    ra_std_dev=sqrt(1/(len(data1)-3)*ra_chi2)
    dec_std_dev=sqrt(1/(len(data1)-3)*dec_chi2)
    fout.write("RA Standard Deviation: "+str(deg_to_time(ra_std_dev))+"\n")
    fout.write("Dec Standard Deviation: "+str(deg_to_dms(dec_std_dev)))
    return RA_eq,Dec_eq

image="H:/SSP/LSPR Images/July8Lick3.fit"
starfile="H:/SSP/LSPR Images/July8Lick3.txt"

def center(obj_img,x_c,y_c,n1=15,n2=15):
    obj_img[where(obj_img<0)]=0
    center=centroid(obj_img[y_c-1-n2:y_c-1+n2+1,x_c-1-n1:x_c-1+n1+1])
    return [[x_c+center[0][0],center[0][1]],[y_c+center[1][0],center[1][1]]]
def mag(v):
    ans=0
    for i in v:
        ans+=i**2
    return sqrt(ans)
def dot_product(v,w):
    ans=0
    for i in range(len(v)):
        ans+=v[i]*w[i]
    return ans

def LaTexify(filename):
    '''
        Turns a file in dlr star format into LaTex.
    '''
    star_file=open(filename)
    resfile=open(filename[:-4]+"_residuals.txt")
    comments=[]
    lines=[]
    linen=1
    residual_line=9
    res_lines=[]
    for line in resfile:
        res_lines.append(line)
    for line in star_file:
        #print linen
        linen+=1
        if line[0]=='%':
            i=0
            while line[i]=='%':
                i+=1
            comments.append(line[i:])
        elif line[0]==' ' or line[0]=='':
            continue
        else:
            parts=line.split(' | ')
            starname=parts[0]

            res_file_line=res_lines[residual_line].split(' | ')
            if (starname!=res_file_line[0]):
                print "Corrupted residual or star file. Please correct."
                return

            RAerror=res_file_line[1]
            RA_error_parts=RAerror[1:-2].split(', ')
            RA_errorh=float(RA_error_parts[0])
            RA_errorm=float(RA_error_parts[1])
            RA_errors=float(RA_error_parts[2])

            RA_errors=RA_errorh*3600.+RA_errorm*60.+RA_errors

            Decerror=res_file_line[2]
            Dec_error_parts=Decerror[1:-2].split(', ')
            Dec_errorh=float(Dec_error_parts[0])
            Dec_errorm=float(Dec_error_parts[1])
            Dec_errors=float(Dec_error_parts[2])
            Dec_errors=Dec_errorh*3600.+Dec_errorm*60.+Dec_errors
            
            
            starx=parts[1]
            stary=parts[2]
            RAh=parts[3]
            RAm=parts[4]
            RAs=parts[5]
            Decd=parts[6]
            Decm=parts[7]
            Decs=parts[8]

            RA_error_str=""
            Dec_error_str=""
            if RA_errors<0.001:
                RA_error_str=str("%.2e"%RA_errors)
                list1=RA_error_str.split('e')
                RA_error_str=list1[0]+"\cdot "+"10^{"+list1[1]+"}$"
            else:
                RA_error_str=str(float("%.2e"%RA_errors))+"$"
                
            if Dec_errors<0.001:
                Dec_error_str=str("%.2e"%Dec_errors)
                list1=Dec_error_str.split('e')
                Dec_error_str=list1[0]+"\cdot "+"10^{"+list1[1]+"}$"
            else:
                Dec_error_str=str(float("%.2e"%Dec_errors))+"$"          
            
            LaTexLine=starname+" & "+starx+" & "+stary+" & "
            LaTexLine=LaTexLine+RAh+"h "+RAm+"m "+RAs+" "
            LaTexLine=LaTexLine+"$\pm "+RA_error_str+"s & "
            LaTexLine=LaTexLine+Decd+"\\degrees \\space "+Decm+"' "+Decs[:-1]
            LaTexLine=LaTexLine+"$\pm "+Dec_error_str+"''"
            LaTexLine=LaTexLine+" \\\\ \\hline"
            lines.append(LaTexLine)
            residual_line+=1
    for line in lines:
        print line
