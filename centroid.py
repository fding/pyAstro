#Centroid program
from __future__ import division
from numpy import *
import pyfits
import numdisplay
from random import *

#Problem 2 b
#The file should have each row seperated by rows in files, and columns seperated by spaces
def P2b(f):
    a=open(f)
    str_input1=a.readline().split(' ')
    str_input2=a.readline().split(' ')
    str_input3=a.readline().split(' ')
    int_array=[]
    row=[]
    for i in str_input1:
        row+=[int(i)]
    int_array.append(row)
    row=[]
    for i in str_input2:
        row+=[int(i)]
    int_array.append(row)
    row=[]
    for i in str_input3:
        row+=[int(i)]
    int_array.append(row)
    return centroid(int_array)
        
#computes centroid, returns [[x, uncertainty],[y,uncertainty]]
#uncertainty is computed by computing the uncertainty of each element in array (equal to squareroot), and doing the appropriate arithmetic on uncertainties.
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

a=array([[0,33,21,33,8],[0,56,51,53,26],[23,120,149,73,18],[55,101,116,50,16],[11,78,26,2,10]])

#Finds background noise by taking a random sample and taking the median.
def background_noise(obj,x,y):
    r=40
    n=70
    sq_list=[]
    for i in range(n):
        pt=obj[y+floor(2*r*(random()-0.5)),x+floor(2*r*(random()-0.5))]
        sq_list.append(pt)
    return sort(sq_list)[floor(len(sq_list)/2)]  
            

def process(obj,flat):
    #Normalize with respect to flat file
    obj_img=pyfits.getdata(obj)
    numdisplay.display(obj_img,z1=2000,z2=40000)
    flat_img=pyfits.getdata(flat)
    flat_img=flat_img/flat_img.mean()
    obj_img=obj_img/flat_img

    x_c=input("Enter x coordinate of asteroid: ")
    y_c=input("Enter y coordinate of asteroid: ")
    #Subtract background noise
    obj_img=obj_img-background_noise(obj_img,x_c,y_c)
    obj_img[where(obj_img<0)]=0
    numdisplay.display(obj_img,z1=100,z2=400)
    #Computes centroid
    n=input("Radius of asteroid-field: ")
    center=centroid(obj_img[y_c-1-n:y_c-1+n+1,x_c-1-n:x_c-1+n+1])
    obj_img*=0.1
    x=round(x_c+center[0][0])
    y=round(y_c+center[1][0])
    obj_img[y-3:y+3,x-3:x+3]=20000
    numdisplay.display(obj_img,z1=100,z2=10000)
    return [[x_c+center[0][0],center[0][1]],[y_c+center[1][0],center[1][1]]]
    
