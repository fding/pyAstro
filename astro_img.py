#Library for image processing
from __future__ import division
from numpy import *
import pyfits
import numdisplay
from random import *
from numpy.random import *
import numpy.linalg as linalg

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

def centroid1(a):
    y_len=len(a)
    x_len=len(a[0])
    x_start=0
    y_start=len(a)
    x=0
    y=0
    s=0
    for i in range(x_len):
        for j in range(y_len):
            s+=a[j,i]
            y+=a[j,i]*(y_start-j)
            x+=a[j,i]*(x_start+i)
    x=x/s
    y=y/s
    return [x,y,s]

a=array([[0,33,21,33,8],[0,56,51,53,26],[23,120,149,73,18],[55,101,116,50,16],[11,78,26,2,10]])

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

#Finds background noise by taking a random sample and taking the median.
def background_noise(obj,n=45):
    sq_list=[]
    i=0
    while i<n:
        x=int(random()*obj.shape[1])
        y=int(random()*obj.shape[0])
        sq_list.append(obj[y,x])
        i+=1
    return sort(sq_list)[floor(len(sq_list)/2)]  

#rebins image from 1 by 1 to factor by factor
def rebin(img_array,factor=4):
    new_img = zeros((1+floor(img_array.shape[0]/factor),1+floor(img_array.shape[1]/factor)))
    new_x=0
    alternate_x=0
    for x in range(0,img_array.shape[1]):
        new_y = 0
        alternate_y = 0
        for y in range(0,img_array.shape[0]):
            #new_x = floor(x/factor)
            #new_y = floor(y/factor)
            new_img[new_y,new_x] += img_array[y,x]
            if alternate_y==factor:
                new_y +=1
                alternate_y = 0
            else:
                alternate_y +=1
        if alternate_x==factor:
            new_x+=1
            alternate_x=0
        else:
            alternate_x+=1
    return new_img

#remove noise by taking average noise of smaller quadrants
def remove_noise(img_array,num_sq=6):
    ymax=img_array.shape[0]
    xmax=img_array.shape[1]
    noise_array=[]
    for i in range(num_sq):
        row=[]
        for j in range(num_sq):
            square=img_array[floor(ymax/num_sq*j):floor((ymax)/num_sq*(j+1)),floor(xmax/num_sq*i):floor((xmax)/num_sq*(i+1))]
            row.append(background_noise(square))
        noise_array.append(row)
    coef_list=[]
    for i in range(0,num_sq-1,2):
        for j in range(0,num_sq-1,2):
            x_b=floor(ymax/num_sq*i)
            y_b=floor(ymax/num_sq*j)
            x1=floor(xmax/num_sq/2)
            x2=floor(xmax/num_sq/2*3)
            y1=floor(ymax/num_sq/2)
            y2=floor(ymax/num_sq/2*3)
            coef=bilinear_interpolation(array([[x1,y1,noise_array[j][i]],[x2,y1,noise_array[j][i+1]],[x1,y2,noise_array[j+1][i]],[x2,y2,noise_array[j+1][i+1]]]))
            square=img_array[floor(ymax/num_sq*j):floor((ymax)/num_sq*(j+2)),floor(xmax/num_sq*i):floor((xmax)/num_sq*(i+2))]
            for x_i in range(square.shape[1]):
                for y_i in range(square.shape[0]):
                    img_array[y_b+y_i,x_b+x_i]=img_array[y_b+y_i,x_b+x_i]-coef[0]*x_i-coef[1]*y_i-coef[2]
    
    return img_array
#computes interquartile range
def iQr(image):
    data=[]
    for row in image:
        for pt in row:
            data.append(pt)
    data=sort(data)
    n=len(data)
    return data[floor(0.75*n)]-data[floor(0.25*n)]
#computes if there are stars in image by testing for randomness (specifically, clumping)
def unrandomness(image):
    image=200*image/image.mean()
    
    y_list=[]
    for row in image:
        y_list.append(sum(row))
    large_indices=where(y_list>array(y_list).mean())[0]
    r1= abs(100*(large_indices.mean()-len(y_list)/2)/len(y_list)) #random_index1 should be 0 for random sample
    r2=abs(100*(large_indices[:len(large_indices)/2].mean()-len(y_list)*0.25)/len(y_list)) #random_index2 should be 0 for random sample
    r3=abs(100*(large_indices[len(large_indices)/2:].mean()-len(y_list)*0.75)/len(y_list)) #random_index3 should be 0 for random sample
    r123=max(r1,r2,r3)
    #r1,r2, and r3 collectively measures clumping in the y direction
    x_list=[]
    for i in range(image.shape[1]):
        x_list.append(sum(image[:,i]))
    large_indices=where(x_list>array(x_list).mean())[0]
    r4= abs(100*(large_indices.mean()-len(x_list)/2)/len(x_list)) #random_index1 should be 0 for random sample
    r5=abs(100*(large_indices[:len(large_indices)/2].mean()-len(x_list)*0.25)/len(x_list)) #random_index2 should be 0 for random sample
    r6=abs(100*(large_indices[len(large_indices)/2:].mean()-len(x_list)*0.75)/len(x_list)) #random_index3 should be 0 for random sample
    r456=max(r4,r5,r6)
    #r4,r5,r6 collectively measures clumping in the x direction
    unrandomness=(r123**2+r456**2)*(image.shape[0]+image.shape[1])
    return unrandomness


Table=[]

def sample(n,m):
    return_val=[]
    
    for i in range(30):
        random_image=rand(n,m)
        #print "Sample", i
        #print "Size: ", n
        #print "Mean: ", random_image.mean()
        #print "Metric: ", contains_signal(random_image)
        return_val.append(unrandomness(random_image))
    return sort(return_val)

def init_Table():
    print "Initializing Library. Please wait..."
    for i in range(2,52,2):
        row=[]
        for j in range(1,50,2):
            statistic=sort(sample(i,j))
            row.append(statistic[int(0.9*29)])
        Table.append(row)
    print "Done"

#init_Table()

def TableLookup(i,j):
    i=int(i/2)
    j=int(j/2)
    if (i>=25):
        if (j>=25):
            return Table[24][24]
        return Table[24][j]
    if (j>=25):
        return Table[i][24]
    return Table[i][j]
            
def contains_signal(image):
    noise_range=iQr(image)
    mean_sig=image.mean()
    img_max=image.max()
    if img_max/(noise_range+1)>4:
        return True
    image=image/image.mean()
    if img_max<0.5:
        return False
    
    #if unrandomness(image)>TableLookup(image.shape[0],image.shape[1]):
    if unrandomness(image)>30000:
        return True
    return False

def flatten(obj_img,flat_img):
    numdisplay.display(obj_img,z1=0,z2=400)
    flat_img=flat_img/flat_img.mean()
    obj_img=obj_img/flat_img
    return obj_img

def center(obj,flat,x_c,y_c,n1=13,n2=13):
    obj_img[where(obj_img<0)]=0
    center=centroid(obj_img[y_c-1-n2:y_c-1+n2+1,x_c-1-n1:x_c-1+n1+1])
    return [[x_c+center[0][0],center[0][1]],[y_c+center[1][0],center[1][1]],[s,u_s]]

def process1(obj_img,num_times=5):
    obj_img=obj_img-background_noise(obj_img,300,300)
    box_list=[] #contain list of top left and bottom right points of each box [x1,y1, x2, y2], not including x2 and y2
    box_list.append([[0,0,obj_img.shape[1]-1,obj_img.shape[0]-1]])
    star_size=24
    x_size=8
    y_size=10
    for i in range(num_times):
        temp_list=[]
        for box in box_list[i]:
            x1=floor((box[2]-box[0])/2.0)
            y1=floor((box[3]-box[1])/2.0)
            subimg1=obj_img[box[1]:box[1]+y1,box[0]:box[0]+x1]
            subimg2=obj_img[box[1]:box[1]+y1,box[0]+x1:box[2]]
            subimg3=obj_img[box[1]+y1:box[3],box[0]:box[0]+x1]
            subimg4=obj_img[box[1]+y1:box[3],box[0]+x1:box[2]]
            if contains_signal(subimg1): #if there is a star in subimg1
                temp_list.append([box[0],box[1],box[0]+x1,box[1]+y1])

            if contains_signal(subimg2):
                temp_list.append([box[0]+x1,box[1],box[2],box[1]+y1])

            if contains_signal(subimg3):
                temp_list.append([box[0],box[1]+y1,box[0]+x1,box[3]])

            if contains_signal(subimg4):
                temp_list.append([box[0]+x1,box[1]+y1,box[2],box[3]])
            box_list.append(temp_list)
    stars=[]
    print "Part 2"
    for box in box_list[num_times-1]:
        pt=centroid1(obj_img[box[1]:box[3],box[0]:box[2]])
        x=pt[0]+box[0]+1
        y=pt[1]+box[1]+1
        stars.append([x,y,pt[2]])
    return stars

def display(brt_pt, size, radius=1):
    img=zeros(size)
    for pt in brt_pt:
        x=round(pt[0])
        y=round(pt[1])
        x_max=min(x+radius,size[1])
        x_min=max(x-radius,1)
        y_max=min(y+radius,size[0])
        y_min=max(y-radius,1)
        img[y_min-1:y_max,x_min-1:x_max]=100
    return img
                
test1=True
if test1:
    img=pyfits.getdata("C:/Python27/series2_master_2.fit")
    #numdisplay.display(image)
    img1=img[378:417,393:436] #contains signal
    img2=img[526:554,236:266] #contains signal
    img3=img[355:385,753:800] #contains signal
    img4=img[295:324,512:536] #contains signal
    img5=img[276:313,328:373] #contains signal
    img9=img[259:322,234:257]

    img6=img[388:411,529:577] #contains no signal
    img7=img[254:272,429:457] #contains no signal
    img8=img[265:284,796:813] #contains no signal
    img10=img[560:620,1000:1050]

    print contains_signal(img1)
    print contains_signal(img2)
    print contains_signal(img3)
    print contains_signal(img4)
    print contains_signal(img5)
    print contains_signal(img9)

    print contains_signal(img6)
    print contains_signal(img7)
    print contains_signal(img8)

test2=False
#image=pyfits.getdata("C:/Python27/series2_master_2.fit")
if test2:
    image=pyfits.getdata("C:/Python27/series2_master_2.fit")
    img1=rebin(image[423:441,377:392]) #contains signal
    img2=rebin(image[526:554,236:266]) #contains signal
    img3=rebin(image[385:389,753:800]) #contains signal
    img4=rebin(image[295:324,512:536]) #contains signal
    img5=rebin(image[304:387,247:327]) #contains signal

    img6=rebin(image[388:411,529:577]) #contains no signal
    img7=rebin(image[254:272,429:457]) #contains no signal
    img8=rebin(image[265:284,796:813]) #contains no signal

    contains_signal(img1)
    contains_signal(img2)
    contains_signal(img3)
    contains_signal(img4)
    contains_signal(img5)

    contains_signal(img6)
    contains_signal(img7)
    contains_signal(img8)
