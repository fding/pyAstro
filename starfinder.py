import random, math
import numpy as np
import pyfits

CUTOFFLIMIT = 6

def gaussian(zscore,variance):
    return 0.5/math.pi/math.sqrt(variance) * math.exp(-0.5*zscore)

def distance(v1, v2):
    return sum([(v1[i]-v2[i])**2 for i in range(2)])

def mmax( x, y):
    return max(int(x),int(y))
    
def mmin( x, y):
  return min(int(x),int(y)) 

def findstars(image):
    count = 0
    stars = []
    for star in quick_findstars(image):
        if count>=50:
            break
        count+=1
        stars.append(star)
    delete_stars = []
    for i in range(len(stars)):
        if i in delete_stars:
            continue
        for j in range(i+1,len(stars)):
            if j in delete_stars:
                continue
            if distance(stars[i][0], stars[j][0])<2*max(stars[i][1],stars[j][1]):
                delete_stars.append(j)

    return_list = []
    for i in range(len(stars)):
        if i not in delete_stars:
            return_list.append(stars[i])
    return return_list


def quick_findstars(image): #np.ndarray[float, ndim=2] image):
    sizex=image.shape[1]
    sizey=image.shape[0]
    
    # Remove noise
    sq_list=[]
    for i in range(2000):
        sq_list.append(image[int(sizey*random.random())][int(sizex*random.random())])
    noise_level=sorted(sq_list)[int(9.4*len(sq_list)/10)]

    for y in range(sizey):
        for x in range(sizex):
            image[y][x]-= noise_level
            if image[y][x]<0:
                image[y][x]=0
    print noise_level
              
    stars = [] #[mean, variance, brightness]
    brightness_table = [[0.0 for x in range(sizex)] for y in range(sizey)]
    i=0
    while True:
        startpoint = (0,0)
        brightestvalue = 0
        # Find the brightest point in the image after removing previously found stars; this should be a new star!
        if i>0:
            for y in range(0,sizey,3):
                for x in range(0,sizex,3):
                    brightness_table[y][x]+=stars[i-1][2]*gaussian(distance((x,y),stars[i-1][0])/stars[i-1][1], stars[i-1][1])
        for y in range(0,sizey,3):
            for x in range(0,sizex,3):
                skip=False
                for star in stars:
                    if distance(star[0],(x,y))<CUTOFFLIMIT*star[1]:
                        skip=True
                        break
                if skip:
                    continue
                if image[y][x]-brightness_table[y][x]>brightestvalue:
                    startpoint=(x,y)
                    brightestvalue=image[y][x]-brightness_table[y][x]
        # No new stars are found
        if brightestvalue<=40:
            break
        
        # Maximum-likelihood gaussian that fits the star
        curmean = startpoint
        #print startpoint
        #print image[120][51]-brightness_table[120][51]
        # A probably-conservative estimate; but if true variance is larger, this procedure will converge
        curvariance = math.sqrt(brightestvalue)/30
        for j in range(2):
            # We want to limit computation to a square around the star
            xmin=mmax(0, curmean[0]-CUTOFFLIMIT*math.sqrt(curvariance))
            xmax=mmin(sizex,curmean[0]+CUTOFFLIMIT*math.sqrt(curvariance))
            ymin=mmax(0, curmean[1]-CUTOFFLIMIT*math.sqrt(curvariance))
            ymax=mmin(sizey,curmean[1]+CUTOFFLIMIT*math.sqrt(curvariance))
            curmean = [0,0]
            curvariance = 0
            curbrightness = 0
            for y in range(ymin,ymax):
                for x in range(xmin,xmax):
                    #cumulative_value = sum([star[2]*gaussian(distance((x,y),star[0])/star[1], star[1]) for star in stars])
                    curbrightness+=(image[y][x])            
            for y in range(ymin,ymax):
                for x in range(xmin,xmax):
                    #cumulative_value = sum([star[2]*gaussian(distance((x,y),star[0])/star[1], star[1]) for star in stars])
                    curmean[0]+=(image[y][x])/float(curbrightness)*x
                    curmean[1]+=(image[y][x])/float(curbrightness)*y
            
            for y in range(ymin,ymax):
                for x in range(xmin,xmax):
                    #cumulative_value = sum([star[2]*gaussian(distance((x,y),star[0])/star[1], star[1]) for star in stars])
                    curvariance+=(x-curmean[0])**2 *(image[y][x])/curbrightness
        stars.append([curmean,curvariance,curbrightness])
        i+=1
        yield stars[i-1]

                
            
