import numpy as np
import math
#a = np.array([5.576091766357422, 8.455062866210938, -0.06940200179815292])
#b = np.array([6.994688034057617, 9.038549423217773, -0.23745399713516235])
#c = np.array([7.181145191192627, 10.318812370300293, -0.020966999232769012])
#d = np.array([8.561141967773438, 10.942671775817871, -0.06447900086641312])
def getNormedVector(a,b):
    return (b-a)/np.linalg.norm(b-a)

def getAngle(a,b):
    return np.rad2deg(np.arccos(np.dot(a/np.linalg.norm(a),b/np.linalg.norm(b))))

def getSign(a,b):
    b = b/np.linalg.norm(b) 
    rec = np.dot(a,b)
    return rec

def getdih(a,b,c,d):
    #v1 = np.subtract(a,b)
    #v2 = np.subtract(b,c)
    #v3 = np.subtract(c,d)
    v1 = getNormedVector(a,b)
    v2 = getNormedVector(b,c)
    v3 = getNormedVector(c,d)
    v1v2 = np.cross(v1,v2)
    v2v3 = np.cross(v2,v3)
    d12 = np.dot(v1v2,v1v2)
    d23 = np.dot(v2v3,v2v3)
    rotdir = getSign(v1,-v2v3)
    #print d12, d23
    #print "sign", rotdir
    rec = getAngle(v1v2,v2v3)
    #print rec
    # if the rotdir is positive, the rotation is clockwise
    if rotdir > 0.0:
        #rec = -rec
        rec = 360.0 - rec
    #costheta = float(np.dot(d12,d23))/float((d12*d23))
    #theta = math.acos(costheta)/math.pi*180.0
    #print costheta, theta
    #return theta
    #if rec < 0.0:
    #    rec = 180.0 + rec
    return rec
#theta = getdih(a,b,c,d)
#print theta


