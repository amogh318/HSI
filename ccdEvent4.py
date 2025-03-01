import math
import numpy as np


class ccdPixelHit():
    def __init__(self, tx, ty):
        self.xpixel = tx
        self.ypixel = ty
        self.edep = 0.0
        self.iRed = 0.0
        self.iBlue = 0.0
        self.iGreen = 0.0

    def GetXPixel(self):
        return self.xpixel

    def GetYPixel(self):
        return self.ypixel

    def GetEdep(self):
        return self.edep

    def SetEdep(self, xx):
        self.edep = xx

    def isSame(self, xx):
        if (xx.GetXPixel() == self.xpixel and xx.GetYPixel() == self.ypixel):
            return True
        else:
            return False


class ccdHoughCell():
    def __init__(self, pix1, pix2):
        self.pixel1 = pix1
        self.pixel2 = pix2
        self.delX = float(pix2.GetXPixel())-float(pix1.GetXPixel())
        self.delY = float(pix2.GetYPixel())-float(pix1.GetYPixel())
        self.slope = (self.delY)/(self.delX+1e-20)
        self.intercept = pix2.GetYPixel()-(self.slope*pix1.GetXPixel())    
        self.theta = math.atan2(-self.delX,self.delY)*180/math.pi
        self.rho = math.fabs(self.intercept)*math.fabs(math.sin(self.theta*math.pi/180))

    def GetPixel1(self):
        return self.pixel1

    def GetPixel2(self):
        return self.pixel2
    
    def GetdelX(self):
        return self.delX
                                          
    def GetdelY(self):
        return self.delY                                   
             
    def GetSlope(self):
        return self.slope

    def GetIntercept(self):
        return self.intercept

    def GetTheta(self):
        return self.theta

    def GetRho(self):
        return self.rho


class ccdEvent():
    def __init__(self, tfname):
        self.eventID = -1
        self.momin = 0
        self.thein = -1000
        self.phiin = -1000
        self.allHits = []
        self.nhits = self.fileRead(tfname)
        self.meanRho = -10.0
        self.sigmaRho = -10.0
        self.meanTheta = -1000.0
        self.sigmaTheta = -1000.0
        self.houghSelHits = []
        #print("evt1:",self.eventID)
        if(self.nhits>1):
            self.houghSelHits = self.ApplyHoughTransform()

    def GetEventID(self):
        return self.eventID
        
    def GetThetaTrue(self):
        return self.thein
    
    def GetPhiTrue(self):
        return self.phiin
            
    def GetMeanRho(self):
        return self.meanRho

    def GetSigmaRho(self):
        return self.sigmaRho

    def GetMeanTheta(self):
        return self.meanTheta

    def GetSigmaTheta(self):
        return self.sigmaTheta

    def AddPixelHit(self, xx):
        self.allHits.append(xx)

    def GetNPixelHits(self):
        return len(self.allHits)

    def GetNHoughPixelHits(self):
        return len(self.houghSelHits)

    def AddHoughSelHit(self, xx):
        isFill = True
        for ihc in self.houghSelHits:
            if (ihc.isSame(xx)):
                isFill = False
                break
        if isFill:
            self.houghSelHits.append(xx)

    def GetPixelHit(self, xx):
        return self.allHits[xx]

    def GetHoughSelHit(self, xx):
        return self.houghSelHits[xx]

    def fileRead(self, fileName):
        inputFile = open(fileName, 'r')
        ilin = 0
        thits = 0
        for xx in inputFile:
            yy = xx.split()
            if (ilin == 0):
                self.eventID = int(yy[0])
                thits = int(yy[2])
            elif (ilin == 1):
                self.pidin = int(yy[1])
                self.momin = float(yy[2])
                self.thein = float(yy[3])
                self.phiin = float(yy[4])
            else:
                ipHit = ccdPixelHit(int(yy[2]), int(yy[3]))
                ipHit.SetEdep(float(yy[4]))
                self.allHits.append(ipHit)
            ilin += 1

        # print(len(iHitList))
        return thits

    def ApplyHoughTransform(self):
        npix = self.GetNPixelHits()
        houghList = []
        thetaList = []
        rhoList = []

        for ij in range(npix-1):
            for jk in range(ij+1, npix):
                iHoughCell = ccdHoughCell(self.GetPixelHit(ij), self.GetPixelHit(jk))
                houghList.append(iHoughCell)
                thetaList.append(iHoughCell.GetTheta())
                rhoList.append(iHoughCell.GetRho())
        if(len(thetaList)>5):
            self.meanRho = np.mean(rhoList)
            self.sigmaRho = math.sqrt(np.var(rhoList))
            self.meanTheta = np.mean(thetaList)
            self.sigmaTheta = math.sqrt(np.var(thetaList))
        else:
            self.meanRho = -10.0
            self.sigmaRho = -10.0
            self.meanTheta = -1000.0
            self.sigmaTheta = -1000.0
            
        for iter in range(100):
            listRho = []
            listTheta = []
            newHoughList = []
            nHoughCell = len(houghList)
            for ihc in houghList:
                if (ihc.GetTheta() < (self.meanTheta+self.sigmaTheta) and ihc.GetTheta() > (self.meanTheta-self.sigmaTheta)
                        and ihc.GetRho() < (self.meanRho + self.sigmaRho) and ihc.GetRho() > (self.meanRho - self.sigmaRho)):
                    newHoughList.append(ihc)
                    listRho.append(ihc.GetRho())
                    listTheta.append(ihc.GetTheta())
            #print("mean, sigma : ", iter)
            #print(self.meanRho, self.sigmaRho)
            #print(self.meanTheta, self.sigmaTheta)
            if (len(newHoughList) < nHoughCell):
                nHoughCell = len(newHoughList)
                if(len(thetaList)>5):
                    self.meanRho = np.mean(rhoList)
                    self.sigmaRho = math.sqrt(np.var(rhoList))
                    self.meanTheta = np.mean(thetaList)
                    self.sigmaTheta = math.sqrt(np.var(thetaList))
                else:
                    self.meanRho = -10.0
                    self.sigmaRho = -10.0
                    self.meanTheta = -1000.0
                    self.sigmaTheta = -1000.0
                houghList = newHoughList
                if (math.fabs(self.sigmaTheta) < 0.0001):
                    self.sigmaTheta = 0.0001
                if (math.fabs(self.sigmaRho) < 0.0001):
                    self.sigmaRho = 0.0001
            else:
                break

        return houghList, rhoList, thetaList

