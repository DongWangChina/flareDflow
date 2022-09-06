from cmath import exp, sqrt
from multiprocessing import Pool
import numpy as np
from scipy.stats import beta
import tarfile
# import os
# from time import time
# from scipy import interpolate

def PDF_Intergrate(cbDict,solFln):
    global theta, z_int, c_int, gz_int, gc_int, gcz_int, solidx
    global z_space, c_space, Src_vals, Yi_vals
    global nScalars, nYis, n_points_z, n_points_c
    global int_pts_z, int_pts_c, int_pts_gz, int_pts_gc, int_pts_gcz

    for solidx in range (cbDict["n_points_h"]):

        z_int = cbDict["z"]
        c_int = cbDict["c"]
        gz_int = cbDict["gz"]
        gc_int = cbDict["gc"]
        gcz_int = cbDict["gcz"]

        nScalars = int(cbDict["nScalars"])
        nYis = int(cbDict["nYis"])

        n_points_z = int(cbDict["n_points_z"])
        n_points_c = int(cbDict["n_points_c"])

        int_pts_z = int(cbDict["int_pts_z"])
        int_pts_c = int(cbDict["int_pts_c"])
        int_pts_gz = int(cbDict["int_pts_gz"])
        int_pts_gc = int(cbDict["int_pts_gc"])
        int_pts_gcz = int(cbDict["int_pts_gcz"])

        theta = 1.0

        tar = tarfile.open(solFln,'r')
        tar.extractall()
        tar.close()

        # read chemTab_01.date
        file_path = 'chemTab_' + '%02d' % (solidx+1) + '.dat'
        z_space = np.zeros(n_points_z)
        c_space = np.zeros(n_points_c)
        Src_vals = np.zeros((n_points_z, n_points_c, nScalars))
        Yi_vals = np.zeros((n_points_z, n_points_c, nYis))
        with open (file_path) as file:
            for i in range(n_points_z):
                for j in range(n_points_c):
                    intp_value = (file.readline()).split()
                    z_space[i] = intp_value[0]
                    c_space[j] = intp_value[1]
                    Src_vals[i,j,:] = intp_value[2 : (nScalars+1+1)]
                    Yi_vals[i,j,:] = intp_value[(nScalars+1+1) :]
                    
        # final integration
        # with Pool(processes=cbDict['n_procs']) as pool:
        #     pool.map(unit_Compute,range(int_pts_z))
        with Pool() as pool:
            pool.map(unit_Compute,range(int_pts_z))
        
        # with Pool(processes=1) as pool:
        #     pool.map(unit_Compute,range(int_pts_z))

        tarPath = 'pdfData.tar'
        with tarfile.open(tarPath, 'a') as tar:
            tar.add('./' + cbDict['output_fln'])
            tar.add('./d2Yeq_table.dat')
            for iz in range(int_pts_z):
                tar.add("./unit" +  '%02d' % (iz+1) + '_h' + '%02d' % (solidx+1) + '.dat')

    return tarPath

def unit_Compute(iz):
    # print("run Unit: "+str(iz))

    unitFile = "unit" +  '%02d' % (iz+1) + '_h' + '%02d' % (solidx+1) + '.dat'
    unit_data = []
    z_loc = locate(z_space,n_points_z,z_int[iz])
    for ic in range(int_pts_c):
        # c_loc = locate(c_space, n_points_c,c_int[ic])
        for igz in range(int_pts_gz):
            for igc in range(int_pts_gc):
                for igcz in range(int_pts_gcz):

                    # print("iz = "+str(iz)+"   ic = "+str(ic)+"  igz = "+str(igz)+\
                    #       "igc = "+str(igc)+"   igcz = "+str(igcz))
                    # print("igz = "+str(igz))
                    # print("igc = "+str(igc))
                    # print("igcz = "+str(igcz))

                    # allocate RAMspace for array
                    yint = np.zeros(nScalars)
                    Yi_int = np.zeros(nYis)

                    # z = 0 or 1
                    if (iz==0 or iz==int_pts_z-1):
                        yint[1] = 0.0
                        yint[2] = 0.0
                        yint[3] = 0.0

                        if (iz==0): z_loc = 0
                        if (iz==int_pts_z-1): z_loc = n_points_z-1

                        yint[4:nScalars] = Src_vals[z_loc,0,4:nScalars]
                        Yi_int[:] = Yi_vals[z_loc,0,:]
                    
                    # c = 0 or 1
                    elif (ic==0 or ic==int_pts_c-1):
                        yint[1] = 0.0
                        yint[2] = 0.0
                        yint[3] = 0.0

                        if (ic==0): c_loc = 0
                        if (ic==int_pts_c-1): c_loc = int_pts_c-1

                        yint[4:nScalars] = Src_vals[0,c_loc,4:nScalars]
                        Yi_int[:] = Yi_vals[0,c_loc,:]

                    # gZ=0 and gc=0
                    elif ((igz==0 or igz==int_pts_gz-1) and (igc==0 or igc==int_pts_gc-1)):
                        yint,Yi_int = delta(z_int[iz],c_int[ic],z_space,c_space,yint,\
                                            Src_vals,Yi_int,Yi_vals)

                    # gZ=0 and gc>0
                    elif (igz==0 or igz==int_pts_gz-1):
                        yint,Yi_int = cbeta(c_int[ic],gc_int[igc],c_space,z_int[iz], \
                                            z_space,yint,Src_vals,Yi_int,Yi_vals)

                    # gZ>0 and gc=0
                    # elif (igz>0 and igc==0):
                    elif (igc==0 or igc==int_pts_gc-1):
                        yint,Yi_int = zbeta(z_int[iz],gz_int[igz],z_space,c_int[ic], \
                                            c_space,yint,Src_vals,Yi_int,Yi_vals)

                    # gZ>0 and gc>0
                    else:
                        c_var = gc_int[igc]*(c_int[ic]*(1.0-c_int[ic]))
                        z_var = gz_int[igz]*(z_int[iz]*(1.0-z_int[iz]))
                        co_var = gcz_int[igcz]*sqrt(c_var)*sqrt(z_var)*0.98

                        yint,Yi_int = int_point(z_int[iz],c_int[ic],c_var,z_var,co_var, \
                                yint,z_space,c_space,Src_vals,Yi_int,Yi_vals,theta)
                    
                    unit_D = [z_int[iz],c_int[ic],gz_int[igz],gc_int[igc],gcz_int[igcz]]
                    for kk in range(1,nScalars):
                        unit_D.append(yint[kk])
                    for kk in range(nYis):
                        unit_D.append(Yi_int[kk])
                    unit_data.append(unit_D)

    np.savetxt(unitFile,unit_data,fmt='%.16E')
    print("Done run Unit: "+str(iz))


#--------------------------------------#
# -----------sub functions-------------#
#--------------------------------------#

#find the lower position of x in array xx(with length n)
def locate(xx,n,x):  
    j = 0
    if (x < xx[0]):
        pass
    elif (x > xx[n-1]):
        j = n-1
    else:
        for ii in range(0,n-2):
            if ((x>xx[ii]) & (x<xx[ii+1])):
                j = ii
    
    return j

# linear interpolation factor
def intfac(xx,xarray,loc_low):
    fac = 0.0
    if (xx<xarray[loc_low]):
        fac = 0.0
    elif (xx>xarray[loc_low+1]):
        fac = 1.0
    else:
        fac = (xx-xarray[loc_low]) / (xarray[loc_low+1] - xarray[loc_low])

    return fac


# delta function
def delta(z_int,c_int,z_space,c_space,y_int,Src_vals,Yi_int,Yi_vals):
    # print("delta")

    z_loc = locate(z_space,n_points_z,z_int)
    z_fac = intfac(z_int,z_space,z_loc)

    c_loc = locate(c_space,n_points_c,c_int)
    c_fac = intfac(c_int,c_space,c_loc)

    y_int[0:nScalars] = (1.0-c_fac) * ( z_fac*Src_vals[z_loc+1,c_loc,0:nScalars] \
                        + (1.0-z_fac) * Src_vals[z_loc,c_loc,0:nScalars] ) \
                    + c_fac * ( z_fac*Src_vals[z_loc+1,c_loc+1,0:nScalars] \
                        + (1.0-c_fac) * Src_vals[z_loc,c_loc+1,0:nScalars] )


    # omega_c/rho
    y_int[1] = y_int[1]/y_int[0]
    # gc_source
    y_int[2] = y_int[1]*c_int*y_int[nScalars-1]
    # gz_source
    y_int[3] = y_int[1]*z_int

    # scalars !number of integral scalars Yi
    Yi_int[0:nYis] = (1.0-c_fac) * (c_fac*Yi_vals[z_loc+1,c_loc,0:nYis] \
                    + (1.0-z_fac) * Yi_vals[z_loc,c_loc,0:nYis] ) \
                + c_fac * (z_fac*Src_vals[z_loc+1,c_loc+1,0:nYis] \
                    + (1.0-z_fac) * Src_vals[z_loc,c_loc+1,0:nYis])

    np.maximum( abs(Yi_int), 0.0, out=Yi_int )
    
    return y_int,Yi_int


# cbeta function
def cbeta(c_int,gc_int,c_space,z_int,z_space,y_int,Src_vals,Yi_int,Yi_vals):
    # print("cbeta")
    Src_vals_int = np.zeros((n_points_c, nScalars))
    Yi_vals_int = np.zeros((n_points_c, nYis))
    cdf001 = np.zeros(n_points_c)
    dsdc = np.zeros((nScalars, n_points_c+1))
    dYdc = np.zeros((nYis, n_points_c+1))
    loc = locate(z_space,n_points_z,z_int)
    fac = intfac(z_int,z_space,loc)

    Src_vals_int[0:n_points_c,0:nScalars] = fac*Src_vals[loc+1,0:n_points_c,0:nScalars] \
                            + ( (1.0-fac)*Src_vals[loc,0:n_points_c,0:nScalars] )
    Yi_vals_int[0:n_points_c,0:nYis] = fac*Yi_vals[loc+1,0:n_points_c,0:nYis] \
                                        + ( (1.0-fac)*Yi_vals[loc,0:n_points_c,0:nYis] )


    alpha_c = c_int*( (1.0/gc_int)-1.0 )
    beta_c = (1.0-c_int) * ( (1.0/gc_int)-1.0 )

    #===== integration by part =====#
    cdf001[1:n_points_c-1] = beta.cdf((c_space[1:(n_points_c-1)]+c_space[2:(n_points_c)])/2.0,
                                        alpha_c,beta_c)
    cdf001[n_points_c-1] = 1.0

    for j in range(nScalars):
        if (j==0 or j>=4):  #1:rho 5:cpe 6:mwt 7:sum_dhfi
            dsdc[j,0:n_points_c-1] = ( Src_vals_int[1:n_points_c,j]-Src_vals_int[0:n_points_c-1,j] ) \
                        / (c_space[1:n_points_c]-c_space[0:n_points_c-1])
        elif (j==1):  #2:omegac
            dsdc[j,0:n_points_c-1] = ( Src_vals_int[1:n_points_c,1]/Src_vals_int[1:n_points_c,0] \
                                        -Src_vals_int[0:n_points_c-1,1]/Src_vals_int[0:n_points_c-1,0] ) \
                                    /(c_space[1:n_points_c]-c_space[0:n_points_c-1])
        elif (j==2): #3:c*omega
            dsdc[j,0:n_points_c-1] = ( c_space[1:n_points_c]*Src_vals_int[1:n_points_c,nScalars-1] \
                                            *Src_vals_int[1:n_points_c,1]/Src_vals_int[1:n_points_c,0] \
                                        - c_space[0:n_points_c-1]*Src_vals_int[0:n_points_c-1,nScalars-1] \
                                            *Src_vals_int[0:n_points_c-1,1]/Src_vals_int[0:n_points_c-1,0] ) \
                                    / (c_space[1:n_points_c]-c_space[0:n_points_c-1])
        else: #4:Z*omegac
            dsdc[j,0:n_points_c-1] = ( z_int*Src_vals_int[1:n_points_c,1]/Src_vals_int[1:n_points_c,0] \
                                        -z_int*Src_vals_int[0:n_points_c-1,1]/Src_vals_int[0:n_points_c-1,0] ) \
                                    / (c_space[1:n_points_c]-c_space[0:n_points_c-1])

        
        dsdc[j,n_points_c-1] = dsdc[j,n_points_c-2]
        dsdc[j,n_points_c] = dsdc[j,0]

        #--------start integration---------#
        y_int[j] = y_int[j]-0.5*np.inner( ( dsdc[j,0:n_points_c-2]*cdf001[0:n_points_c-2] 
                                            + dsdc[j,1:n_points_c-1]*cdf001[1:n_points_c-1] ), \
                                        ( 0.5*(c_space[2:n_points_c]+c_space[1:n_points_c-1]) \
                                            - 0.5*(c_space[1:n_points_c-1]+c_space[0:n_points_c-2]) ) )

    y_int[0:nScalars] = y_int[0:nScalars]-dsdc[0:nScalars,n_points_c-1]*cdf001[n_points_c-1] \
                *( c_space[n_points_c-1]-c_space[n_points_c-2] )/2.0 \
            - dsdc[0:nScalars,n_points_c]*cdf001[0]*(c_space[1]-c_space[0])/2.0 \
            + Src_vals_int[n_points_c-1,0:nScalars]

    for j in range(nYis):
        #compute derivatives
        dYdc[j,0:n_points_c-1] = (Yi_vals_int[1:n_points_c,j]-Yi_vals_int[0:n_points_c-1,j] ) \
                                / (c_space[1:n_points_c]-c_space[0:n_points_c-1])

        dYdc[j,n_points_c-1] = dYdc[j,n_points_c-2]
        dYdc[j,n_points_c] = dYdc[j,0]

        #--------start integration---------#
        Yi_int[j] = Yi_int[j]-0.5*np.inner( ( dYdc[j,0:n_points_c-2]*cdf001[0:n_points_c-2] \
                                                +dYdc[j,1:n_points_c-1]*cdf001[1:n_points_c-1] ), \
                                            ( (c_space[2:n_points_c]+c_space[1:n_points_c-1])*0.5 \
                                                - (c_space[1:n_points_c-1]+c_space[0:n_points_c-2])*0.5 ) )

    Yi_int[0:nYis] = Yi_int[0:nYis]-dYdc[0:nYis,n_points_c-1]*cdf001[n_points_c-1] \
                    *(c_space[n_points_c-1]-c_space[n_points_c-2])/2.0 \
                - dYdc[0:nYis,n_points_c]*cdf001[0]*( c_space[1]-c_space[0] )/2.0 \
                + Yi_vals_int[n_points_c-1,0:nYis]

    np.maximum( abs(Yi_int), 0.0, out=Yi_int )
    
    return y_int,Yi_int


def zbeta(z_int,gz_int,z_space,c_int,c_space,y_int,Src_vals,Yi_int,Yi_vals):
    # print("zbeta")
    dYdc = np.zeros((nYis,n_points_z+1))
    cdf001 = np.zeros(n_points_z)
    dsdc = np.zeros((nScalars,n_points_z+1))
    Src_vals_int = np.zeros((n_points_z,nScalars))
    Yi_vals_int = np.zeros((n_points_z,nScalars))

    loc = locate(c_space,n_points_c,c_int)
    fac = intfac(c_int,c_space,loc)

    # interp sc_vals
    Src_vals_int[0:n_points_z,0:nScalars] = fac*Src_vals[0:n_points_z,loc+1,0:nScalars] \
                                            + (1.0-fac)*Src_vals[0:n_points_z,loc,0:nScalars]
    Yi_vals_int[0:n_points_z,0:nYis] = fac*Yi_vals[0:n_points_z,loc+1,0:nYis] \
                                        + (1.0-fac)*Yi_vals[0:n_points_z,loc,0:nYis]


    alpha_z = z_int*(1.0/gz_int-1.0)
    beta_z = (1.0-z_int)*(1.0/gz_int-1.0)

    #===== integration by part =====#
    cdf001[1:n_points_z-1] = beta.cdf((z_space[1:(n_points_z-1)]+z_space[2:(n_points_z)])/2.0,
                                        alpha_z,beta_z)
    cdf001[n_points_z-1] = 1.0

    for j in range(nScalars):
        if (j==0 or j>=4): #1:density 5:cpe 6:mwt 7:sum_dhfi
            dsdc[j,0:n_points_z-1] = ( Src_vals_int[1:n_points_z,j]-Src_vals_int[0:n_points_z-1,j] ) \
                                    /(z_space[1:n_points_z]-z_space[0:n_points_z-1])
        
        elif (j==1): #2:omegac
            dsdc[j,0:n_points_z-1] = ( Src_vals_int[1:n_points_z,1]/Src_vals_int[1:n_points_z,0] \
                                        -Src_vals_int[0:n_points_z-1,1]/Src_vals_int[0:n_points_z-1,0] ) \
                                    / (z_space[1:n_points_z]-z_space[0:n_points_z-1])

        elif (j==2): #3:c*omegac
            dsdc[j,0:n_points_z-1] = (c_int*Src_vals_int[1:n_points_z,nScalars-1] \
                                            *Src_vals_int[1:n_points_z,1]/Src_vals_int[1:n_points_z,0] \
                                        - c_int*Src_vals_int[0:n_points_z-1,nScalars-1] \
                                            *Src_vals_int[0:n_points_z-1,1]/Src_vals_int[0:n_points_z-1,0]) \
                                    / (z_space[1:n_points_z]-z_space[0:n_points_z-1])
        
        else: #4:Z*omegac
            dsdc[j,0:n_points_z-1] = ( z_space[1:n_points_z]*Src_vals_int[1:n_points_z,1] \
                                            /Src_vals_int[1:n_points_z,0] \
                                        -z_space[0:n_points_z-1]*Src_vals_int[0:n_points_z-1,1] \
                                            /Src_vals_int[0:n_points_z-1,0] ) \
                                    / (z_space[1:n_points_z]-z_space[0:n_points_z-1])


        # get derivative
        dsdc[j,n_points_z-1] = dsdc[j,n_points_z-2]
        dsdc[j,n_points_z] = dsdc[j,0]

        # #start integration
        y_int[j] = y_int[j]-0.5*np.inner( ( dsdc[j,0:n_points_z-2]*cdf001[0:n_points_z-2] \
                                            + dsdc[j,1:n_points_z-1]*cdf001[1:n_points_z-1] ), \
                                        ( 0.5*(z_space[2:n_points_z]+z_space[1:n_points_z-1]) \
                                            - 0.5*(z_space[1:n_points_z-1]+z_space[0:n_points_z-2]) ) )
        
    y_int[0:nScalars] = y_int[0:nScalars]-dsdc[0:nScalars,n_points_z-1]*cdf001[n_points_z-1] \
                    *( z_space[n_points_z-1] - z_space[n_points_z-2] )/2.0 \
                - dsdc[0:nScalars,n_points_z]*cdf001[0]*(z_space[1]-z_space[0])/2.0 \
                + Src_vals_int[n_points_z-1,0:nScalars]

    for j in range(nYis):
        dYdc[j,0:n_points_z-1] = (Yi_vals_int[1:n_points_z,j]-Yi_vals_int[0:n_points_z-1,j]) \
                                / (z_space[1:n_points_z]-z_space[0:n_points_z-1])
        
        dYdc[j,n_points_z-1] = dYdc[j,n_points_z-2]
        dYdc[j,n_points_z] = dYdc[j,0]

        #start integration
        Yi_int[j] = Yi_int[j]-0.5*np.inner( ( dYdc[j,0:n_points_z-2]*cdf001[0:n_points_z-2] \
                                                +dYdc[j,1:n_points_z-1]*cdf001[1:n_points_z-1] ), \
                                            ( 0.5*(z_space[2:n_points_z]+z_space[1:n_points_z-1]) \
                                                - 0.5*(z_space[1:n_points_z-1]-z_space[0:n_points_z-2]) ) )

    Yi_int[0:nYis] = Yi_int[0:nYis]-dYdc[0:nYis,n_points_z-1]*cdf001[n_points_z-1] \
                    *(z_space[n_points_z-1]-z_space[n_points_z-2])/2.0 \
                - dYdc[0:nYis,n_points_z]*cdf001[0]*(z_space[1]-z_space[0])/2.0 \
                + Yi_vals_int[n_points_z-1,0:nYis]

    np.maximum( abs(Yi_int), 0.0, out=Yi_int )

    return y_int,Yi_int


# def plackettFunc(c_point,z_point,alpha_z,beta_z,alpha_c,beta_c,theta,CDF_C,CDF_Z,plaF):
#     CDF_Z = beta.cdf(z_point,alpha_z,beta_z)
#     CDF_C = beta.cdf(c_point,alpha_c,beta_c)

#     S = 1.0+(theta-1.0)*(CDF_Z+CDF_C)

#     if (theta < 1.0e-5):
#         th = 1.0e-5
#     else:
#         th = theta
    
#     if ( abs(1.0-th) < 1.0e-5 ):
#         plaF = 1.0
#     else:
#         plaF = th*( 1.0+(th-1.0)*(CDF_Z+CDF_C-2.0*(CDF_Z*CDF_C)) ) \
#                 /( ( S*S-(4.0*th*(th-1.0)*CDF_Z*CDF_C) )**(3.0/2.0) )

#     return CDF_C,CDF_Z,plaF

def plackettFunc(theta,CDF_C,CDF_Z,plaF):
    S = 1.0+(theta-1.0)*(CDF_Z+CDF_C)

    if (theta < 1.0e-5):
        th = 1.0e-5
    else:
        th = theta
    
    if ( abs(1.0-th) < 1.0e-5 ):
        plaF = 1.0
    else:
        plaF = th*( 1.0+(th-1.0)*(CDF_Z+CDF_C-2.0*(CDF_Z*CDF_C)) ) \
                /( ( S*S-(4.0*th*(th-1.0)*CDF_Z*CDF_C) )**(3.0/2.0) )

    return plaF


def int_point(z_mean,c_mean,c_var,z_var,co_var,y_int,z_space,c_space,Src_vals,Yi_int,Yi_vals,theta):
    # print("int_point")
    CDF_C = np.zeros(n_points_c+1)
    CDF_Z = np.zeros(n_points_z+1)
    # plaF = np.zeros((n_points_z+1, n_points_c+1))
    #note: 不带copula的时候plaF都等于1.0
    plaF = np.ones((n_points_z+1, n_points_c+1))
    dPsidc = np.zeros((n_points_z+1, n_points_c+1, nScalars))
    Psi = np.zeros((n_points_z+1, n_points_c+1, nScalars))
    Q_int = np.zeros((n_points_z+1, nScalars))
    dQdz = np.zeros((n_points_z+1, nScalars))
    YiPsi = np.zeros((n_points_z+1, n_points_c+1, nYis))
    YiQ_int = np.zeros((n_points_z+1, nYis))
    dYiQdz = np.zeros((n_points_z+1, nYis))
    dYiPsidc = np.zeros((n_points_z+1, n_points_c+1, nYis))

    # rho = co_var/(sqrt(z_var)*sqrt(c_var))

    alpha_z = z_mean*( ( (z_mean*(1.0-z_mean))/z_var ) -1.0 )
    alpha_c = c_mean*( ( (c_mean*(1.0-c_mean))/c_var ) -1.0 )
    beta_z = (1.0-z_mean)*( ( (z_mean*(1.0-z_mean))/z_var ) -1.0 )
    beta_c = (1.0-c_mean)*( ( (c_mean*(1.0-c_mean))/c_var ) -1.0 )

    if (z_mean == 1.0):
        print("z_mean = 1.0")
    if (c_mean == 1.0):
        print("z_mean = 1.0")
    if (z_var==1.0):
        print("z_var = 0.0")
    if (c_var == 0.0):
        print("c_var = 0.0")

    #-----assign values for Psi and obtain CDFs-------

    CDF_Z[1:n_points_z-1] = beta.cdf((z_space[1:n_points_z-1]+z_space[2:n_points_z])/2.0, \
                                    alpha_z,beta_z)
    CDF_Z[n_points_z-1] = beta.cdf((z_space[n_points_z-2]+3.0*z_space[n_points_z-1])/4.0, \
                                    alpha_z,beta_z)
    CDF_Z[n_points_z] = beta.cdf((3.0*z_space[0]+z_space[1])/4.0,alpha_z,beta_z)

    CDF_C[1:n_points_c-1] = beta.cdf((c_space[1:n_points_c-1]+c_space[2:n_points_c])/2.0, \
                                    alpha_c,beta_c)
    CDF_C[n_points_c-1] = beta.cdf((c_space[n_points_c-2]+3.0*c_space[n_points_c-1])/4.0, \
                                    alpha_c,beta_c)
    CDF_C[n_points_c] = beta.cdf((3.0*c_space[0]+c_space[1])/4.0,alpha_c,beta_c)


    for k in range(4,nScalars):
        Psi[0:n_points_z,0:n_points_c,k] = Src_vals[:,:,k]*plaF[0:n_points_z,0:n_points_c]
    
    Psi[0:n_points_z,0:n_points_c,1] = Src_vals[:,:,1]*plaF[0:n_points_z,0:n_points_c]/Src_vals[:,:,0]
    Psi[0:n_points_z,0:n_points_c,2] = Src_vals[:,:,1]*plaF[0:n_points_z,0:n_points_c] \
                                        / Src_vals[:,:,0]*Src_vals[:,:,nScalars-1]
    for j in range(n_points_c-1):
        Psi[0:n_points_z,j,2] = Psi[0:n_points_z,j,2]*c_space[j]
    Psi[0:n_points_z,0:n_points_c,3] = Src_vals[:,:,1]*plaF[0:n_points_z,0:n_points_c]/Src_vals[:,:,0]
    for i in range(n_points_z-1):
        Psi[i,0:n_points_c,3] = Psi[i,0:n_points_c,3]*z_space[i]

    for k in range(nYis):
        YiPsi[0:n_points_z,0:n_points_c,k] = Yi_vals[:,:,k]*plaF[0:n_points_z,0:n_points_c]


    # #-----2:omegac 3:c*omegac 4:z*omegac-------#
    # #-----5:Cp 6:MW 7:formEnthal 8:T------#
    for j in range(n_points_c-1):
        dPsidc[0:n_points_z,j,:] = (Psi[0:n_points_z,j+1,:]-Psi[0:n_points_z,j,:]) \
                                    /(c_space[j+1]-c_space[j])
    j = n_points_c-1
    dPsidc[:,j,:] = dPsidc[:,j-1,:]  #CZ:half rectangle
    #first point in c space, dPsidc(1:n_points_z,n_points_c+1)
    j = n_points_c
    dPsidc[:,j,:] = dPsidc[:,0,:]  #CZ:half rectangle

    # for k in range(nScalars):
    for j in range(n_points_c-2):
        Q_int[0:n_points_z,0:nScalars] = Q_int[0:n_points_z,0:nScalars]-0.5*(dPsidc[0:n_points_z,j,0:nScalars] \
                                            *CDF_C[j]+dPsidc[0:n_points_z,j+1,0:nScalars]*CDF_C[j+1]) \
                                        * ((c_space[j+2]+c_space[j+1])*0.5 - (c_space[j+1]+c_space[j])*0.5)

    Q_int[0:n_points_z,0:nScalars] = Q_int[0:n_points_z,0:nScalars]-dPsidc[0:n_points_z,n_points_c-1,0:nScalars] \
                                        *CDF_C[n_points_c-1]*(c_space[n_points_c-1]-c_space[n_points_c-1])/2.0 \
                                    - dPsidc[0:n_points_z,n_points_c,0:nScalars]*CDF_C[n_points_c]*(c_space[1] \
                                        -c_space[0])/2.0 + Psi[0:n_points_z,n_points_c-1,0:nScalars]

    for i in range(n_points_z-1):
        dQdz[i,0:nScalars] = ((Q_int[i+1,0:nScalars]-Q_int[i,0:nScalars])) \
                                /(z_space[i+1]-z_space[i])

    # last point in z space, dQdz(n_points_z)
    i = n_points_z-1
    dQdz[i,0:nScalars] = dQdz[0,0:nScalars]
    #     #     y_int[k] = y_int[k]-0.5*(dQdz[i,k]*CDF_Z[i]+dQdz[i+1,k]*CDF_Z[i+1]) \
    #     #                 *( (z_space[i+2]+z_space[i+1])*0.5-(z_space[i+1]+z_space[i])*0.5 )
    for k in range(nScalars):
        y_int[k] = y_int[k]-0.5*np.inner( (dQdz[0:n_points_z-2,k]*CDF_Z[0:n_points_z-2] \
                                            +dQdz[1:n_points_z-1,k]*CDF_Z[1:n_points_z-1]), \
                                        ( (z_space[2:n_points_z]+z_space[1:n_points_z-1])*0.5 \
                                            -(z_space[1:n_points_z-1]+z_space[0:n_points_z-2])*0.5 ) )


    y_int[0:nScalars] = y_int[0:nScalars]+dQdz[n_points_z-1,0:nScalars]*CDF_Z[n_points_z-1] \
                    *(z_space[n_points_z-1]-z_space[n_points_z-2])/2.0 \
                + dQdz[n_points_z,0:nScalars]*CDF_Z[n_points_z]*(z_space[1]-z_space[0])/2.0 \
                + Q_int[n_points_z-1,0:nScalars]

    

    #-----Yis------#
    for k in range(nYis):
        dYiPsidc[0:n_points_z,0:n_points_c-1,k] = (YiPsi[0:n_points_z,1:n_points_c,k] \
                                                        -YiPsi[0:n_points_z,0:n_points_c-1,k]) \
                                                    /(c_space[1:n_points_c]-c_space[0:n_points_c-1])
    j = n_points_c-1
    dYiPsidc[0:n_points_z,j,0:nYis] = dYiPsidc[0:n_points_z,j-1,0:nYis]
    j = n_points_c
    dYiPsidc[0:n_points_z,j,0:nYis] = dYiPsidc[0:n_points_z,0,0:nYis]

    for j in range(n_points_c-2):
        YiQ_int[0:n_points_z,0:nYis] = YiQ_int[0:n_points_z,0:nYis] \
                                        -0.5*(dYiPsidc[0:n_points_z,j,0:nYis]*CDF_C[j] \
                                            +dYiPsidc[0:n_points_z,j+1,0:nYis]*CDF_C[j+1]) \
                                        * ((c_space[j+2]+c_space[j+1])*0.5 \
                                            -(c_space[j+1]+c_space[j])*0.5)

    YiQ_int[0:n_points_z,0:nYis] = YiQ_int[0:n_points_z,0:nYis]-dYiPsidc[0:n_points_z,n_points_c-1,0:nYis] \
                                        *CDF_C[n_points_c]*(c_space[n_points_c-1]-c_space[n_points_c-2])/2.0 \
                                    - dYiPsidc[0:n_points_z,n_points_c,0:nYis]*CDF_C[n_points_c] \
                                        *(c_space[1]-c_space[0])/2.0 \
                                    + YiPsi[0:n_points_z,n_points_c-1,0:nYis]

    for i in range(n_points_z-1):
        dYiQdz[i,0:nYis] = ( YiQ_int[i+1,0:nYis]-YiQ_int[i,0:nYis] )/(z_space[i+1]-z_space[i])
    # last point in z space, dYiQdz(n_points_z)
    i = n_points_z-1
    dYiQdz[i,0:nYis] = dYiQdz[i-1,0:nYis]
    # first point in z space, dYiQdz(n_points_z+1)
    i = n_points_z
    dYiQdz[i,0:nYis] = dYiQdz[0,0:nYis]

    # start integration
    for k in range(nYis):
        Yi_int[k] = Yi_int[k]-0.5*np.inner( (dYiQdz[0:n_points_z-2,k]*CDF_Z[0:n_points_z-2] \
                                                +dYiQdz[1:n_points_z-1,k]*CDF_Z[1:n_points_z-1]), \
                                            ((z_space[2:n_points_z]+z_space[1:n_points_z-1])*0.5 \
                                                -(z_space[1:n_points_z-1]+z_space[0:n_points_z-2])*0.5) )

    Yi_int[0:nYis] = Yi_int[0:nYis]+dYiQdz[n_points_z-1,0:nYis]*CDF_Z[n_points_z-1] \
                        *(z_space[n_points_z-1]-z_space[n_points_z-2])/2.0 \
                    + dYiQdz[n_points_z,0:nYis]*CDF_Z[n_points_z] \
                        *(z_space[1]-z_space[0])/2.0 \
                    +YiQ_int[n_points_z-1,0:nYis]


    np.maximum( abs(Yi_int), 0.0, out=Yi_int )
    
    return y_int,Yi_int
