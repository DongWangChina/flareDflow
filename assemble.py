import tarfile
import numpy as np

def assemble(cbDict,solFln):
    tar = tarfile.open(solFln,'r')
    tar.extractall()
    tar.close()

    for i in range(cbDict["n_points_h"]):
        #n_points_h=1
        print('Reading unit' + '%02d' % 1 + '_h' + '%02d' % (i+1) + ' ... \n')
        M = np.loadtxt("./" + 'unit01_h' + '%02d' % (i+1) + '.dat')

        for j in range(1,cbDict["int_pts_z"]):
            print('Reading unit' + '%02d' % (j+1) + '_h' + '%02d' % (i+1) + ' ... \n')
            tmp = np.loadtxt('./' + 'unit' + '%02d' % (j+1) + '_h' + '%02d' % (i+1) + '.dat')
            M = np.insert(tmp,0,M,axis=0) 

    # remove unwanted columns - h(12),qdot(13),Yc_max(14)
    n_column = np.shape(tmp)[1]
    if(cbDict['scaled_PV']): 
        rm_list = [n_column-6,n_column-5,n_column-4] 
    else:
        rm_list = [n_column-6,n_column-5]
    MM = np.delete(M,rm_list,axis=1)

    # write assembled table
    fln = './' +  cbDict['output_fln'] 
    print('Writing assembled table ...')
    with open(fln,'a') as strfile:
        #写入混合物分数z
        strfile.write(str(cbDict['int_pts_z']) + '\n')
        np.savetxt(strfile,cbDict['z'],fmt='%.5E',delimiter='\t') 

        #写入进程变量c
        strfile.write(str(cbDict['int_pts_c']) + '\n')
        np.savetxt(strfile,cbDict['c'],fmt='%.5E',delimiter='\t')
        
        #写入混合物分数z的方差gz
        strfile.write(str(cbDict['int_pts_gz']) + '\n')
        np.savetxt(strfile,cbDict['gz'],fmt='%.5E',delimiter='\t')

        #写入进程变量c的方差gc
        strfile.write(str(cbDict['int_pts_gc']) + '\n')
        np.savetxt(strfile,cbDict['gc'],fmt='%.5E',delimiter='\t')

        #写入混合物分数z和进程变量c的协方差gcz
        strfile.write(str(cbDict['int_pts_gcz']) + '\n')
        np.savetxt(strfile,cbDict['gcz'],fmt='%.5E',delimiter='\t')  

        #写入MM的第5列及以后的数据
        strfile.write(str(MM.shape[1]-5-cbDict['nYis']) + '\t' +
                        str(cbDict['nYis']) + '\n')
        np.savetxt(strfile,MM[:,5:],fmt='%.5E',delimiter='\t') 

        #写入截断后的混合物分数的概率密度分布
        if(cbDict['scaled_PV']):
            d2Yeq_table = np.loadtxt('./d2Yeq_table.dat')
            np.savetxt(strfile,d2Yeq_table,fmt='%.5E')
    strfile.close()

    print("\n Done writing flare.tbl")

    return fln