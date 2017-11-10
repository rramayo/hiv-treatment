# -*- coding:utf-8 -*-

#*************************************************************
#****       0: Saudáveis susceptíveis (S)                 ****
#****       1: Infectadas produtivas (A)                  ****
#****       2: Saudáveis RT-Inibidas (Srti)               ****
#****       3: Saudáveis P-Inibidas (Spi)                 ****
#****       4: Saudáveis resistentes (SR)                 ****
#****       5: Infectadas não produtivas (B)              ****
#****       6: Mortas (M)                                 ****
#*************************************************************

# Module principal
  #integer     iseed,seed,idelay(0:1100,0:1100),idelayv2(0:1100,0:1100)
  #integer     celulas(0:1100,0:1100),celulas1(0:1100,0:1100),virus(0:1100,0:1100),virus1(0:1100,0:1100)
  #real        iviz(0:8),Ds(0:10000),Da(0:10000),Dvi(0:10000),Dvni(0:10000),Db(0:10000),Dan(0:1100),Dvm(0:10000),Dm(0:10000),Dinf(0:10000),Dv(0:10000)
  #real        Dsq(0:10000),Daq(0:10000),Dbq(0:10000),Dviq(0:10000),Dmq(0:10000),Dinfq(0:10000),Danq(0:10000),Dvq(0:10000)
  #real        DDs(0:10000),DDa(0:10000),DDan(0:10000),DDvni(0:10000),DDb(0:10000),DDm(0:10000),Dvmq(0:10000),DDvm(0:10000),DDv(0:10000)
  #real        DDinf(0:10000),Dvniq(0:10000),DDvi(0:10000)
  #integer     idelayv(0:1100,0:1100),i1,j1,i,j,k1,j2,l1
  #real        d,m,f,a,y,pHIV,pnasce,pinf,preg,preg_inf,preg_virus,preinf,ppi,prti,pb,p0rti,y4,y5,y6,y7,y8,y9,y10,y11,y1,y2,y3,y12,y13,y14,y15,y16,y17
  #integer     iDs,iDa,iDb,iDvi,iDvm,iDinf,iDm,iDvni,iDvq
  #integer     b,L,iadelay,ivdelay,IRB,IRA,t0,t1,t2,t3,t4,t5,tmax,Nsamp,it,it1,it2,k,it3
# end Module principal

import os
from random import random
from multiprocessing import Process


def therapy_hiv(L=None, iadelay=None, pHIV=None, pnasce=None,
                pinf=None, IRA=None, IRB=None, tmax=None,
                Nsamp=None, t0=None, p0rti=None, p0pi=None, b=None):
    # para ser calculado uma única vez.

    flag_d = 0

    Ds = {}
    Da = {}
    Db = {}
    Dm = {}
    Dan = {}
    Dinf = {}
    Dsq = {}
    Daq = {}
    Dbq = {}
    Dmq = {}
    Danq = {}
    Dinfq = {}
    DDs = {}
    DDa = {}
    DDb = {}
    DDm = {}
    DDinf = {}

    celulas = {}
    celulas1 = {}
    idelay = {}
    iviz = {}

    if not os.path.exists('./data'):
        os.mkdir('./data')

    if (p0pi > 0) and (p0rti > 0):
        path = './data/combined'
        dir_matrix = path + '/matrix'
        dir_mean = path + '/mean'

        if not os.path.exists(path):
            os.mkdir(path)
        if not os.path.exists(dir_mean):
            os.mkdir(dir_mean)
        if not os.path.exists(dir_matrix):
            os.mkdir(dir_matrix)

    elif (p0pi > 0) and (p0rti == 0):
        path = './data/monotherapy_PI'
        dir_matrix = path + '/matrix'
        dir_mean = path + '/mean'

        if not os.path.exists(path):
            os.mkdir(path)
        if not os.path.exists(dir_mean):
            os.mkdir(dir_mean)
        if not os.path.exists(dir_matrix):
            os.mkdir(dir_matrix)

    elif (p0pi == 0) and (p0rti > 0):
        path = './data/monotherapy_RTI'
        dir_matrix = path + '/matrix'
        dir_mean = path + '/mean'

        if not os.path.exists(path):
            os.mkdir(path)
        if not os.path.exists(dir_mean):
            os.mkdir(dir_mean)
        if not os.path.exists(dir_matrix):
            os.mkdir(dir_matrix)

    for it in range(0, tmax + 1):
        Ds[it] = 0.0
        Da[it] = 0.0
        Db[it] = 0.0
        Dm[it] = 0.0
        Dan[it] = 0.0
        Dinf[it] = 0.0
        Dsq[it] = 0.0
        Daq[it] = 0.0
        Dbq[it] = 0.0
        Dmq[it] = 0.0
        Danq[it] = 0.0
        Dinfq[it] = 0.0
        DDs[it] = 0.0
        DDa[it] = 0.0
        DDb[it] = 0.0
        DDm[it] = 0.0
        DDinf[it] = 0.0

    for k1 in range(0, Nsamp):
        d = flag_d

        iDs = 0
        iDa = 0
        for i in range(-1, L + 1):
            for j in range(-1, L + 1):
                celulas[(i, j)] = 0
                celulas1[(i, j)] = 0

                if (-1 < i < L) and (-1 < j < L):
                    y = random()
                    if (y > pHIV):
                        iDs = iDs + 1
                    else:
                        celulas[(i, j)] = 1
                        iDa = iDa + 1
                    idelay[(i, j)] = 0

        Ds[0] = Ds[0] + (1.0*iDs)/(1.0*L*L)
        Da[0] = Da[0] + (1.0*iDa)/(1.0*L*L)
        Dinf[0] = Da[0]
        Dsq[0] = Dsq[0] + (1.0*iDs/(1.0*L*L))**2
        Daq[0] = Daq[0] + (1.0*iDa/(1.0*L*L))**2
        Dinfq[0] = Dinfq[0] + (1.0*iDa/(1.0*L*L))**2

        s_matrix = open(dir_matrix + '/s-Nsamp-{}-tmax-{}-t0-{}-pHIV-{}-pnasce-{}-pinf-{}-p0rti-{}-p0pi-{}-L-{}.dat'.format(Nsamp, tmax, t0, str(pHIV).replace('.','_'), str(pnasce).replace('.','_'), str(pinf).replace('.','_'), str(p0rti).replace('.','_'), str(p0pi).replace('.','_'), L), 'a')
        i_matrix = open(dir_matrix + '/i-Nsamp-{}-tmax-{}-t0-{}-pHIV-{}-pnasce-{}-pinf-{}-p0rti-{}-p0pi-{}-L-{}.dat'.format(Nsamp, tmax, t0, str(pHIV).replace('.','_'), str(pnasce).replace('.','_'), str(pinf).replace('.','_'), str(p0rti).replace('.','_'), str(p0pi).replace('.','_'), L), 'a')
        m_matrix = open(dir_matrix + '/m-Nsamp-{}-tmax-{}-t0-{}-pHIV-{}-pnasce-{}-pinf-{}-p0rti-{}-p0pi-{}-L-{}.dat'.format(Nsamp, tmax, t0, str(pHIV).replace('.','_'), str(pnasce).replace('.','_'), str(pinf).replace('.','_'), str(p0rti).replace('.','_'), str(p0pi).replace('.','_'), L), 'a')

        for it in range(1, tmax+1):
            # print(it)

            for j2 in range(0, L + 1):
                celulas[(0, j2)] = celulas[(L, j2)]
                celulas[(L, j2)] = celulas[(1, j2)]
                celulas[(j2, 0)] = celulas[(j2, L)]
                celulas[(j2, L)] = celulas[(j2, 1)]
                celulas[(0, 0)] = celulas[(0, L)]
                celulas[(0, L)] = celulas[(0, 1)]
                celulas[(L, 0)] = celulas[(L, L)]
                celulas[(L, L)] = celulas[(1, L)]

            iDs = 0
            iDa = 0
            iDb = 0
            iDm = 0

            for k in range(0, L):
                for j in range(0, L):
                    contaa = 0
                    contab = 0
                    if (celulas[(k, j)] == 0): #***********células S******************
                        iDs = iDs + 1
                        iviz[1] = celulas[(k - 1, j - 1)]
                        iviz[2] = celulas[(k - 1, j)]
                        iviz[3] = celulas[(k - 1, j + 1)]
                        iviz[4] = celulas[(k, j - 1)]
                        iviz[5] = celulas[(k, j + 1)]
                        iviz[6] = celulas[(k + 1, j - 1)]
                        iviz[7] = celulas[(k + 1, j)]
                        iviz[8] = celulas[(k + 1, j + 1)]
                        for m in range(1, 9):
                            if (iviz[m] == 1):
                                contaa = contaa + 1
                            if (iviz[m] == 5):
                                contab = contab + 1
                        if(it < t0) :
                            if ((contaa >= IRA) or (contab >= IRB)):
                                celulas1[(k, j)] = 1
                            else:
                                celulas1[(k, j)] = 0
                        else:
                            if d == 0:
                                d = len([cel for pos, cel in celulas.items() \
                                         if (cel in [1, 5]) and ((-1 < pos[0] < L) and (-1 < pos[1] < L))])/(L*L)
                                print('d: {}'.format(d)) #para calcular d

                            if ((contaa >= IRA) or (contab >= IRB)):
                                y1=random()
                                y2=random()
                                if ((y1 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y2 >= p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                    celulas1[(k,j)]=1
                                else:
                                    if ((y1 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y2 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                        celulas1[(k,j)]=3
                                    else:
                                        if ((y1 < p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y2 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                            celulas1[(k,j)]=4
                                        else:
                                            celulas1[(k,j)]=2
                            else:
                                y3=random()
                                y4=random()
                                if ((y3 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y4 >= p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                    celulas1[(k,j)]=0
                                else:
                                    if ((y3 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y4 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                        celulas1[(k,j)]=3
                                    else:
                                        if ((y3 < p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y4 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                            celulas1[(k,j)]=4
                                        else:
                                            celulas1[(k,j)]=2
                    if  (celulas[(k,j)] == 1): #***********células A******************
                        iDa = iDa + 1
                        idelay[(k,j)] = idelay[(k,j)] + 1
                        if (idelay[(k,j)] == iadelay):
                            celulas1[(k,j)] = 5
                            idelay[(k,j)] = 0
                        else:
                            celulas1[(k,j)]=1
                    if (celulas[(k,j)] == 5): #***********células B******************
                        iDb=iDb+1
                        celulas1[(k,j)]=6
                    if  (celulas[(k,j)] == 6): #***********células M******************
                        iDm=iDm+1
                        if (it < t0) :
                            y12=random()
                            y13=random()
                            if (y12 < pnasce) : #3
                                if (y13 < pinf) : #4
                                    celulas1[(k,j)]=1
                                else: #4
                                    celulas1[(k,j)]=0
                            else: #3
                                celulas1[(k,j)]=6
                        else:
                            if d == 0:
                                d = len([cel for pos, cel in celulas.items() \
                                         if (cel in [1, 5]) and ((-1 < pos[0] < L) and (-1 < pos[1] < L))])/(L*L)
                                print('d: {}'.format(d)) #para calcular d

                            y14=random()
                            y15=random()
                            y16=random()
                            if (y14 < pnasce) : #5
                                if (y15 < pinf) : #6
                                    if (y16 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) : #7
                                        celulas1[(k,j)]=1
                                    else:
                                        celulas1[(k,j)]=2
                                else: #6
                                    celulas1[(k,j)]=0
                            else: #5
                                celulas1[(k,j)]=6
                    if (celulas[(k,j)] == 2): #**********células RTI-Inibidas*********
                        iDs=iDs+1
                        iviz[1]=celulas[(k-1,j-1)]
                        iviz[2]=celulas[(k-1,j)]
                        iviz[3]=celulas[(k-1,j+1)]
                        iviz[4]=celulas[(k,j-1)]
                        iviz[5]=celulas[(k,j+1)]
                        iviz[6]=celulas[(k+1,j-1)]
                        iviz[7]=celulas[(k+1,j)]
                        iviz[8]=celulas[(k+1,j+1)]
                        for m in range(1, 9):
                            if (iviz[m] == 1):
                                contaa = contaa + 1
                            if (iviz[m] == 5):
                                contab = contab + 1
                        if (it < t0):
                            if ((contaa >= IRA) or (contab >= IRB)):
                                celulas1[(k,j)]=1
                            else:
                                celulas1[(k,j)]=0
                        else:
                            if d == 0:
                                d = len([cel for pos, cel in celulas.items() \
                                         if (cel in [1, 5]) and ((-1 < pos[0] < L) and (-1 < pos[1] < L))])/(L*L)
                                print('d: {}'.format(d)) #para calcular d

                            if ((contaa >= IRA) or (contab >= IRB)):
                                y5=random()
                                y6=random()
                                if ((y5 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y6 >= p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                    celulas1[(k,j)]=1
                                else:
                                    if ((y5 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y6 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                        celulas1[(k,j)]=4
                                    else:
                                        if ((y5 < p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y6 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                            celulas1[(k,j)]=4
                                        else:
                                            celulas1[(k,j)]=2
                            else:
                                y7=random()
                                y8=random()
                                if ((y7 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y8 >= p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                    celulas1[(k,j)]=0
                                else:
                                    if ((y7 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y8 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                        celulas1[(k,j)]=4
                                    else:
                                        if ((y7 < p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y8 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                            celulas1[(k,j)]=4
                                        else:
                                            celulas1[(k,j)]=2
                    if (celulas[(k,j)] == 3) : #***********células PI-Inibidas*********
                        iDs=iDs+1
                        iviz[1]=celulas[(k-1,j-1)]
                        iviz[2]=celulas[(k-1,j)]
                        iviz[3]=celulas[(k-1,j+1)]
                        iviz[4]=celulas[(k,j-1)]
                        iviz[5]=celulas[(k,j+1)]
                        iviz[6]=celulas[(k+1,j-1)]
                        iviz[7]=celulas[(k+1,j)]
                        iviz[8]=celulas[(k+1,j+1)]
                        for m in range(1, 9):
                            if (iviz[m] == 1):
                                contaa=contaa+1
                            if (iviz[m] == 5):
                                contab=contab+1
                        if (it < t0) :
                            if ((contaa >= IRA) or (contab >= IRB)) :
                                celulas1[(k,j)]=1
                            else:
                                celulas1[(k,j)]=0
                        else:
                            if d == 0:
                                d = len([cel for pos, cel in celulas.items() \
                                         if (cel in [1, 5]) and ((-1 < pos[0] < L) and (-1 < pos[1] < L))])/(L*L)
                                print('d: {}'.format(d)) #para calcular d

                            if ((contaa >= IRA) or (contab >= IRB)) :
                                y9=random()
                                y10=random()
                                if ((y9 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y10 >= p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                    celulas1[(k,j)]=1
                                else:
                                    if ((y9 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y10 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                        celulas1[(k,j)]=3
                                    else:
                                        if ((y9 < p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y10 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                            celulas1[(k,j)]=4
                                        else:
                                            celulas1[(k,j)]=4
                            else:
                                y11=random()
                                y17=random()
                                if ((y11 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y17 >= p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                    celulas1[(k,j)]=0
                                else:
                                    if ((y11 >= p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y17 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                        celulas1[(k,j)]=3
                                    else:
                                        if ((y11 < p0rti*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2) and (y17 < p0pi*(((1.0*(iDa+iDb)/(1.0*L*L))/d)**(2*d/(1-d)))*((1-(1.0*(iDa+iDb)/(1.0*L*L)))/(1-d))**2)) :
                                            celulas1[(k,j)]=4
                                        else:
                                            celulas1[(k,j)]=4
                    if (celulas[(k,j)] == 4): #*************células SR***************
                        iDs = iDs + 1
                        celulas1[(k,j)] = 0

            celulas = celulas1.copy()

            ds_ = (1.0*iDs)/(1.0*L*L)
            da_ = (1.0*iDa)/(1.0*L*L)
            db_ = (1.0*iDb)/(1.0*L*L)
            dm_ = (1.0*iDm)/(1.0*L*L)
            dinf_ = (da_ + db_)
            dsq_ = (1.0*iDs/(1.0*L*L))**2
            daq_ = (1.0*iDa/(1.0*L*L))**2
            dbq_ = (1.0*iDb/(1.0*L*L))**2
            dmq_ = (1.0*iDm/(1.0*L*L))**2
            dinfq_ = (1.0*iDa/(1.0*L*L)+1.0*iDb/(1.0*L*L))**2

            s_matrix.write(str(ds_) + "\t")
            i_matrix.write(str(dinf_) + "\t")
            m_matrix.write(str(dm_) + "\t")

            Ds[it] = Ds[it] + ds_
            Da[it] = Da[it] + da_
            Db[it] = Db[it] + db_
            Dm[it] = Dm[it] + dm_
            Dinf[it] = Dinf[it] + dinf_
            Dsq[it] = Dsq[it] + dsq_
            Daq[it] = Daq[it] + daq_
            Dbq[it] = Dbq[it] + dbq_
            Dmq[it] = Dmq[it] + dmq_
            Dinfq[it] = Dinfq[it] + dinfq_

        s_matrix.write('\n')
        i_matrix.write('\n')
        m_matrix.write('\n')
        s_matrix.close()
        i_matrix.close()
        m_matrix.close()

    s_mean = open(dir_mean + '/s-Nsamp-{}-tmax-{}-t0-{}-pHIV-{}-pnasce-{}-pinf-{}-p0rti-{}-p0pi-{}-L-{}.dat'.format(Nsamp, tmax, t0, str(pHIV).replace('.','_'), str(pnasce).replace('.','_'), str(pinf).replace('.','_'), str(p0rti).replace('.','_'), str(p0pi).replace('.','_'), L), 'w')
    i_mean = open(dir_mean + '/i-Nsamp-{}-tmax-{}-t0-{}-pHIV-{}-pnasce-{}-pinf-{}-p0rti-{}-p0pi-{}-L-{}.dat'.format(Nsamp, tmax, t0, str(pHIV).replace('.','_'), str(pnasce).replace('.','_'), str(pinf).replace('.','_'), str(p0rti).replace('.','_'), str(p0pi).replace('.','_'), L), 'w')
    m_mean = open(dir_mean + '/m-Nsamp-{}-tmax-{}-t0-{}-pHIV-{}-pnasce-{}-pinf-{}-p0rti-{}-p0pi-{}-L-{}.dat'.format(Nsamp, tmax, t0, str(pHIV).replace('.','_'), str(pnasce).replace('.','_'), str(pinf).replace('.','_'), str(p0rti).replace('.','_'), str(p0pi).replace('.','_'), L), 'w')

    for it in range(0, tmax + 1):
        Ds[it]=Ds[it]/(1.0*Nsamp)
        Da[it]=Da[it]/(1.0*Nsamp)
        Db[it]=Db[it]/(1.0*Nsamp)
        Dm[it]=Dm[it]/(1.0*Nsamp)
        Dinf[it]=Dinf[it]/(1.0*Nsamp)
        Dsq[it]=Dsq[it]/(1.0*Nsamp)
        Daq[it]=Daq[it]/(1.0*Nsamp)
        Dbq[it]=Dbq[it]/(1.0*Nsamp)
        Dmq[it]=Dmq[it]/(1.0*Nsamp)
        Dinfq[it]=Dinfq[it]/(1.0*Nsamp)
        DDs[it]=(abs(Dsq[it]-Ds[it]**2))**0.5
        DDa[it]=(abs(Daq[it]-Da[it]**2))**0.5
        DDb[it]=(abs(Dbq[it]-Db[it]**2))**0.5
        DDm[it]=(abs(Dmq[it]-Dm[it]**2))**0.5
        DDinf[it]=(abs(Dinfq[it]-Dinf[it]**2))**0.5

        s_mean.write(str(it)+"\t"+str(Ds[it])+"\t"+str(DDs[it])+"\n")
        i_mean.write(str(it)+"\t"+str(Dinf[it])+"\t"+str(DDinf[it])+"\n")
        m_mean.write(str(it)+"\t"+str(Dm[it])+"\t"+str(DDm[it])+"\n")

    s_mean.close()
    i_mean.close()
    m_mean.close()

    print('--- End ---')


def extract_tasks(jobs, threads):
    tasks = []
    while jobs:
        tmp = []
        for i in range(threads):
            if jobs:
                tmp.append(jobs.pop())
            else:
                break
        if tmp:
            tasks.append(tmp)
    return tasks


if __name__ == '__main__':
    keywords_args = {
        'L': 700,
        'iadelay': 4,
        'pHIV': 0.05,
        'pnasce': 0.99,
        'pinf': 0.00001,
        'IRA': 1,
        'IRB': 4,
        'tmax': 1000,
        'Nsamp': 1,
        't0': 300,
        'p0rti': 0.7,
        'p0pi': 0.7,
        'b': 2
    }

    N = 100
    M = 1000
    B = 1000000

    for pHIV in [50]: #range(10, 160, 5):
        pHIV = pHIV/M

        for pnasce in [99]: #range(90, 100, 1):
            pnasce = pnasce/N

            for pinf in [10]: #range(100, 0, -10):
                pinf = pinf/B

                for t0 in [100]: #range(100, 600, 100):

                  for p0pi in [70]: #range(0, N + 1):
                        p0pi = p0pi/N

                        jobs = []

                        for p0rti in [70]: #range(0, N + 1):
                            p0rti = p0rti/N

                            temp = {
                            'pHIV': pHIV,
                            'pnasce': pnasce,
                            'pinf': pinf,
                            't0': t0,
                            'p0pi': p0pi,
                            'p0rti': p0rti
                            }

                            kwargs_ = keywords_args.copy()
                            kwargs_.update(temp)

                            jobs.append(Process(target=therapy_hiv, kwargs=kwargs_))

                        threads = os.cpu_count()
                        tasks = extract_tasks(jobs, threads)

                        for task in tasks:
                            for t in task:
                                t.start()
                            for t in task:
                                t.join()
