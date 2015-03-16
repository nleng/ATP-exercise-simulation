#!/usr/bin/env python

import numpy as np
import scipy.integrate as itg

do_exercise = True

max_atp_usage = 300.
# determines the mitochondrial ATP synthesis, healthy controls have 1.0, CFS patients between 0.6-1.0
CFS_factor = 1.

spezies_num = 21
F = 0.096484
R = 0.0083143
T = 298.
Z = np.log(10.)*R*T/F
u = 0.861
nA = 2.5

N_t = 2.97
Q_t = 1.35
c_t = 0.27
deltaT_flxues = 0.001
t_global = -deltaT_flxues  # to write the first value at t=0

ATPusageFactor = 1.
ATPusageFactorStart = 1.
ATPusageFactorMid = 1.
ATPusageFactorEnd = 1.
countala = 0
global_summa = 0.
t_fluxes = []
ATPusageFactor_array = [[], []]
v_O2_array = []
fluxes = [[] for i in range(spezies_num)]
x0 = [0. for i in range(spezies_num)]
x = [x0[i] for i in range(spezies_num)]

ATP_per_NADH = 2.5
NADH_per_pyruvate = 5.6 
ATP_per_2pyruvate = 28.
ATP_per_glucose = 31.
ATP_factor = ATP_per_2pyruvate/ATP_per_glucose

exercise_times = [0.5, 3., 90.]

x0[0] = 5.9979809687956491    # ATP,   
x0[1] = 0.015763420370545012    # freeADP
x0[2] = 4.6399558783491064e-05    # freeAMP
x0[3] = 7.0706813976285634 # ATP_in
x0[4] = 4.490202659365945    # Pi
x0[5] = 16.306561467204542  # Pi_in
x0[6] = 0.06787796999618452    # pyruvate
x0[7] = 3.6369374057978301e-05  # H_in
x0[8] = 0.0001007050433454901    # H, ph-value about 7.0
x0[9] = 1.0240515334433469    # lactate
x0[10] = 2.8293186023714365  # ADP_in
x0[11] = 0.027124566250814747    # IMP
x0[12] = 9.0045425235287695   # creatine
x0[13] = 20.495457476471135    # phosphocreatine
x0[14] = 0.    # purine_loss
x0[15] = 0.0    # NOT USED right now
x0[16] = 0. # blood H efflux
x0[17] = 0. # lactate efflux
x0[18] = 0.46893794871872252 # # NAD
x0[19] = 0.80732835310224171 # U_total = 1.35
x0[20] = 0.040342382907871691 # c_total = 0.27

x = [x0[i] for i in range(spezies_num)]

kLK1 = 7.919796484126334e-05
kLK2 = 0.038
kDH = 37.689000420846966
kC1 = 0.11240931952860435
kC3 = 0.5212238962876438
kC4 = 1951.9793836559998
KmO = 0.12
kSN = 3.0342937497807068
kEX = 3.5079409969651185
KmADP = 0.0035
kPi = 9452.467034848663
kfCK = 1925890.0
kbCK = 1.16017642402
kUT = 0.6865
KmA = 0.15
kAK = 862100.0
kbAK = 340097.172332
kGlyco = 2.756269152398911
kf_lac = 950000.0
kb_lac = 6.33994277742
k_lac_eff = 0.06
AMPd_K = 1.0
kAD = 683.0
kAS = 7.48701711831
kPL = 2.18
kIS = 0.0485
Km_GTP = 0.1
Km_IMP = 0.3
Km_AMP = 3.057
KmATP_IS = 0.014
KiADP_IS = 0.01
KiAMP_IS = 0.092
KiIMP_IS = 0.18

def getPH(protons):
    return -np.log10(protons/1000.)

def getProtons(pH):
    return pow(10., -pH)*1000.

def ATP_m(ATP):
    return ATP - ATP/(1. + 4./0.024)

def ATP_f(ATP):
    return ATP/(1. + 4./0.024)

def ATP_fi(ATP):
    return ATP/(1. + 4./0.017)

def ADP_m(ADP):
    return ADP - ADP/(1. + 4./0.347)

def ADP_f(ADP):
    return ADP/(1. + 0.38/0.347)

def ADP_fi(ADP):
    return ADP/(1. + 0.38/0.282)

def deltaP(H, H_in):
    pH_diff = getPH(H_in)-getPH(H)
    deltapH = Z*pH_diff
    deltaP = 1./(1.-u) * deltapH
    return deltaP

def v_EFF(H):
    return 13.1*(7. - getPH(H))

def v_LK(H, H_in):
    return kLK1 * (np.exp(kLK2 * deltaP(H, H_in)) - 1.)

def v_DH(NADH, pyruvate):
    NAD = N_t - NADH
    return kDH*pyruvate/(1.+100.*NADH/NAD)**0.8

def v_C1(NADH, UQH2, H, H_in):
    NAD = N_t - NADH
    UQ = Q_t - UQH2
    EmN = -320. + Z/2. * np.log10(NAD/NADH)
    EmU = 85. + Z/2. * np.log10(UQ/UQH2)
    deltaE_C1 = EmU - EmN - deltaP(H, H_in) * 4./2.
    return kC1*deltaE_C1

def v_C3(c2, UQH2, H, H_in):
    UQ = Q_t - UQH2
    c3 = c_t - c2
    Emc = 250. + Z * np.log10(c3/c2)
    EmU = 85. + Z/2. * np.log10(UQ/UQH2)
    deltaE_C3 = Emc - EmU - deltaP(H, H_in) * (4.-2.*u)/2.
    return kC3 * deltaE_C3

def v_C4(c2, H, H_in):
    c3 = c_t - c2
    Emc = 250. + Z * np.log10(c3/c2)
    Ema = Emc + deltaP(H, H_in) * (2.+2.*u)/2.
    A32 = 10.**((Ema-540.)/Z)
    a2 = 0.135 /(1.+A32)
    return kC4 * a2 * c2 /(1.+KmO/0.24)

def v_SN(ATP, ADP, Pi, H, H_in):
    deltaGp = 31.9/F + Z * np.log10(1000. * ATP/(ADP*Pi))
    deltaG = nA * deltaP(H, H_in) - deltaGp
    gamma = 10.**(deltaG/Z)
    return kSN * (gamma - 1.)/(gamma + 1.)

def v_EX(ATP, ADP, ATP_in, ADP_in, H, H_in):
    deltapH = Z*(getPH(H_in)-getPH(H))
    deltaPhi = -(deltaP(H, H_in) - deltapH)
    phi = -0.35*deltaPhi
    phi_in = 0.65*deltaPhi
    ATP = ATP_f(ATP)
    ADP = ADP_f(ADP)
    ATP_in = ATP_fi(ATP_in)
    ADP_in = ADP_fi(ADP_in)
    return kEX*(ADP/(ADP + ATP*10.**(-phi/Z)) - ADP_in/(ADP_in + ATP_in*10.**(-phi_in/Z)))*1./(1. + KmADP/ADP)


def r_H2PO4(H):
    return 1./(1. + 10.**(getPH(H)-6.8)) 

def v_Pi(Pi, Pi_in, H, H_in):
    H2PO4 = Pi*r_H2PO4(H)
    H2PO4_in = Pi_in*r_H2PO4(H_in)
    return kPi*(H2PO4*H - H2PO4_in*H_in)

def v_CK(ATP, ADP, PCr, Cr, H):
    return kfCK * ADP * PCr* H - kbCK * ATP * Cr

def v_UT(ATP):
    return kUT/(1. + KmA/ATP)
    #return kUT2*ATP


def v_AK(ATP, ADP, AMP):
    return kAK * ADP_f(ADP) * ADP_m(ADP) - kbAK * ATP_m(ATP) * AMP


def v_glyco(ADP, H):
    return kGlyco*ADP*x0[8]/H

def v_lac(pyruvate, lactate, H):
    return kf_lac*pyruvate*H - kb_lac*lactate

def v_lac_eff(lactate):
    return k_lac_eff*(lactate - 1.)

def v_AD(ADP, AMP, Pi, H, ATPuse): # nimmt noch ein H weg!!! siehe 00KONSTANTEN.txt # inhibitor ATP und activator ADP, paper: AMPdesaminase_ADP_activation_ATP_inhibitor
    myosin_factor = 1. + ATPuse/30.
    if getPH(H) > 6.5 and getPH(H) <= 7.5:
        return kAD * AMP /(AMP + AMPd_K*Pi/x0[4]/ADP*x0[1]/myosin_factor) * (7.5-getPH(H))
    else:
        return kAD * AMP /(AMP + AMPd_K*Pi/x0[4]/ADP*x0[1]/myosin_factor) * (7.5-6.5)

def v_AS(ATP, AMP, IMP):     # effektiv: IMP + ATP -> AMP + ADP + Pi (eigentlich GTP->GDP), inhibiert durch AMP
    GTP = 0.148*ATP  # traut1994physiological
    return kAS*GTP * IMP / (1. + GTP/Km_GTP) / (1. + IMP/Km_IMP + AMP/Km_AMP)

def v_PL(IMP):
    return kPL * (IMP-x0[11])

def v_IS(ATP, ADP, AMP, IMP):
    if (x0[0]+x0[1]+x0[2])+x0[11]-ATP-ADP-AMP-IMP>1.e-3:
        return kIS * ATP / (1. + ATP/KmATP_IS + ADP/KiADP_IS) / (1. + AMP/KiAMP_IS + IMP/KiIMP_IS)  
    return 0.

def printFuxes():
    print "v_DH: ", v_DH(x0[18], x0[6])
    print "v_C1: ", v_C1(x0[18], x0[19], x0[8], x0[7])
    print "v_C3: ", v_C3(x0[20], x0[19], x0[8], x0[7])
    print "v_C4 (*2): ", v_C4(x0[20], x0[8], x0[7])*2.
    print "v_SN: ", v_SN(x0[3], x0[10], x0[5], x0[8], x0[7])
    print "v_EX: ", v_EX(x0[0], x0[1], x0[3], x0[10], x0[8], x0[7])
    print "v_Pi: ", v_Pi(x0[4], x0[5], x0[8], x0[7])
    print "v_UT: ", v_UT(x0[0])
    print "v_glyco+v_SN: ", 1.5*v_glyco(x0[1], x0[8]) + v_SN(x0[3], x0[10], x0[5], x0[8], x0[7])
    print "v_AK: ", v_AK(x0[0], x0[1], x0[2])
    print "v_CK: ", v_CK(x0[0],x0[1], x0[13], x0[12], x0[8])
    print "v_LK: ", v_LK(x0[8], x0[7])
    print "v_AS: ", v_AS(x0[0], x0[2], x0[11])
    print "v_AD: ", v_AD(x0[1], x0[2], x0[4], x0[8], 1.)
    print "v_IS: ", v_IS(x0[0], x0[1], x0[2], x0[11])

    print "v_SN/v_glyco: ", v_SN(x0[3], x0[10], x0[5], x0[8], x0[7])/v_glyco(x0[1], x0[8])
    print "v_DH/v_glyco: ", v_DH(x0[18], x0[6])/v_glyco(x0[1], x0[8])

def systemala(x, t):
    global t_global, ATPusageFactor, countala
    countala += 1
    dxdt = [0. for i in range(spezies_num)]
    if do_exercise:
        transitionTime = 5./60.
        transitionTime_production = 25./60.
        if t<exercise_times[0]:
            ATPusageFactor = ATPusageFactorStart
            production_factor = ATPusageFactor
        elif t<exercise_times[1]:
            ATPusageFactor = ATPusageFactorStart - (ATPusageFactorStart - ATPusageFactorMid) * (t-exercise_times[0])/(exercise_times[1]-exercise_times[0])
            production_factor = ATPusageFactor
        elif t > exercise_times[1]:
            ATPusageFactor = ATPusageFactorEnd + (ATPusageFactorMid - ATPusageFactorEnd) * np.exp(-(t-exercise_times[1])/transitionTime)
            production_factor = ATPusageFactorEnd + (ATPusageFactorMid - ATPusageFactorEnd) * np.exp(-(t-exercise_times[1])/transitionTime_production)

    else:
        ATPusageFactor = 1.
        production_factor = 1.
    glycoFactor = (production_factor)**1.15 
    oxPhosFactor = (production_factor)**0.62*CFS_factor 
        
    vs_DH = v_DH(x[18], x[6]) * oxPhosFactor
    vs_C1 = v_C1(x[18], x[19], x[8], x[7]) * oxPhosFactor
    vs_C3 = v_C3(x[20], x[19], x[8], x[7]) * oxPhosFactor
    vs_C4 = v_C4(x[20], x[8], x[7]) * oxPhosFactor
    vs_SN = v_SN(x[3], x[10], x[5], x[8], x[7]) * oxPhosFactor
    vs_EX = v_EX(x[0], x[1], x[3], x[10], x[8], x[7]) * oxPhosFactor
    vs_Pi = v_Pi(x[4], x[5], x[8], x[7]) * oxPhosFactor
    vs_LK = v_LK(x[8], x[7])
    vs_glyco = v_glyco(x[1], x[8]) * glycoFactor
    vs_AK = v_AK(x[0], x[1], x[2])
    vs_CK = v_CK(x[0],x[1], x[13], x[12], x[8])# *0.01
    vs_EFF = v_EFF(x[8])
    vs_lac = v_lac(x[6], x[9], x[8])
    vs_lac_eff = v_lac_eff(x[9])
    vs_AS = v_AS(x[0], x[2], x[11])
    vs_AD = v_AD(x[1], x[2], x[4], x[8], ATPusageFactor)
    vs_IS = v_IS(x[0], x[1], x[2], x[11])
    vs_PL = 0.
    if ATPusageFactor>1. and t<0.5:
        vs_PL = v_PL(x[11])
    vs_UT = v_UT(x[0]) * ATPusageFactor
    
    dxdt[0] = vs_EX - vs_UT + vs_AK + vs_CK + 1.5*vs_glyco - vs_AS - 5.*vs_IS # - vs_PS
    dxdt[1] = - vs_EX + vs_UT - 2*vs_AK - vs_CK - 1.5*vs_glyco + vs_AS + 4.*vs_IS
    dxdt[2] = vs_AK + vs_AS - vs_AD + vs_IS # vs_PS  # vs_IS
    dxdt[3] = 15.*(vs_SN - vs_EX)
    dxdt[4] = vs_UT - vs_Pi - 1.5*vs_glyco + vs_AS + 6.*vs_IS # + 2.*vs_PS
    dxdt[5] = 15.*(vs_Pi - vs_SN)
    dxdt[6] = vs_glyco - vs_DH/NADH_per_pyruvate - vs_lac
    dxdt[7] =  - 15.*(2.*(2.+2.*u)*vs_C4 + (4.-2.*u)*vs_C3 + 4.*vs_C1 - nA*vs_SN - u*vs_EX - (1.-u)*vs_Pi -vs_LK - (1.-r_H2PO4(x[7]))*vs_Pi + (1.-r_H2PO4(x[7]))*vs_SN)
    dxdt[8] = 2.*(2.+2.*u)*vs_C4 + (4.-2.*u)*vs_C3 + 4.*vs_C1 - nA*vs_SN - u*vs_EX - (1.-u)*vs_Pi -vs_LK - vs_CK - vs_EFF + (0.5+1.5*r_H2PO4(x[8]))*(vs_glyco-vs_DH/NADH_per_pyruvate) - (1.-r_H2PO4(x[8]))*vs_Pi + (1.-r_H2PO4(x[8]))*(vs_UT+6.*vs_IS) +(2.-r_H2PO4(x[8]))*vs_AS- vs_lac - vs_AD # eigentlich +(0.5+1.5*r_H2PO4(x[8]))*(vs_glyco -vs_DH/5.6)
    dxdt[9] = vs_lac - vs_lac_eff
    dxdt[10] = - dxdt[3]
    dxdt[11] = - vs_AS + vs_AD - vs_PL + vs_IS
    dxdt[12] = vs_CK
    dxdt[13] = - vs_CK
    dxdt[14] = vs_PL
    dxdt[15] = 0. # NOT USED right now
    dxdt[16] = vs_EFF
    dxdt[17] = vs_lac_eff
    dxdt[18] = 15.*(vs_DH - vs_C1)/5. # /5. ist buffering capacity for NAD
    dxdt[19] = 15.*(vs_C1 - vs_C3)
    dxdt[20] = 15.*(vs_C3 - 2.*vs_C4)
    
    r_buffi = 22./np.log(10.)/x[7]
    r_buffe = 25./np.log(10.)/x[8]

    dxdt[7] /= r_buffi
    dxdt[8] /= r_buffe
    dxdt[16] /= r_buffe
    if t>t_global+deltaT_flxues:
        t_global = t
        t_fluxes.append(t)
        fluxes[0].append(-vs_UT) # ATPuse
        fluxes[1].append(vs_EX) # oxPhos
        fluxes[2].append(1.5*vs_DH/NADH_per_pyruvate) # vs_DH/5. # aerobic glycolysis
        fluxes[3].append(1.5*vs_glyco) # glycolysis
        fluxes[4].append(vs_CK) #append(creatine_ATP) # creatine
        fluxes[5].append(vs_AK) # adenylate kinase
        fluxes[6].append(-vs_AS) # AMPsynthesis
        fluxes[7].append(-5.*vs_IS) # IMPsynthesis
        ATPusageFactor_array[0].append(ATPusageFactor)
        ATPusageFactor_array[1].append(production_factor)
        v_O2_array.append(vs_C4)
    return dxdt

flux_0 = []
for i in range(len(x0)):
    if x0[i] != 0.:
        flux_0.append(systemala(x0, 0.)[i]/x0[i])
    else:
        flux_0.append(systemala(x0, 0.)[i])
print flux_0

def simy(time_start, time_stop, deltaT):
    global x
    t = [time_start + deltaT*i for i in range(int((time_stop-time_start)/deltaT))]
    Y = map(list, zip(*itg.odeint(systemala, x, t, mxstep=2000, h0=1.e-7))) 
    x = [Y[i][-1] for i in range(spezies_num)]
    print "xEnd:", x
    return Y, t

def simyExercise(exercise_timala):
    global ATPusageFactorStart, ATPusageFactorMid, ATPusageFactorEnd, exercise_times, do_exercise
    exercise_times = exercise_timala
    ATPusageFactorExercise = max_atp_usage
    ATPusageFactorStart = ATPusageFactorExercise
    if ATPusageFactorExercise >= 1.7:
        ATPusageFactorMid = ATPusageFactorExercise*0.6
    else:
        ATPusageFactorMid = 1.
    ATPusageFactorEnd = 1.
    Y, t = simy(0., exercise_timala[2], 0.005)
    return Y, t
