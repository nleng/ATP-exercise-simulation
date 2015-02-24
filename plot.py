#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import ATPmodel as k

dpi = 150
lw = 1.8

k.printFuxes()

def plotta(t_plt, x_plt):
    linestyles = ['-', ':', '-.', ':', ':', '5--', '6-.', '7-', '8-', '--', '10-.', '--', '-.', '-', '14-.', '15:', '16-', '17--']
    colors = ['blue', 'purple', 'red', 'purple', 'purple', '5', '6', '7', '8', 'green', '10', 'green', 'red', 'blue', '14', '15', '16', '17']
    labels = ["ATP [mM]", r"ADP/ADP$_0$", r"AMP/AMP$_0$/10", "ATP_in", "Pi", "Pi_in", "pyruvat", "pH mito", "pH cyto", "lactate", "ADP_in", "IMP [mM]", "Cr", "PCr", "purine loss", "PRPP", "blood pH", "blood lactate", "NADH/NAD", "UQH2", "c2"]
    ATPlabels = ["ATP [mM]", r"ADP/ADP$_0$", r"AMP/AMP$_0$/10", "IMP [mM]"]
    inLabels = ["ATP_in", "ADP_in", "Pi_in", "pH mito", "pH cyto", "pH test", "NADH/NAD", "UQH2", "c2"]
    out_label = ["Cr", "PCr", "Pi", "lactate"]
    
    takeEvery = 1
    
    print "plotting..."
    
    print 'deltaT (original): ', t_plt[1], t_plt[100]
    t_plt = t_plt[0::takeEvery]
    print 'deltaT: ', t_plt[1], t_plt[10]
    for i in range(len(x_plt)):
        x_plt[i] = x_plt[i][0::takeEvery]
    for i in range(len(x_plt[1])):
        x_plt[1][i] /= k.x0[1]
    for i in range(len(x_plt[2])):
        x_plt[2][i] /= k.x0[2]*10.
    for i in range(len(x_plt[8])):
        x_plt[8][i] = k.getPH(x_plt[8][i])
    for i in range(len(x_plt[7])):
        x_plt[7][i] = k.getPH(x_plt[7][i])
    
    for i in range(len(x_plt[18])):
        x_plt[18][i] /= (k.N_t-x_plt[18][i])
    for i in range(len(x_plt[19])):
        x_plt[19][i] /= k.Q_t
    for i in range(len(x_plt[20])):
        x_plt[20][i] /= k.c_t
    
    print "ATP: min, at 30s \t", min(x_plt[0]), x_plt[0][100]
    print "ADP: max , at 30s \t", max(x_plt[1]), x_plt[1][100]
    print "AMP: max , at 30s \t", 10.*max(x_plt[2]), 10.*x_plt[2][100]
    print "IMP: max , at 30s \t", max(x_plt[11]), x_plt[11][100]
    print "end ATP: ", x_plt[0][-1]
    print "end PRPP: ", x_plt[15][-1]
    print "end IMP: ", x_plt[11][-1]
    print "purine loss: ", x_plt[14][-1]
    print "PH: min, at 30s \t\t\t", min(x_plt[8]), x_plt[8][100], x_plt[8][200]
    print "lactate+pyr: max, at 30s, 6.5min \t", max(x_plt[9][0:600])+max(x_plt[6]), x_plt[9][100]+x_plt[6][100]
    print "lactate efflux: total, 30s, 60s \t", x_plt[17][-1], x_plt[17][100], x_plt[17][200]
    print "H efflux: ", x_plt[16][-1]
    print "min PCr: ", min(x_plt[13])
    print "max Cr: ", max(x_plt[12])
    
    maxaArray = []
    while len(maxaArray)<k.spezies_num:
        maxama = 0.0
        maxaI = 0
        for i in range(k.spezies_num):
            if x_plt[i][0] >= maxama and i not in maxaArray:
                maxama = x_plt[i][0]
                maxaI = i
        maxaArray.append(maxaI)
        
    #---------- plot inLabels ----------#
    fig = plt.figure(figsize=(6, 4.5))
    for i in maxaArray:
        if labels[i] in inLabels:
            plt.plot(t_plt, x_plt[i], label = labels[i], linewidth=lw)
    plt.ylabel("concentration [mM]", fontsize=20)
    plt.xlabel("time [min]", fontsize=20)
    plt.legend()
    plt.tight_layout()
    plt.savefig("cons_in.png", dpi=300)
    
    #---------- plot the rest ----------#
    fig = plt.figure(figsize=(6, 4.5))
    for i in maxaArray:
        if labels[i] not in ATPlabels and labels[i] not in inLabels:
            plt.plot(t_plt, x_plt[i], label = labels[i], linewidth=lw)
    plt.ylabel("concentration [mM]", fontsize=20)
    plt.xlabel("time [min]", fontsize=20)
    plt.legend()
    plt.tight_layout()
    plt.savefig("cons_out2.png", dpi=300)
    
    #---------- plot ATPlabels ----------#
    fig = plt.figure(figsize=(6, 4.5))
    for i in maxaArray:
        if labels[i] in ATPlabels:
            plt.plot(t_plt, x_plt[i], label = labels[i], linewidth=lw, linestyle=linestyles[i], color=colors[i])
    leg = plt.legend(loc=1, labelspacing=0.1, borderpad=0.1, borderaxespad=0.2)
    leg.get_frame().set_linewidth(0.0)
    for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)
    plt.xlim(0., 1.5)
    plt.ylim(0., 6.15)
    plt.xticks([0, 0.5, 1., 1.5])
    plt.yticks([0, 2, 4, 6])
    plt.ylabel("concentration")
    plt.xlabel("time [min]")
    plt.tight_layout()
    plt.savefig("cons_atp.png", dpi=dpi)
    
    #---------- plot out_label ----------#
    fig = plt.figure(figsize=(6, 4.5))
    for i in maxaArray:
        if labels[i] in out_label:
            plt.plot(t_plt, x_plt[i], label = labels[i], linewidth=lw, linestyle=linestyles[i], color=colors[i])
    leg = plt.legend(loc=1, labelspacing=0.1, borderpad=0.1, borderaxespad=0.2)
    leg.get_frame().set_linewidth(0.0)
    for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)
    plt.xlim(0., 3.)
    plt.ylim(0., 33.)
    plt.xticks([0, 1, 2, 3])
    plt.yticks([0, 10, 20, 30])
    plt.ylabel("concentration [mM]")
    plt.xlabel("time [min]")
    plt.tight_layout()
    plt.savefig("cons_out.png", dpi=dpi)
    
    #---------- integrate fluxes over 30 sec ----------#
    oxPhos_integrated, glyco_integrated, glyco_anaerobic_integrated, CK_integrated = 0., 0., 0., 0.
    for i in range(1, len(k.t_fluxes)):
        dT = k.t_fluxes[i]-k.t_fluxes[i-1]
        oxPhos_integrated += k.fluxes[1][i]*dT
        glyco_integrated += k.fluxes[2][i]*dT
        glyco_anaerobic_integrated +=k.fluxes[3][i]*dT
        CK_integrated += k.fluxes[4][i]*dT
        if k.t_fluxes[i]>0.5:
            print i, k.t_fluxes[i]
            break
    print "end exercise t:", k.t_fluxes[i], k.fluxes[4][i], k.fluxes[4][i+1], k.fluxes[4][i+2], "FLUXES: ", oxPhos_integrated, glyco_integrated, glyco_anaerobic_integrated, CK_integrated
    
    #---------- flxues plot normal ----------#
    fluxLabel = ["ATP usage", "oxPhos", "aerobic glycolysis", "glycolysis", "creatine kinase", "adenylate kinase", "AMP synthesis", "IMP synthesis"]
    print "max ATP usage: ", min(k.fluxes[0])
    print "mean ATP usage: ", 0.8*min(k.fluxes[0])
    plt.figure("Fluxes")
    for i in range(7):
        plt.plot(k.t_fluxes, k.fluxes[i], label = fluxLabel[i], lw=lw)
    plt.legend()
    plt.xlim(0., 1.5)
    maxRatio = 0.
    for i in range(len(k.fluxes[0])):
        tmpRatio = (k.fluxes[2][i]+k.fluxes[3][i])/k.fluxes[1][i]
        if tmpRatio > maxRatio:
            maxRatio = tmpRatio
    print "maximum glyco/oxPhos ratio: ", maxRatio

    #---------- fluxes plot as area ----------#
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    fluxLabelsOut = ["ATP usage", "AMP synthesis", "IMP synthesis", "", ""]
    fluxColorsOut = ["red", "cyan", "olive", "blue", "darkorange"]
    fluxLabelsOut = ["ATP usage","", ""]
    fluxColorsOut = ["red", "blue", "darkorange"]
    fluxLabelsIn = ["oxPhos", "glycolysis", "creatine kinase", "adenylate kinase"]
    fluxColorsIn = ["green", "purple", "blue", "darkorange"] # , "blueviolet"
    CK_exercise, CK_recovery, AK_exercise, AK_recovery = [], [], [], []
    for bb in range(len(k.fluxes[4])):
        if k.fluxes[4][bb] > 0.:
            CK_exercise.append(k.fluxes[4][bb])
            CK_recovery.append(0.)
        else:
            CK_exercise.append(0.)
            CK_recovery.append(k.fluxes[4][bb])
    for bb in range(len(k.fluxes[5])):
        if k.fluxes[5][bb] > 0.:
            AK_exercise.append(k.fluxes[5][bb])
            AK_recovery.append(0.)
        else:
            AK_exercise.append(0.) 
            AK_recovery.append(k.fluxes[5][bb])
    yOut = [k.fluxes[0], k.fluxes[6], k.fluxes[7], CK_recovery, AK_recovery]
    yOut = [k.fluxes[0], CK_recovery, AK_recovery]
    yIn = [k.fluxes[1], k.fluxes[3], CK_exercise, AK_exercise]
    print min(k.fluxes[5])
    for i in range(1, len(yOut)):
        for j in range(len(yOut[i])):
            yOut[i][j] += yOut[i-1][j]


    
    for i in range(1, len(yIn)):
        for j in range(len(yIn[i])):
            yIn[i][j] += yIn[i-1][j]
    y0 = [0. for i in range(len(yIn[0]))]
    for i in reversed(range(1, len(yIn))):
        ax.fill_between(k.t_fluxes, yIn[i], yIn[i-1], facecolor=fluxColorsIn[i], alpha=0.5)
        plt.plot(k.t_fluxes, yIn[i], label = fluxLabelsIn[i], color=fluxColorsIn[i], linewidth=0., alpha=0.5)
    ax.fill_between(k.t_fluxes, yIn[0], y0, facecolor=fluxColorsIn[0], alpha=0.5)
    plt.plot(k.t_fluxes, yIn[0], label = fluxLabelsIn[0], color=fluxColorsIn[0], linewidth=0., alpha=0.5)
    
    y0 = [0. for i in range(len(yOut[0]))]
    for i in reversed(range(1, len(yOut))):
        ax.fill_between(k.t_fluxes, yOut[i], yOut[i-1], facecolor=fluxColorsOut[i], alpha=0.5)
        plt.plot(k.t_fluxes, yOut[i], label = fluxLabelsOut[i], color=fluxColorsOut[i], linewidth=0., alpha=0.5)
    ax.fill_between(k.t_fluxes, yOut[0], y0, facecolor=fluxColorsOut[0], alpha=0.5)
    plt.plot(k.t_fluxes, yOut[0], label = fluxLabelsOut[0], color=fluxColorsOut[0], linewidth=0., alpha=0.5)

    leg =plt.legend(loc=1, labelspacing=0.1, borderpad=0.1, borderaxespad=0.2)
    leg.get_frame().set_linewidth(0.0)
    for legobj in leg.legendHandles:
            legobj.set_linewidth(5.0)
    plt.xticks([0., 0.5, 1., 1.5])
    plt.xlabel("Time [min]")
    ytickas = [-180., -90., 0., 90., 180.]
    plt.yticks(ytickas, [round(yt/60.,2) for yt in ytickas])
    plt.ylabel("Flux [mM/s]")

    plt.xlim(0., 1.5)
    plt.ylim(-205, 205)
    #plt.ylim(-5, 155)
    plt.tight_layout()
    plt.savefig("fluxes.png", dpi=dpi)
    
    #------- pH value plot -------#
    fig = plt.figure(figsize=(6, 4))
    plt.plot(t_plt, x_plt[8], lw=2., color="black")
    plt.xlim(0., 10.)
    plt.ylim(6.3, 7.)
    plt.xticks([0, 1, 2,4,6,8, 10])
    plt.yticks([6.4, 6.6, 6.8, 7.])
    plt.ylabel("pH")
    plt.xlabel("time [min]")
    plt.tight_layout()
    plt.savefig("pH_time.png", dpi=dpi)
    
    #------- ATPusageFactor plot -------#
    fig = plt.figure(figsize=(6, 4))
    plt.plot(k.t_fluxes, k.ATPusageFactor_array[0], lw=2., label="ATP consumption factor")
    plt.plot(k.t_fluxes, k.ATPusageFactor_array[1], lw=2., label="ATP synthesis factor")
    leg = plt.legend(loc=1, labelspacing=0.1, borderpad=0.1, borderaxespad=0.2)
    leg.get_frame().set_linewidth(0.0)
    for legobj in leg.legendHandles:
            legobj.set_linewidth(5.0)
    plt.xlim(0., 1.5)
    plt.xticks([0., 0.5, 1., 1.5])
    plt.yticks([0, 100, 200, 300])
    plt.ylabel("amplification factor")
    plt.xlabel("time [min]")
    plt.tight_layout()
    plt.savefig("ATPusage.png", dpi=dpi)
    
    plt.show()

maximal_time = 0.
half_time = 0.5
YY, tt = k.simyExercise([maximal_time, maximal_time+half_time, maximal_time+half_time + 10.*60.]) # 90. 600. 720.
plotta(tt, YY)

