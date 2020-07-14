import numpy as np
import matplotlib.pyplot as plt
import Library as Lib 
import Outgassing as Out
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import math

params = {'text.usetex':True,'font.size':16,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black']
Linestyles = ['-', '--', ':', 'dashdot','-', '--', ':', 'dashdot','-', '--', ':', 'dashdot']

class System(): 
    def __init__(self, Setup, Material, Solute, Version):
        self.Name = Setup
        self.Material = Material
        self.Solute = Solute
        self.Version = Version

        self.TimeScale = 'Day'

        self.Constraints = []
        self.ConstraintIndex = []
        
        self.Diffusion = Lib.Material.get(Material).get(Solute).get('Diffusion Constant')
        self.Solubility = Lib.Material.get(Material).get(Solute).get('Solubility')
        self.ActivationEnergy = Lib.Material.get(Material).get(Solute).get('Activation Energy')
        self.Abundance =  Lib.Gas.get(Solute).get('Abundance in Air')
        self.MolarMass =  Lib.Gas.get(Solute).get('Molar Mass')

        self.XeMass = Lib.System.get(Setup).get('Xenon Mass')
        self.Volume = Lib.System.get(Setup).get(Material).get(Version).get('Volume')
        self.Area = Lib.System.get(Setup).get(Material).get(Version).get('Area')
        self.Thickness = Lib.System.get(Setup).get(Material).get(Version).get('Thickness')

    def Print(self):
        Attributes = vars(self)
        print('\n '.join("%s: %s" % item for item in Attributes.items()))

def GetTimeStamps(Points, Spacing): 
    X = []
    for ii, x in enumerate(Points): 
        if ii == len(Points)-1: break 
        Time = np.linspace(Points[ii], Points[ii+1], int((Points[ii+1] - Points[ii]) / Spacing + 1))
        X.append(Time)
    return X

def PlotImpuritiesVsTime(Data, XRange=0, YRange=0, XTicks=0, TimeScale='Hour', Size=(8,6)):
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    plt.xlabel('Time [%s]' % TimeScale, fontsize=16)
    plt.ylabel(r'Total Number of Impurities', fontsize=16)
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    
    for jj, data in enumerate(Data):
        for ii,(X,Y) in enumerate(zip(data.Time, data.Impurities)):
            plt.plot(X, Y, label=data.Labels[ii], color=colors[ii], linewidth=2.0, linestyle=Linestyles[jj])

    if XTicks == 0:
        XTicks = int(np.max([np.max(x.Time) for x in Data ])/10)
    else:
        XTicks = XTicks 
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(XTicks))
    ax.grid(b=True, which='major', color='grey', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    
    if XRange == 0:
        plt.xlim(0,10)
        # plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    else:
        plt.xlim(XRange[0], XRange[1])
    if YRange == 0:
        YMaxExp = np.max([np.max(x.Impurities) for x in Data ])
        YMaxExp = math.ceil(math.log10(YMaxExp))
        plt.ylim(ymin=1E10, ymax = np.power(10.0,YMaxExp))
    else:
        plt.ylim(YRange[0], YRange[1])
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()

def PlotFlowRateVsTime(Data, XRange=0, YRange=0, XTicks=0, TimeScale='Hour', Size=(8,6)):
    fig = plt.figure(figsize=Size)
    ax = fig.gca()
    plt.xlabel('Time [%s]' % TimeScale, fontsize=16)
    plt.ylabel(r'Outgassing Rate [mBar$\,\cdot\,$Liter/s]', fontsize=16)
    # plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    
    for jj,data in enumerate(Data):
        # for ii,(X,Y) in enumerate(zip(data.Time, data.FlowRate)):
        #     print(ii,X,Y)

        plt.plot(data.Time, data.FlowRate, label='Data', color=colors[jj], linewidth=2.0, linestyle=Linestyles[jj])

    if XTicks == 0:
        XTicks = int(np.max([np.max(x.Time) for x in Data ])/10)
    else:
        XTicks = XTicks 
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(XTicks))
    ax.grid(b=True, which='major', color='grey', alpha=0.3, linestyle='--')
    ax.grid(b=True, which='minor', color='grey', alpha=0.3, linestyle=':')

    # if len(Data)>1: 
    #     YValues = [x.FlowRate for x in Data]
    #     YValues = [item for sublist in YValues for item in sublist]
    #     YValues = [item for sublist in YValues for item in sublist]
    #     YMaxExp = np.max(np.unique(YValues))
    # else:
    #     YValue = [item for sublist in Data[0].FlowRate for item in sublist]
    #     YValue = np.ravel(YValue)
    #     YMaxExp = np.max(YValue)

    # YMaxExp = math.ceil(math.log10(YMaxExp))
    if XRange == 0:
        plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    else:
        plt.xlim(XRange[0], XRange[1])
    if YRange == 0:
        plt.ylim(ymin=1E-12, ymax = 1*10**YMaxExp)
    else:
        plt.ylim(YRange[0], YRange[1])

    ax.legend(loc='upper right', fontsize=12, framealpha=1.0)
    fig.tight_layout()

def DoModelling(Systems, Labels, Temperature, Time, TimeScale, Constraints):
    for ii, System in enumerate(Systems): 
        # Define the different temperatures for which to calculate outgassing
        System.Temp = Temperature

        # Calculate the diffusion constants for the above defined temperatures using the Arrhenius equation
        System.DiffConstants = Out.GetDiffTemp(System, Temperatures=System.Temp, TimeScale=TimeScale)  
        fig = plt.figure(figsize=(8,6))  
        plt.plot(Time, System.DiffConstants)
        plt.xlabel('Time [Day]')
        plt.ylabel('Diff [cm$^2$/%s]' % TimeScale)
        plt.yscale('log')
        plt.tight_layout()
        # plt.savefig('temp.pdf')

        # Get the initial number of impurities from model parameters. Can define '#', 'ppm,ppb,ppt' or 'Mass' for units 
        System.InitialImpurities = Out.GetInitialImpurities(System, '#')

        # Forward above defined labels, that will be later used for the plot legend.
        System.Labels = Labels[ii]

        System.Time = Time

        System.Constraints = Constraints[ii]

        # Calculate the number of impurities left in the sample vs time by using the solution to the diffusion equation.
        System.Impurities = Out.GetImpuritiesVsTime(Data=System, TimeScale=TimeScale)
       
        # Calculate the outgassing rate vs time using Fick's 1st law 
        System.FlowRate = Out.GetFlowRateVsTime(Data=System, Units='mBar Liter', TimeScale=TimeScale)

        # This will output all important parameters that have been imported into each model
        # System.Print()

def GetLabels(Systems, Temperature): 
    Labels = []
    for ii, System in enumerate(Systems): 
        label = [] 
        for jj, temp in enumerate(Temperature):
            first = '%s, $ d= %s \, \mathrm{cm}, \, T = %d \, \mathrm{K}, \, E_a = %s\,\mathrm{eV}$' % (System.Version, System.Thickness, temp[0], System.ActivationEnergy)
            label.append(first)
        Labels.append(label)
    return Labels 


if __name__ == '__main__':
    # Define parameters for each model. Details about which options are available are in Library.py
    # Units for all parameters are defined in Library.py.

    # S1 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='EXO-Teflon')
    S1 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='Stock-Teflon Thick')
    S1 = System(Setup='EXO-200', Material='Teflon', Solute='Oxygen', Version='EXO-Teflon')

    
    # Bundle above defined models together to do calculations for all of them 
    Systems = [S1]
    Time = GetTimeStamps(Points=[0,10,100,1000], Spacing=0.01)

    Time = [0, 12, 24, 48, 288, 768]
    Temperature = [293.15, 293.15, 200, 200, 164, 164]
    TimeScale = 'Hour'

    Time = [0, 1, 2, 3, 11, 31]
    Temperature = [293.15, 293.15, 200, 200, 164, 164]
    TimeScale = 'Day'

    # Time = [x * 3600*24 for x in [0, 1, 2, 4, 5, 6]]
    # Temperature = [293.15, 293.15, 200, 200, 164, 164]
    # TimeScale = 'Second'

    Temp = []
    TTime = []
    import scipy.interpolate
    for ii in range(len(Temperature)-1):
        X = [Time[ii],Time[ii+1]]
        Y = [Temperature[ii],Temperature[ii+1]]
        y_interp = scipy.interpolate.interp1d(X,Y)
        x_interp = np.arange(Time[ii],Time[ii+1],0.01)
        TTime.append(x_interp)
        Temp.append(y_interp(x_interp))
  

    Time = np.array([item for sublist in TTime for item in sublist])
    Temperature = np.array([item for sublist in Temp for item in sublist])
    fig = plt.figure(figsize=(8,6))
    plt.plot(Time, Temperature)
    plt.xlabel('Time [Day]')
    plt.ylabel('Temperature [K]')
    plt.tight_layout()

    Labels = [['EXO-200 Teflon', 'EXO-200 Teflon2', 'EXO-200 Teflon3']]

    

    # Labels = GetLabels(Systems, Temperature)

    Constraints = [[1E6,1E6,1E6], [1E6,1E6,1E6], [1E6,1E6,1E6], [1E6,1E6,1E6]]

    # Executing all steps to get impuritie numbers and outgassing rate as a function of time 
    DoModelling(Systems, Labels, Temperature, Time, TimeScale, Constraints)

    PlotFlowRateVsTime(Systems, XRange=[0,300], YRange=[1E-25, 1E0], XTicks=50, TimeScale=TimeScale, Size=(10,6))
    # Xval = np.linspace(0,1000,100)
    # plt.plot(Xval, [6.8E-8]*len(Xval), label='EXO-200 Steady-State Leak Rate', color='k', linewidth=2.0, linestyle='-')
    # plt.legend(loc='lower right')
    plt.savefig('outgassing_rate_new.pdf')
    import tikzplotlib
    width = 6.50127*2.54
    width = 6.0*2.54
    tikzplotlib.save("outgassing_rate_new.tex", figurewidth='%.2fcm' % width, figureheight='%.2fcm' % (width/2.0) )
    # plt.show()