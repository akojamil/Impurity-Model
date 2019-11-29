import numpy as np
import matplotlib.pyplot as plt
import Library as Lib 
import Outgassing as Out
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import math

params = {'text.usetex':True,'font.size':12,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black']
Linestyles = ['-', '--', ':', 'dashdot','-', '--', ':', 'dashdot','-', '--', ':', 'dashdot']

class System(): 
    def __init__(self, Setup, Material, Solute, Version):
        self.Name = Setup
        self.Material = Material
        self.Solute = Solute
        self.Version = Version

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

def GetTimeStamps(Points, Spacing, TimeScale): 
    X = []
    for ii, x in enumerate(Points): 
        if ii == len(Points)-1: break 
        Time = np.linspace(Points[ii], Points[ii+1], int((Points[ii+1] - Points[ii]) / Spacing + 1))
        X.append(Time)
    return np.array(X)

def PlotImpuritiesVsTime(Data, XRange=0, YRange=0, XTicks=0):
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    plt.xlabel('Time [hours]', fontsize=16)
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
    ax.legend(loc='upper right', fontsize=10)
    fig.tight_layout()

def PlotFlowRateVsTime(Data, XRange=0, YRange=0, XTicks=0):
    print(Data)
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    plt.xlabel('Time [Hours]', fontsize=16)
    plt.ylabel(r'Outgassing Rate [mBar$\,\cdot\,$Liter/s]', fontsize=16)
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    

    for jj,data in enumerate(Data):
        for ii,(X,Y) in enumerate(zip(data.Time, data.FlowRate)):
            plt.plot(X, Y, label=data.Labels[ii], color=colors[ii], linewidth=2.0, linestyle=Linestyles[jj])

    if XTicks == 0:
        XTicks = int(np.max([np.max(x.Time) for x in Data ])/10)
    else:
        XTicks = XTicks 
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(XTicks))
    ax.grid(b=True, which='major', color='grey', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')

    if len(Data)>1: 
        YValues = [x.FlowRate for x in Data]
        YValues = [item for sublist in YValues for item in sublist]
        YValues = [item for sublist in YValues for item in sublist]
        YMaxExp = np.max(np.unique(YValues))
    else:
        YValue = [item for sublist in Data[0].FlowRate for item in sublist]
        YValue = np.ravel(YValue)
        YMaxExp = np.max(YValue)

    # YMaxExp = np.max(Data)
    YMaxExp = math.ceil(math.log10(YMaxExp))
    if XRange == 0:
        plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    else:
        plt.xlim(XRange[0], XRange[1])
    if YRange == 0:
        plt.ylim(ymin=1E-12, ymax = 1*10**YMaxExp)
    else:
        plt.ylim(YRange[0], YRange[1])

    ax.legend(loc='upper right', fontsize=10)
    fig.tight_layout()

def DoModelling(Systems, Labels, Temperature, Time, TimeScale, Constraints):
    for ii, System in enumerate(Systems): 
        # Define the different temperatures for which to calculate outgassing
        System.Temp = Temperature

        # Calculate the diffusion constants for the above defined temperatures using the Arrhenius equation
        System.DiffConstants = Out.GetDiffTemp(System, Temperatures=System.Temp)    

        # Get the initial number of impurities from model parameters. Can define '#', 'ppm,ppb,ppt' or 'Mass' for units 
        System.InitialImpurities = Out.GetInitialImpurities(System, '#')

        # Forward above defined labels, that will be later used for the plot legend.
        System.Labels = Labels[ii]

        System.Time = Time

        System.Constraints = Constraints[ii]

        # Calculate the number of impurities left in the sample vs time by using the solution to the diffusion equation.
        System.Impurities = Out.GetImpuritiesVsTime(Data=System, TimeScale=TimeScale)

        # Calculate the outgassing rate vs time using Fick's 1st law 
        System.FlowRate = Out.GetFlowRateVsTime(Data=System, Units='#', TimeScale=TimeScale)

        

        # This will output all important parameters that have been imported into each model
        # System.Print()

def GetLabels(Systems, Temperature): 
    Labels = []
    for ii, System in enumerate(Systems): 
        label = [] 
        for jj, temp in enumerate(Temperature):
            first = '%s, $ d= %s \, \mathrm{cm}, T = %d \, \mathrm{K}, \, E_a = %s\,\mathrm{eV}$' % (System.Version, System.Thickness, temp, System.ActivationEnergy)
            label.append(first)
        Labels.append(label)
    return Labels 


if __name__ == '__main__':
    # Define parameters for each model. Details about which options are available are in Library.py
    # Units for all parameters are defined in Library.py.
    # S1 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='EXO-Teflon')
    S2 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='Stock-Teflon')
    S3 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='Stock-Teflon Thick')
    S3 = System(Setup='EXO-200', Material='Teflon', Solute='Oxygen', Version='EXO-Teflon')
    S2 = System(Setup='EXO-200', Material='Teflon Yale', Solute='Oxygen', Version='EXO-Teflon')
    
    # Bundle above defined models together to do calculations for all of them 
    Systems = [S3, S2]

    # Define the different temperatures for which to calculate outgassing
    Temperature = [290, 293.2, 164]

    # Labels = [['EXO-200 Teflon', 'EXO-200 Teflon2', 'EXO-200 Teflon3'], ['Stock Room Teflon'], ['Columbia Setup Teflon']]
    Labels = [['EXO-200 Teflon', 'EXO-200 Teflon2', 'EXO-200 Teflon3']]

    # Define time scale over which to calculate and plot the outgassing rate and impurity concentrations
    TimeScale = 'Days'

    # Sets time regions for the above defined temperature values
    # Points goes from 0 to first value, first value to second value and so on. Ex: [0,10,15] gives two regions, 0-10 and 10-15
    # TimeScale defines in which units to plot the x-Axis later.
    Time = GetTimeStamps(Points=[0,1,2,10], Spacing=0.01, TimeScale=TimeScale)

    Labels = GetLabels(Systems, Temperature)

    Constraints = [[1000,10,100], [1000,10,100]]

    # Executing all steps to get impuritie numbers and outgassing rate as a function of time 
    DoModelling(Systems, Labels, Temperature, Time, TimeScale, Constraints)

    # Plotting the impurities and outgassing rates for all above defined models and saving the plots. 
    PlotImpuritiesVsTime(Systems, XRange=[0,10], YRange=[4E14, 5E20], XTicks=1)
    # plt.savefig('impurities.pdf')

    # Systems.Time /= 24.0
    PlotFlowRateVsTime(Systems, XRange=[0,10], YRange=[1E10, 1E21], XTicks=1)
    # plt.savefig('outgassing_rate.pdf')
    plt.show()
