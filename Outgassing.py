import numpy as np 
import matplotlib.pyplot as plt


IdealGasMolarVolume = 22.4 # liter/mol
Avogadro = 6.022E23 # particles/mol
Boltzmann = 8.6173303E-5 # eV/K
BoltzmannL = 83.14 # mBar*Liter/(mol*K)

def DoTimeConversion(TimeScale): 
    if TimeScale == 'Second': return 1
    if TimeScale == 'Minute': return 60.0
    if TimeScale == 'Hour': return 3600.0
    if TimeScale == 'Day': return 3600.0*24.0
    if TimeScale == 'Week': return 3600.0*24.0*7.0

def GetImpuritiesVsTime(Data, TimeScale='Hour'):
    Impurities = []
    InitialImpurities = Data.InitialImpurities
    NewImpurities = []
    for ii, (Time,Diff) in enumerate(zip(Data.Time, Data.DiffConstants)):
        Impurities.append(SolveDiffusionEquation(Time, Diff, Data.Thickness, InitialImpurities))
        InitialImpurities = Impurities[-1]
        NewImpurities.append(InitialImpurities)
    # Impurities = SolveDiffusionEquation(Data.Time, Data.DiffConstants, Data.Thickness,Data.InitialImpurities)
    # Data.InitialImpurities = 
    return np.array(Impurities)

def SolveDiffusionEquation(Time, Diff, Thickness, Conc): 
    value = 0.0 
    for N in range(0,20,1): 
        factor1 = 1.0/((2.0*N+1.0)**2)
        factor2 = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*Time)
        comp = factor1*factor2 * (Conc*8.0*Thickness)/(2.0*np.pi**2)
        value = value + comp
    return value 

def GetFlowRateVsTime(Data, Units='#', TimeScale='Hour'): 
    # print(Units)
    # print(TimeScale)
    FlowRate = []
    # print(Data.InitialImpurities)
    InitialConcentration = Data.Impurities/(Data.Volume*1E3)
    # print(InitialConcentration)
    for ii, (Time,Diff) in enumerate(zip(Data.Time, Data.DiffConstants)):
        FlowRate.append(GetFlowRate(Time, Diff*DoTimeConversion(TimeScale), Data.Thickness, InitialConcentration[ii], Data.Area))
        # InitialConcentration = InitialConcentration - FlowRate[ii]
    for ii,x in enumerate(FlowRate): 
        # if ii % 100 == 0: 
        #     input()
        print(ii, Data.Time[ii], FlowRate[ii], Data.Impurities[ii], Data.DiffConstants[ii])
    # print(FlowRate)
    # quit()
    FlowRate = np.array(FlowRate) / DoTimeConversion(TimeScale) * Data.Area
    print()
    # test = 4*InitialConcentration*Data.DiffConstants/Data.Thickness * np.exp(- (np.pi/Data.Thickness)**2 * Data.DiffConstants * Data.Time )
    # plt.plot(Data.Time, test)
    # plt.xlabel('Time [Days]')
    # plt.ylabel('Flow Rate [cm$^2$/%s]' % TimeScale)
    # plt.yscale('log')
    # plt.tight_layout()
    # plt.savefig('flow.pdf')
    # plt.show()
    # plt.close()
    # quit()

    # print()
    fig = plt.figure(figsize=(8,6))
    plt.plot(Data.Time, FlowRate)
    plt.xlabel('Time [Days]')
    plt.ylabel('Flow Rate [cm$^2$/%s]' % TimeScale)
    plt.yscale('log')
    plt.tight_layout()
    plt.xlim(0,2)
    plt.ylim(1E-9,1E20)
    # plt.savefig('flow.pdf')

    if Units == '#': pass 
    if Units == 'mBar Liter': 
        FlowRate = FlowRate * (BoltzmannL * Data.Temp) / Avogadro
    
    fig = plt.figure(figsize=(8,6))
    plt.plot(Data.Time, FlowRate)
    plt.xlabel('Time [Days]')
    plt.ylabel('Flow Rate [cm$^2$/%s]' % TimeScale)
    plt.yscale('log')
    plt.tight_layout()
    # plt.savefig('flow_mbar.pdf')
    plt.show()
    plt.close()
    return np.asarray(FlowRate)

def GetFlowRate(Time, Diff, Thickness, Conc, Area): 
    value = 0.0
    for N in range(0,1000,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*Time)
        value = value + factor*(4.0*Conc*Diff)/Thickness
    value = value * Area 
    return value

def GetPartsConversion(Units):
    if Units == 'ppm': return 1E6
    if Units == 'ppb': return 1E9
    if Units == 'ppt': return 1E12

def GetInitialImpurities(Data, Units): 
    ImpurityVolume = Data.Volume * Data.Solubility * Data.Abundance
    ImpurityMass = ImpurityVolume/IdealGasMolarVolume * Data.MolarMass
    if Units == 'Mass':
        pass
    if 'pp' in Units: 
        ImpurityMass = ImpurityMass / Data.XeMass * GetPartsConversion(Units)
    if Units == '#': 
        ImpurityMass = ImpurityMass / Data.MolarMass * Avogadro
    return ImpurityMass

def GetDiffTemp(Data, Temperatures, TimeScale): 
    DiffTemp = []
    for Temp in Temperatures: 
        Diff = Data.Diffusion * np.exp(Data.ActivationEnergy/Boltzmann * ((1.0/293.15) - (1.0/Temp)))
        DiffTemp.append(Diff)
    return np.asarray(DiffTemp) * DoTimeConversion(TimeScale)

    Time = [0, 12, 24, 48, 288, 768]
    Temperature = [293.15, 293.15, 200, 200, 164, 164]
    TimeScale = 'Hour'

    # Time = [0, 0.5, 1, 2, 11, 31]
    # Temperature = [293.15, 293.15, 200, 200, 164, 164]
    # TimeScale = 'Days'

    Time = [0, 1, 20, 100, 200, 300]
    Temperature = [293.15, 293.15, 200, 200, 164, 164]