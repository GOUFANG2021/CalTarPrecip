# Code to model the ion equilibrium in wine and find the supersaturation ratio # of calcium tartrate.
# Written by Jack Muir 2022 - 2024
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import copy

def IonArrays():
    # Extract data from Excel and format. Not set up for every possible ion
    # combination e.g. K2Tar. 
    Names_A = Ion_Data['Anion Names'].dropna()
    Names_A = Names_A.values.tolist()
    z_A = Ion_Data['Anion Charges'].dropna()
    z_A = z_A.values.tolist()
    sigma_A = Size_Data.iloc[:, 2].dropna()
    sigma_A = sigma_A.values.tolist()
    Names_C = Ion_Data['Cation Names'].dropna()
    Names_C = Names_C.values.tolist()
    z_C = Ion_Data['Cation Charges'].dropna()
    z_C = z_C.values.tolist()
    sigma_C = Size_Data.iloc[0, :].dropna()
    sigma_C = sigma_C.values.tolist()
    Mw_Array = Conc_Data['Molar Mass (g/mol)']
    Mw_Array = Mw_Array.values.tolist()
    C_T_Array_set = Conc_Data['Wine Concentration (mol/L)']
    C_T_Array_set = C_T_Array_set.values.tolist()
    C_T_Array_Crystals = [0]*len(C_T_Array_set)
    for cc in range(0, len(Crystal_T_Array_o), 1):
            if np.isnan(Crystal_T_Array_o[cc]) == False:
                C_T_Array_Crystals[cc] = Crystal_T_Array_o[cc]          
    Tot_Names = Conc_Data['Symbol']
    Tot_Names = Tot_Names.values.tolist()
    sigma_ip = Size_Data.iloc[1::, 3::].dropna().T
    sigma_ip = sigma_ip.values.tolist()
    global nC_0
    global nA_0
    nT_0 = len(C_T_Array_o)
    nA_0 = len(Symbol_A)
    nC_0 = len(Symbol_C)
    Names_ip = [[0]*nA_0 for i in range(nC_0)]
    z_ip = [[0]*nA_0 for i in range(nC_0)]
    Names_ip_all = [[0]*nA_0 for i in range(nC_0)]
    global symbol_0
    symbol_0 = []
    for ii in range(0, nT_0, 1):
        if C_T_Array_o[ii] <= 0:
            Mw_Array[ii] = 'remove'
            C_T_Array_set[ii] = 'remove'
            C_T_Array_Crystals[ii] = 'remove'
            Tot_Names[ii] = 'remove'
            symbol_0.append(C_T_symbol[ii])
    for aa in range(0, nA_0, 1):
        for cc in range(0, nC_0, 1):
            if Symbol_A[aa] in symbol_0:
                Names_A[aa] = 'remove'
                z_A[aa] = 'remove'
                sigma_A[aa] = 'remove'
                sigma_ip[cc][aa] = 'remove'
                z_ip[cc][aa] = 'remove'
                Names_ip[cc][aa] = 'remove'
                Names_ip_all[cc][aa] = 'remove'
            elif Symbol_C[cc] in symbol_0:
                Names_C[cc] = 'remove'
                z_C[cc] = 'remove'
                sigma_C[cc] = 'remove'
                sigma_ip[cc][aa] = 'remove'
                z_ip[cc][aa] = 'remove'
                Names_ip[cc][aa] = 'remove'
                Names_ip_all[cc][aa] = 'remove'
            else:
                Names_ip[cc][aa] = Names_C[cc] + Names_A[aa]
                # Change any instances of 'HH' to 'H2' and 'HOH' to 'H2O' etc.
                # Should be added to if more ion pair combinations were used.
                Names_ip[cc][aa] = Names_ip[cc][aa].replace("HH2", "H3")
                Names_ip[cc][aa] = Names_ip[cc][aa].replace("HH", "H2")
                Names_ip[cc][aa] = Names_ip[cc][aa].replace("HOH", "H2O")
                Names_ip_all[cc][aa] = Names_ip[cc][aa]
                z_ip[cc][aa] = z_C[cc] + z_A[aa]
                if Names_ip[cc][aa] in Names_A:
                    Names_ip[cc][aa] = 'duplicate'
                    z_ip[cc][aa] = 'duplicate'
                    sigma_ip[cc][aa] = 'duplicate'
                if Names_ip[cc][aa] == 'H2O':
                    z_ip[cc][aa] = 'H2O'
                    sigma_ip[cc][aa] = 'H2O'
    z_ip = [i for sublist in z_ip for i in sublist]
    Names_ip = [i for sublist in Names_ip for i in sublist]
    sigma_ip = [i for sublist in sigma_ip for i in sublist]
    # Remove any compounds with a concentration of 0 mol/L.
    remove = 'remove'
    Names_A = [i for i in Names_A if i != remove]
    z_A = [i for i in z_A if i != remove]
    sigma_A = [i for i in sigma_A if i != remove]
    Names_C = [i for i in Names_C if i != remove]
    z_C = [i for i in z_C if i != remove]
    sigma_C = [i for i in sigma_C if i != remove]
    Mw_Array = [i for i in Mw_Array if i != remove]
    C_T_Array_set = [i for i in C_T_Array_set if i != remove]
    C_T_Array_Crystals = [i for i in C_T_Array_Crystals if i != remove]
    Tot_Names = [i for i in Tot_Names if i != remove]
    sigma_ip = [i for i in sigma_ip if i != remove]
    z_ip = [i for i in z_ip if i != remove]
    Names_ip = [i for i in Names_ip if i != remove]
    for ii in range(0, nC_0, 1):
        Names_ip_all[ii] = [i for i in Names_ip_all[ii] if i != remove]
    Names_ip_all = [x for x in Names_ip_all if x != []]
    n_A = len(Names_A)
    n_C = len(Names_C)
    # Remove any duplicates.
    duplicate = "duplicate"
    Names_ip = [i for i in Names_ip if i != duplicate]
    z_ip = [i for i in z_ip if i != duplicate]
    sigma_ip = [i for i in sigma_ip if i != duplicate]
    # Remove H2O
    Names_ip.remove('H2O')
    z_ip.remove('H2O')
    sigma_ip.remove('H2O')
    z = z_A + z_C + z_ip
    n_T = len(z)
    Names_all = Names_A + Names_C + Names_ip
    sigma = sigma_A + sigma_C + sigma_ip
    Results = [Names_A, n_A, Names_C, n_C, Names_ip, Names_ip_all, z, n_T,
               Names_all, sigma, Mw_Array, C_T_Array_set, Tot_Names, 
               C_T_Array_Crystals]
    return Results

def InitialGuesses(C_T_Array):
    # Gives initial concentration guesses (µmol/L) for Newton's method based
    # on total concentrations. 
    x0 = np.zeros([n_T, 1])
    for ii in range(0, len(Tot_Names), 1):
        for jj in range(0, n_T, 1):
            if (Tot_Names[ii] in Names_all[jj]) and (x0[jj] == 0):
                x0[jj] = C_T_Array[ii]
            #  Make initial guess the limiting concentration. 
            elif (Tot_Names[ii] in Names_all[jj]):
                x0[jj] = min(x0[jj], C_T_Array[ii])
            # H and OH are assumed to be plentiful. 
            elif Names_all[jj] == 'H':
                x0[jj] = 0.1 * 1e6
            elif Names_all[jj] == 'OH':
                x0[jj] = 0.1 * 1e6
    return x0

def Bounds(C_T_Array):
    # Gives upper bounds (µmol/L) for Newton's method based on total
    # concentrations. 
    ub_ions = np.zeros([n_T, 1])
    for ii in range(0, len(Tot_Names), 1):
        for jj in range(0, n_T, 1):
            if (Tot_Names[ii] in Names_all[jj]) and (ub_ions[jj] == 0):
                ub_ions[jj] = C_T_Array[ii]
            # Use the limiting concentration to set the bound.
            elif (Tot_Names[ii] in Names_all[jj]):
                ub_ions[jj] = min(ub_ions[jj], C_T_Array[ii])
            # H and OH are assumed to be plentiful.
            elif Names_all[jj] == 'H':
                ub_ions[jj] = 1e20  
            elif Names_all[jj] == 'OH':
                ub_ions[jj] = 1e20  
    # Bounds for the precipitates based on limiting concentration.  
    ub_ppt = np.zeros([n_ppt, 1])
    for ii in range(0, len(Tot_Names), 1):
        for jj in range(0, n_ppt, 1):
            if (Tot_Names[ii] in ppt_names[jj]) and (ub_ppt[jj] == 0):
                ub_ppt[jj] = C_T_Array[ii]
            elif (Tot_Names[ii] in ppt_names[jj]):
                ub_ppt[jj] = max(ub_ppt[jj], C_T_Array[ii])
    ub = np.append(ub_ions, [ub_ppt])
    ub = np.array([ub]).reshape(len(ub), 1)
    return ub

def WaterConc(C_EtOH, C_T_solutes):
    # Find concentration of water using a simple Newton’s method. Set up for
    # 25 °C, ethanol mass fractions of 0 - 0.63 (0 - 12 mol/L), and solute
    # mass fractions of 0 - 0.2 (0 - 3.8 mol/L).   
    C_EtOH = C_EtOH/1e6 #(mol/L)
    C_T_solutes = C_T_solutes/1e6 #(mol/L)
    # Molar mass (g/mol)
    Mw_H2O = 18.016
    Mw_EtOH = 46.07
    # Density (g/L)
    rho_EtOH = C_EtOH*Mw_EtOH
    rho_solutes = np.sum(C_T_solutes*Mw_solutes*1e6)
    C_H2O_Guess = 55.5  #(mol/L)
    # Function f(C_H2O) = 0
    # Fitted to empirical data taken from Perry and Green (1997) and
    # Galleguillos et al. (2003)
    def f(C_H2O): return (37.5*np.sum(C_T_solutes) +
                          (54.17*((C_EtOH*Mw_EtOH)/(C_H2O*Mw_H2O))**6
                           - 296.56*((C_EtOH*Mw_EtOH)/(C_H2O*Mw_H2O))**5
                           + 623.43*((C_EtOH*Mw_EtOH)/(C_H2O*Mw_H2O))**4
                           - 632.20*((C_EtOH*Mw_EtOH)/(C_H2O*Mw_H2O))**3
                           + 342.29*((C_EtOH*Mw_EtOH)/(C_H2O*Mw_H2O))**2
                           - 178.62*((C_EtOH*Mw_EtOH)/(C_H2O*Mw_H2O)) + 
                           997.08) - rho_EtOH - rho_solutes)/Mw_H2O - C_H2O
    # Derivative of function f
    def d(C_H2O): return (-6*54.17 * (C_EtOH*Mw_EtOH/Mw_H2O)**6 * C_H2O**-7
                          + 5*296.56 * (C_EtOH*Mw_EtOH/Mw_H2O)**5 * C_H2O**-6
                          - 4*623.43 * (C_EtOH*Mw_EtOH/Mw_H2O)**4 * C_H2O**-5
                          + 3*632.20 * (C_EtOH*Mw_EtOH/Mw_H2O)**3 * C_H2O**-4
                          - 2*342.29 * (C_EtOH*Mw_EtOH/Mw_H2O)**2 * C_H2O**-3
                          + 178.62 * (C_EtOH*Mw_EtOH/Mw_H2O)
                          * C_H2O**-2)/Mw_H2O - 1
    C_H2O = SimpleNewtons(C_H2O_Guess, f, d)  #(mol/L)
    C_H2O = C_H2O*1e6  #(µmol/L)
    return C_H2O

def SimpleNewtons(x, f, d):
    # Simple Newton's method
    # Takes initial guess x, function f, and derivative d.
    x_i = x     # Initial Guess
    N = 50      # Max iterations
    tol = 1e-6  # Tolerance
    approxArray = np.zeros(N+1)
    approxArray[0] = x_i
    errorArray = np.zeros(N+1)
    errorArray[0] = abs(f(x_i))
    changeInApproxArray = np.zeros(N)
    for nn in range(0, N, 1):
        x = x - f(x)/d(x)
        approxArray[nn + 1] = x
        errorArray[nn + 1] = abs(f(x))
        changeInApproxArray[nn] = abs(approxArray[nn + 1] - approxArray[nn])
        if (changeInApproxArray[nn] < tol) and (errorArray[nn+1] < tol):
            break
    return x

def AssocConstants(C_EtOH):
    # Get association constants at 20 - 25 °C (L/mol)
    Ka_EtOH_0 = Ka_Data.iloc[:, 2::].values.tolist()
    C_EtOH = C_EtOH/1e6  #(mol/L)
    Sym_H = Symbol_C.index('H')  
    # Estimate effect of ethanol for species with no specific data.  
    Ka = 10**(0.0938*C_EtOH + np.log10(Ka_EtOH_0))
    # Estimate effect of ethanol for species with specific data from Usseglio-
    # Tomasset & Bosia (1978)   
    for ii in range(0, len(Symbol_A_ion), 1):
        if Symbol_A_ion[ii] == 'H2Cit-':
            Ka[ii, Sym_H] = 68.98*C_EtOH**2 + 174.33*C_EtOH + 1431.7
        elif Symbol_A_ion[ii] == 'HCit2-':
            Ka[ii, Sym_H] = 2272.9*C_EtOH**2 + 10164*C_EtOH + 53599
        elif Symbol_A_ion[ii] == 'HTar-':
            Ka[ii, Sym_H] = 51.604*C_EtOH**2 + 139.44*C_EtOH + 1197.5
        elif Symbol_A_ion[ii] == 'Tar2-':
            Ka[ii, Sym_H] = 1551.8*C_EtOH**2 + 3626.2*C_EtOH + 24663
        elif Symbol_A_ion[ii] == 'HMal-':
            Ka[ii, Sym_H] = 132.76*C_EtOH**2 + 388.47*C_EtOH + 2994
        elif Symbol_A_ion[ii] == 'Mal2-':
            Ka[ii, Sym_H] = 8026.6*C_EtOH**2 + 23590*C_EtOH + 126724
        elif Symbol_A_ion[ii] == 'HSuc-':
            Ka[ii, Sym_H] = 856.81*C_EtOH**2 + 2019.4*C_EtOH + 16208
        elif Symbol_A_ion[ii] == 'Suc2-':
            Ka[ii, Sym_H] = 35048*C_EtOH**2 + 50694*C_EtOH + 430860
        elif Symbol_A_ion[ii] == 'Lac-':
            Ka[ii, Sym_H] = 358.61*C_EtOH**2 + 1000.3*C_EtOH + 7807.3
        elif Symbol_A_ion[ii] == 'Ace-':
            Ka[ii, Sym_H] = 2394.6*C_EtOH**2 + 4492*C_EtOH + 57265
    Ka = Ka.tolist()
    # Remove any that have total concentrations of 0.
    for aa in range(0, nA_0, 1):
        for cc in range(0, nC_0, 1):
            if Symbol_A[aa] in symbol_0:
                Ka[aa][cc] = 'remove'
            elif Symbol_C[cc] in symbol_0:
                Ka[aa][cc] = 'remove'
    Ka = [i for sublist in Ka for i in sublist]
    remove = "remove"
    Ka = [i for i in Ka if i != remove]
    Ka = np.array(Ka).reshape(n_A, n_C)
    return Ka

def H2O_EtOH_Dielectric(C_H2O, C_EtOH):
    # Gets dielectric constant for binary ethanol and water mixture at 25 °C.
    # Uses data from Usseglio-Tomasset and Bosia (1978) for 0 – 30 wt% EtOH.
    C_m_H2O = C_H2O * Mw_H2O
    C_m_EtOH = C_EtOH * Mw_EtOH
    mass_frac_EtOH = (C_m_EtOH)/(C_m_EtOH + C_m_H2O)
    D = -3.9474*mass_frac_EtOH**2 - 56.942*mass_frac_EtOH + 78.54
    return D

def Dielectric_Constant(C_H2O, C_EtOH, C_solutes):
    # Gets solution dielectric constant based on Zuber et al. (2014).
    # Critical properties from Perry and Green (1997).
    vc_H2O = 56        # Critical volume of pure water (cm^3/mol)
    Zc_H2O = 0.228     # Critical compressibility factor of pure water (-)
    Tc_H2O = 647.13    # Critical temperature of pure water (K)
    vc_EtOH = 168      # Critical volume of pure ethanol (cm^3/mol)
    Zc_EtOH = 0.240    # Critical compressibility factor of pure ethanol (-)
    Tc_EtOH = 513.92   # Critical temperature of pure ethanol (K)
    # Parameters from Zuber et al. (2014)
    ac_H2O = 2.6    
    aa_H2O = 7.89   
    ac_EtOH = 9.04  
    aa_EtOH = 21.36 
    C_T = C_H2O + C_EtOH + np.sum(C_solutes) #(µmol/L)
    x_H2O = C_H2O/C_T     
    x_EtOH = C_EtOH/C_T   
    x_solutes = C_solutes/C_T 
    Tr_H2O = T_K/Tc_H2O
    Tr_EtOH = T_K/Tc_EtOH
    # Molar volume (cm^2/mol) from Rackett equation (Vetere, 1992).
    v_H2O = vc_H2O * Zc_H2O**((1-Tr_H2O)**(2/7))
    v_EtOH = vc_EtOH * Zc_EtOH**((1-Tr_EtOH)**(2/7))
    # Volume fraction
    phi_H2O = (x_H2O*v_H2O) / (x_H2O*v_H2O + x_EtOH*v_EtOH)
    phi_EtOH = (x_EtOH*v_EtOH) / (x_H2O*v_H2O + x_EtOH*v_EtOH)
    epsilon_solvents = H2O_EtOH_Dielectric(C_H2O, C_EtOH)
    AnionSum = 0
    CationSum = 0
    for ii in range(0, n_T, 1):
        if z[ii] < 0:
            sum_ii = x_solutes[ii]*(aa_H2O*phi_H2O + aa_EtOH*phi_EtOH)
            AnionSum = AnionSum + sum_ii
        elif z[ii] > 0:
            sum_ii = x_solutes[ii]*(ac_H2O*phi_H2O + ac_EtOH*phi_EtOH)
            CationSum = CationSum + sum_ii
    epsilon_r = epsilon_solvents/(1 + AnionSum + CationSum)
    return epsilon_r

def ActCoeff(z, sigma, rho, C_H2O, C_EtOH, C_solutes):
    # Get activity coefficients using MSA method
    # Electrostatic contribution
    e_0 = 1.6e-19          # Electron charge (Coulomb)
    KB = 1.38e-23          # Boltzmann constant (J/K)
    epsilon_0 = 8.85e-12   # Permittivity of a vacuum (Farads/m)
    epsilon_r = Dielectric_Constant(C_H2O, C_EtOH, C_solutes)
    L_B = (e_0**2) / (4*math.pi*epsilon_0*epsilon_r*KB*T_K)
    kappa = np.sqrt(4*math.pi*L_B*sum(rho*z**2))
    Gamma = np.sqrt(kappa**2/4)
    Delta = 1 - (math.pi/6*sum(sum(rho*sigma**3)))
    Omega = 1 + ((math.pi/2*Delta)*sum((rho*sigma**3)/(1+Gamma*sigma)))
    xi = (math.pi/(2*Omega*Delta))*sum((rho*sigma*z)/(1+Gamma*sigma))
    part1 = (Gamma*z**2/(1 + Gamma*sigma))
    part2 = (2*z - xi*sigma**2)/(1 + Gamma*sigma) + (xi*sigma**2)/3
    part3 = (xi*sigma)*part2**2
    lny_es = -L_B * (part1 + part3)
    # Hard sphere contribution
    X0 = (math.pi/6)*sum(rho*sigma**0)
    X1 = (math.pi/6)*sum(rho*sigma**1)
    X2 = (math.pi/6)*sum(rho*sigma**2)
    X3 = (math.pi/6)*sum(rho*sigma**3)
    F1 = (3*X2)/(1-X3)
    F2 = (3*X1/(1-X3))+(3*X2**2/(X3*(1-X3)**2))+(3*X2**2*np.log(1-X3)/X3**2)
    F3 = (X0-X2**3/X3**2)*(1/(1-X3))+((3*X1*X2-(X2**3/X3**2))/(1-X3)**2)\
    + (2*X2**3/(X3*(1-X3)**3))+(3*X2**2*np.log(1-X3)/X3**3)
    lny_hs = -np.log(1-X3) + sigma*F1 + sigma**2*F2 + sigma**3*F3
    # Electrostatic contribution + hard sphere contribution
    y = np.exp(lny_es + lny_hs)
    return y

def CaTar_Ksp(A, C_EtOH, pH):
    # Estimate solubility product using data from Curvelo-Garcia (1987)
    # for CaTar at 18 °C for ethanol contents of 10 - 20% and pH
    # ranging from 2.5 - 4.5.
    A_H = (A[H_Index])/1e6  #(mol/L)
    pH = -np.log10(A_H)  
    Mw_EtOH = 46.07 #(g/mol)
    # Density from Perry and Green (1997)
    Density_EtOH_25 = 785.06  # g_EtOH/L_EtOH at 25 °C
    Density_EtOH_20 = 789.34  # g_EtOH/L_EtOH at 20 °C
    if T_C == 25:
        Density = Density_EtOH_25
    elif T_C == 20:
        Density = Density_EtOH_20
    else:
        print('Error: Temperature not in range')
    C_EtOH = C_EtOH/1e6  # mol/L
    EtOH_g = C_EtOH*Mw_EtOH  # g_EtOH/L_Solution
    EtOH_ABV = EtOH_g/Density * 100  # L_EtOH/L_Solution
    Gradient = 0.0027*EtOH_ABV**2 + 0.0045*EtOH_ABV - 0.1348
    Intercept_pH = 2.28
    Intercept_Ksp = 6.76
    pKsp = Gradient * (pH - Intercept_pH) + Intercept_Ksp
    Ksp = 10**(-pKsp) #(mol/L)^2
    return Ksp

def Equations(x, param):
    C_T_Array = param[0]   #(µmol/L)
    C_H2O = param[1]       #(µmol/L)
    C_EtOH = param[2]      #(µmol/L)
    C_Sugar = param[3]     #(µmol/L)
    Flag = param[4]        
    if Flag == False:
        f = np.zeros([n_T+n_ppt, 1])
    elif Flag == True:
        C_T_Array[S1_T_Index] = x[n_T+n_ppt]
        f = np.zeros([n_T+n_ppt+1, 1])
    K = AssocConstants(C_EtOH)/1e6 #(L/µmol)
    Conc_Array = np.append(x[0:n_T], [C_EtOH, C_Sugar])   #(µmol/L)
    Conc_Array = Conc_Array.reshape(len(Conc_Array), 1)   #(µmol/L)
    rho = Conc_Array/1e6*N_A*1000                         #(number/m3)
    C_solutes = np.append(x[0:n_T], [C_Sugar])            #(µmol/L)
    C_solutes = C_solutes.reshape(len(C_solutes), 1)      #(µmol/L)
    y = ActCoeff(z, sigma, rho, C_H2O, C_EtOH, C_solutes)
    # Activities (µmol/L)
    A = y[0:n_T]*x[0:n_T] 
    A_A = A[0:n_A]      
    A_C = A[n_A:n_A+n_C]
    A_ip = A[n_A+n_C:n_T]
    y_H2O = 1 # Assumed
    A_H2O = y_H2O * C_H2O
    A_ip = np.append(A_ip, A_H2O)
    Current_pH = -np.log10((A[H_Index])/1e6)
    n = 0
    a = 0
    # Solubility products (µmol/L)^2
    Ksp = Precipitate_Data['Solubility Product Constant (mol/L)2']*(1e6)**2 
    Ksp = np.array([Ksp]).T
    # Solubility product for CaTar (µmol/L)^2. 
    for ii in range(0, n_ppt, 1):
        if ppt_names[ii] == 'CaTar':
            Ksp[ii] = CaTar_Ksp(A, C_EtOH, Current_pH)*(1e6)**2
    v = Precipitate_Data['v']
    v = np.array([v]).T
    Supersat_ratio = x[n_T:n_T+n_ppt]
    # Equilibrium equations: f = A_ionpair - Ka(A_cation * A_anion) = 0
    for cc in range(0, n_C, 1):
        for aa in range(0, n_A, 1):
            if Names_ip_all[cc][aa] in Names_A:
                for ii in range(0, n_A, 1):
                    if Names_A[ii] == Names_ip_all[cc][aa]:
                        f[n] = (A_A[ii] - K[aa, cc]*A_C[cc]*A_A[aa])
                        n = n + 1    
            elif Names_ip_all[cc][aa] == 'H2O':
                logKw = -4787.3/T_K - 7.1332*np.log10(T_K) - 0.010365*T_K\
                + 22.801
                Kw = 10**logKw #(L/mol)
                Kw_intrinsic = ((A_H2O/1e6)/Kw) #(L/mol)
                Kw_intrinsic = (10**(0.0938*C_EtOH/1e6\
                                + np.log10(Kw_intrinsic)))/1e6 #(L/µmol)
                f[n] = (A_H2O - Kw_intrinsic*A_C[cc]*A_A[aa])
                n = n + 1    
            else:
                f[n] = (A_ip[a] - K[aa, cc]*A_C[cc]*A_A[aa])
                n = n + 1    
                a = a + 1    
    # Electroneutrality: sum(z_i*C_i) = 0
    Scale = sum(x[0:n_T])
    f[n_A*n_C] = sum(z[0:n_T]*x[0:n_T])/Scale
    # Component balances: f = total conc - conc in all different forms = 0 
    for ii in range(0, len(C_T_Array), 1):
        All_Forms = 0
        for jj in range(0, n_T, 1):
            if Tot_Names[ii] in Names_all[jj]:
                All_Forms = All_Forms + x[jj]
        Scale = C_T_Array[ii]
        f[n_A*n_C + ii+1] = (C_T_Array[ii] - All_Forms)/Scale
    # Solubility equations: f = (IAP/Ksp)^(1/v) - Supersaturation ratio = 0
    for ii in range(0, n_ppt, 1):
        ppt_an_index = Names_all.index(ppt_anion_names[ii])
        ppt_cat_index = Names_all.index(ppt_cation_names[ii])
        f[n_T + ii] = (A[ppt_cat_index]*A[ppt_an_index]/Ksp[ii])**(1/v[ii])\
        - Supersat_ratio[ii]
    # pH adjustment: f = pH - Measured pH = -log10(A_H+) - Measured pH = 0
    if Flag == True:
        f[n_T+n_ppt] = (-np.log10(A[H_Index]/1e6) - Measured_pH)
    return f

def Jacobian(x, param):
    # Gives Jacobian matrix.
    x_length = len(x)
    J = np.zeros((x_length, x_length))
    fx = Equations(x, param)
    for ii in range(0, x_length, 1):
        x_copy = copy.deepcopy(x[ii])
        # Perturb x_ii with a small change
        dx_ii = 1e-7 * x[ii]
        # If dx_ii is less than 1e-9, then set it to 1e-9
        if abs(dx_ii) < 1e-9:
            dx_ii = 1e-9
        x[ii] = x[ii] + dx_ii
        # Evaluate f(x) with the perturbed value
        fx_ii = Equations(x, param)
        Jac = (fx_ii - fx)/dx_ii
        J[:, ii] = Jac[:, 0]
        # Restore x_ii to its original value.
        x[ii] = x_copy
    return J, fx

def NewtonsMethod(x0, options, param):
    # Newton's method. 
    maxits = options[0]   
    xTol = options[1]     
    fTol = options[2]     
    lb = options[3]       
    ub = options[4]       
    run = options[5]   
    Wine = options[6]
    x = x0  
    for ii in range(0, maxits, 1):
        [J, fx] = Jacobian(x, param)
        fx = np.reshape(fx, [len(fx),1])
        if Scaling == True:
            # Modifications for automatic scaling
            D1diag = 1/np.sqrt(np.max(abs(J), axis=1))
            D2diag = 1/np.sqrt(np.max(abs(J), axis=0))
            D1 = np.diag(D1diag)
            D2 = np.diag(D2diag)
            D1JD2 = np.matmul(D1, np.matmul(J, D2))
            D1fx = D1diag*fx
            D2delta_x = np.linalg.solve(D1JD2, D1fx)
            delta_x = D2diag*D2delta_x  
        else:
            delta_x = np.linalg.solve(J, fx)
        x = x - 1 * delta_x
        # Check norm(h) and norm(fx) for convergence
        if (np.linalg.norm(delta_x, np.inf) < xTol) and (np.linalg.norm(fx,
            np.inf) < fTol):
            # print('Newton\'s method converged for ' + Wine +
            #       ', run number = ' + run)
            status = True
            break
        elif ii == (maxits - 1):
            # print('Newton\'s method did not converge for ' + Wine +
            #       ', run number = ' + run)
            status = False
        # Check upper and lower bounds
        for i in range(0, len(x), 1):
            if x[i] < lb[i]:
                x[i] = lb[i]
                # print('Lower bound hit')
            if x[i] > ub[i]:
                x[i] = ub[i]
                # print('Upper bound hit')
    return x, status

def PrecipitationRate(Supersat_ratio, Precipitates):
    # Gives the precipitation rate.
    # Convert rate constants to units of min^-1.
    k_prec = Precipitate_Data['Precipitation Rate Constant (hr-1)']/60
    k_prec = np.array([k_prec]).T
    n_prec = Precipitate_Data['Order of the Precipitation Reaction']
    n_prec = np.array([n_prec]).T
    k_diss = Precipitate_Data['Dissolution Rate Constant (hr-1)']/60 
    k_diss = np.array([k_diss]).T
    n_diss = Precipitate_Data['Order of the Dissolution Reaction']  
    n_diss = np.array([n_diss]).T
    sigma = Supersat_ratio - 1 # Relative supersaturation
    r_prec = np.zeros([n_ppt, 1])
    r_diss = np.zeros([n_ppt, 1])
    for jj in range(0, n_ppt, 1):
        if sigma[jj] > 0:
            # Solution supersaturated so precipitation will occur.
            # No dissolution.
            r_prec[jj] = k_prec[jj] * Precipitates[jj] * sigma[jj]**n_prec[jj]
            r_diss[jj] = 0
        elif Precipitates[jj] > 0:
            # Solution not supersaturated so no precipitation.
            # Precipitate already exists, so dissolution will occur.
            r_diss[jj] = k_diss[jj] * Precipitates[jj] * sigma[jj]**n_diss[jj]
            r_prec[jj] = 0
        else:
            # Solution is not supersaturated so no precipitation.
            # No precipitate already exists, so no dissolution will occur.
            r_prec[jj] = 0
            r_diss[jj] = 0
    # Total rate of precipitation = precipitation + dissolution ((µmol/L)/min)
    dCdt = r_prec + r_diss
    return dCdt

def func(t, R_P, ii):
    # Finds the current values at time t. 
    for nn in range(0, len(C_T_Array), 1):
        All_Forms = 0
        for jj in range(0, n_ppt, 1):
            if Tot_Names[nn] in ppt_names[jj]:
                All_Forms = All_Forms + R_P[jj]
        # Subtract precipitated species from the solution
        C_T_Array[nn] = C_T_Array_set[nn] - All_Forms       #(µmol/L)
    C_T_solutes = np.append(C_T_Array, [C_Sugar])           #(µmol/L)
    C_T_solutes = C_T_solutes.reshape(len(C_T_solutes), 1)  #(µmol/L)
    C_H2O = WaterConc(C_EtOH, C_T_solutes)                  #(µmol/L)
    x0 = InitialGuesses(C_T_Array) #(µmol/L)
    Supersat_ratio_guess = np.ones([n_ppt, 1])
    x0 = np.append(x0, Supersat_ratio_guess)
    x0 = np.array([x0]).reshape(len(x0), 1)
    if ii > 0:
        x0 = R_C[ii-1, :]
        x0 = np.reshape(x0, [n_T+n_ppt, 1])
    param = [C_T_Array, C_H2O, C_EtOH, C_Sugar, False]
    lb = np.ones([n_T+n_ppt, 1])*1e-20  # Lower bounds (µmol/L)
    ub = Bounds(C_T_Array)              # Upper bounds (µmol/L)
    options = [100, 1e-9, 1e-9, lb, ub, '1', Wine[0]]
    x, status = NewtonsMethod(x0, options, param)
    Conc_Array = np.append(x[0:n_T], [C_EtOH, C_Sugar])  #(µmol/L)
    Conc_Array = Conc_Array.reshape(len(Conc_Array), 1)  #(µmol/L)
    rho = Conc_Array/1e6*N_A*1000                        #(number/m3)
    C_solutes = np.append(x[0:n_T], [C_Sugar])           #(µmol/L)
    C_solutes = C_solutes.reshape(len(C_solutes), 1)     #(µmol/L)
    y = ActCoeff(z, sigma, rho, C_H2O, C_EtOH, C_solutes)
    R_y[ii, :] = y[0:n_T, 0]
    R_C[ii, :] = x[0::, 0] #(µmol/L)
    R_CT[ii+1, :] = C_T_Array[0:n_Tot, 0] #(µmol/L)
    Precipitates = R_P #(µmol/L)
    Supersat_ratio = x[n_T: n_T+n_ppt]
    dCdt = PrecipitationRate(Supersat_ratio, Precipitates)
    return dCdt

def ODE_Solver(func, t_values, y0, ODE_options):
    # Solve the ODE with Euler's or Heun's method.
    ODE_Method = ODE_options[0]
    t_step = t_values[2]
    n_steps = t_values[3]
    t = t_values[4]
    y = np.zeros([n_steps+1, n_ppt])
    y[0, :] = y0
    if ODE_Method == 'Eulers':
        for ii in range(0, n_steps, 1):
            dydt1 = func(t[ii], y[ii, :], ii)
            y[ii+1, :] = y[ii, :] + t_step*dydt1.T
            print('Completed ODE Step for t = ' + str(t[ii+1]) + ' min')
    elif ODE_Method == 'Heuns':
        for ii in range(0, n_steps, 1):
            dydt1 = func(t[ii], y[ii, :], ii)
            y[ii+1, :] = y[ii, :] + t_step*dydt1.T
            dydt2 = func(t[ii+1], y[ii+1, :], ii)
            y[ii+1, :] = y[ii, :] + t_step/2*(dydt1.T + dydt2.T)
            print('Completed ODE Step for t = ' + str(t[ii+1]) + ' min')
    else:
        print('ODE_Method not recognised')
    return y

def EtOH_Conversion(EtOH_ABV, T_C):
    # Takes ethanol vol% and converts it to mol/L for temperature of 20 – 25 
    Mw_EtOH = 46.07 #(g/mol)
    # Density of pure ethanol from Perry and Green (1997)
    Density_EtOH_25 = 785.06  # g_EtOH/L_EtOH at 25 °C
    Density_EtOH_20 = 789.34  # g_EtOH/L_EtOH at 20 °C
    if T_C == 25:
        Density = Density_EtOH_25
    elif T_C == 20:
        Density = Density_EtOH_20
    else:
        print('Error: Temperature not in range')
    # EtOH_ABV = [Volume pure ethanol]/[Volume Solution] * 100
    EtOH_g = EtOH_ABV/100 * Density  # g_EtOH/L_Solution Volume
    C_EtOH = EtOH_g/Mw_EtOH  # mol_EtOH/L_Solution
    return C_EtOH
##############################################################################
                              ### OPTIONS ###
# Code description
Description = 'Find the supersaturation ratio for wines W1 - W24.'
# Excel output file name
Output_Name = 'Wine Results.xlsx'
# Set to True to produce an Excel File with the results.
Excel_f = False
# Wine number(s) 
# W_A = np.array([1,5,6,7,8,9,24,2,3,10,13,15,16,17,18,19,21,22,23,12,20])
# Use this array to just do one wine at a time.
W_A = np.array([1])
# If True then sugar is removed from the system.
Sugar_Change = True
# Set to True to add in the concentrations from crystals into the system or
# False to exclude them.
Crystals = False
# Set to True to turn on automatic scaling.
Scaling = False
# Can set the model to 'Reduced' to change the Add_1 array, which will 
# allow it to solve a bit faster. It will still solve when this isn't used. 
Model = 'Reduced'
##############################################################################
T_C = 25                # Temperature (°C)
T_K = T_C + 273.15      # Temperature (K)
N_A = 6.022e23          # Avogadro's number (mol^-1)
z_EtOH = 0
z_Sugar = 0
sigma_non_dissociating = [0.45, 0.88]  #(nm)
# Molar Masses (g/µmol)
Mw_H2O = 18.016/1e6
Mw_EtOH = 46.07/1e6
Mw_Sugar = 180.16/1e6 
Wine_num = len(W_A)
Wine = [0]*Wine_num
W_y = [0]*Wine_num
W_C = [0]*Wine_num
W_pH = [0]*Wine_num
W_pKsp = [0]*Wine_num
W_SS = [0]*Wine_num
W_C_T_Array = [0]*Wine_num
W_Convergence = [0]*Wine_num
# Sequential search array for NaOH or HCl (µmol/L) added to the 'no pH
# adjustment' run. 
if Model == 'Reduced': 
    Add_1 = np.array([10000,5000,15000,20000,0,25000,30000,
                      35000,40000,45000,50000])
else: 
    Add_1 = np.array([0,10000,5000,15000,20000,25000,30000,
                      35000,40000,45000,50000])
# Sequential search array for initial NaOH or HCl guesses (µmol/L)
Add_2 = np.array([20000, 30000,40000,10000,50000,60000,
                  70000,80000,90000,100000,0])
if Crystals == True:
    L = 4
elif Crystals == False:
    L = 2
for ii in range(0, Wine_num, 1): 
    Wine[ii] = 'W' + str(W_A[ii])
    Filename = 'Wine Data.xlsx'
    Excel_data = pd.read_excel(Filename, sheet_name=['Ion_Data', 'Ka_Data',
                               'Size_Data', 'Conc_Data', 
                               'Precipitate_Data', 'Other_Data'])
    Ion_Data = Excel_data.get('Ion_Data')
    Ka_Data = Excel_data.get('Ka_Data')
    Size_Data = Excel_data.get('Size_Data')
    Conc_Data = Excel_data.get('Conc_Data')
    Precipitate_Data = Excel_data.get('Precipitate_Data')
    Other_Data = Excel_data.get('Other_Data')
    C_T_Array_o = Conc_Data['Wine Concentration (mol/L)'].values.tolist()
    Crystal_T_Array_o = Conc_Data['Crystal Concentration (mol/L)']\
    .values.tolist()
    C_T_symbol = Conc_Data['Symbol'].values.tolist()
    Symbol_A = Ion_Data['Anion Symbol'].dropna()
    Symbol_A = Symbol_A.values.tolist()
    Symbol_C = Ion_Data['Cation Symbol'].dropna()
    Symbol_C = Symbol_C.values.tolist()
    Symbol_A_ion = Ka_Data.iloc[:, 1]
    Symbol_A_ion = Symbol_A_ion.values.tolist()
    Ion_Arrays = IonArrays()
    Names_A = Ion_Arrays[0]     
    n_A = Ion_Arrays[1]         
    Names_C = Ion_Arrays[2]     
    n_C = Ion_Arrays[3]        
    Names_ip = Ion_Arrays[4]    
    Names_ip_all = Ion_Arrays[5]
    z = Ion_Arrays[6]
    n_T = Ion_Arrays[7]
    Names_all = Ion_Arrays[8]
    sigma = Ion_Arrays[9]
    Mw_Array = Ion_Arrays[10]      #(g/mol)
    C_T_Array_set = Ion_Arrays[11] #(mol/L)
    Tot_Names = Ion_Arrays[12]
    C_T_Array_Crystals = Ion_Arrays[13]     
    n_ppt = len(Precipitate_Data)
    H_Index = Names_all.index('H')
    z = np.append(z, [z_EtOH, z_Sugar])
    z = np.array([z]).reshape(len(z), 1)
    sigma = 1e-9*(np.array([sigma + sigma_non_dissociating])).T  #(m)
    # Molar masses (g/µmol)
    Mw_Array = np.array([Mw_Array]).reshape(len(Mw_Array), 1)/1e6   
    Mw_solutes = np.append(Mw_Array, [Mw_Sugar])                    
    Mw_solutes = np.array([Mw_solutes]).reshape(len(Mw_solutes), 1) 
    # Total concentrations (µmol/L)
    C_T_Array_set = np.array(C_T_Array_set)*1e6 # Not changed by each run 
    C_T_Array = np.array([C_T_Array_set]).T     # Changed by each run
    ppt_anion_names = Precipitate_Data['Anion Names']
    ppt_cation_names = Precipitate_Data['Cation Names']
    ppt_names = ppt_cation_names + ppt_anion_names
    EtOH_ABV = Other_Data['Ethanol Content (vol%)']
    EtOH_ABV = float(EtOH_ABV[0])
    C_EtOH = EtOH_Conversion(EtOH_ABV, T_C)*1e6 #(µmol/L)
    C_Sugar = Other_Data['Sugar Content (g/L)']
    C_Sugar = float(C_Sugar[0])
    C_Sugar = C_Sugar/Mw_Sugar #(µmol/L)
    if Sugar_Change == True:
        C_Sugar = 0
    Measured_pH = Other_Data['pH']
    Measured_pH = round(float(Measured_pH[0]), 2)
    R_y = np.zeros([L, n_T])      
    R_C = np.zeros([L, n_T])      
    R_SS = np.zeros([L, n_ppt])   
    R_pKsp = np.zeros([L, n_ppt])  
    R_pH = np.zeros([L, 1])       
    R_C_T_Array = np.zeros([L, len(C_T_Array)]) 
    R_Convergence = ['Did not converge']*L 
    status2 = False
    for jj in range(0, L, 1):
        status = False
        if jj == 1:
            Flag = True
        else:
            Flag = False
        C_T_Array[:, 0] = C_T_Array_set          
        # Add the crystals to the system on certain runs of Newton's method.
        if jj > 1:
            C_T_Array = C_T_Array + np.array([C_T_Array_Crystals]).T*1e6
        if Flag == True:
            # Compare the model and measured pH to decide if NaOH or HCl
            # should be added to adjust pH.
            if R_pH[0, :] > Measured_pH:
                Spec1 = 'Cl'
                S1_name = 'Chloride'
                S1_T_Index = Tot_Names.index(Spec1)
            elif R_pH[0, :] < Measured_pH:
                Spec1 = 'Na'
                S1_name = 'Sodium'
                S1_T_Index = Tot_Names.index(Spec1)
            else:
                print('No pH adjustment needed')        
        # On the second run, go through a sequential search of NaOH or HCl
        # guesses to find values that allow the code to converge.
        if jj == 1: 
            for aa in range(0,len(Add_1),1):
                # The purpose of this loop is to get a good starting place 
                # for the pH adjustment
                if status2 == True: 
                    break
                Flag = False
                status = False
                x0 = x0_1
                C_S1_1 = C_T_Array_set[S1_T_Index] + Add_1[aa]   
                C_T_Array[S1_T_Index] = C_S1_1
                C_T_solutes = np.append(C_T_Array, [C_Sugar]) #(µmol/L)
                C_T_solutes = C_T_solutes.reshape(len(C_T_solutes), 1) 
                C_H2O = WaterConc(C_EtOH, C_T_solutes) #(µmol/L)
                x0 = InitialGuesses(C_T_Array) #(µmol/L)
                Supersat_ratio_guess = np.ones([n_ppt, 1])
                x0 = np.append(x0, Supersat_ratio_guess)
                x0 = np.array([x0]).reshape(len(x0), 1)
                lb = np.ones([n_T+n_ppt, 1])*1e-20  # Lower bounds (µmol/L)
                ub = Bounds(C_T_Array)              # Upper bounds (µmol/L)
                run = '2.' + str(aa+1)
                param = [C_T_Array, C_H2O, C_EtOH, C_Sugar, Flag]
                # Max iterations, x tolerance, f tolerance, lower bounds,
                # upper bounds, current run number, wine label
                options = [100, 1e-9, 1e-9, lb, ub, run, Wine[ii]] 
                # Error trapping                 
                try: 
                    xsoln, status = NewtonsMethod(x0, options, param)#(µmol/L)
                except:    
                    print('Error occurred, but other runs will continue.')
                    pass  
                # Only proceed to the next loop if Newton's method converged. 
                if status == True: 
                    x0_2 = xsoln
                    for ss in range(0,len(Add_2),1):
                        status = False
                        # Initial guess for the concentration of NaOH or HCl
                        # (S1) used to adjust pH is based on total amount
                        # measured in wine + added S1
                        Flag = True
                        C_S1 = C_T_Array_set[S1_T_Index] + Add_2[ss]  
                        C_T_Array[S1_T_Index] = C_S1
                        x0 = np.append(x0_2, C_S1)
                        x0 = np.array([x0]).reshape(len(x0), 1)
                        lb = np.ones([n_T+n_ppt+1, 1])*1e-20  
                        lb[n_T+n_ppt] = C_T_Array_set[S1_T_Index]
                        ub = Bounds(C_T_Array)
                        # Upper bound for the total amount of NaOH or HCl used
                        # to adjust the pH. Currently set to 0.1 mol/L.
                        ub_S1 = 0.1*1e6 #(µmol/L)
                        ub = np.append(ub, ub_S1)
                        ub = np.array([ub]).reshape(len(ub), 1)
                        run = '2.' + str(aa+1) + '.' + str(ss+1)
                        C_T_solutes = np.append(C_T_Array, [C_Sugar])
                        C_T_solutes = C_T_solutes.reshape(len(C_T_solutes), 1)
                        C_H2O = WaterConc(C_EtOH, C_T_solutes)
                        param = [C_T_Array, C_H2O, C_EtOH, C_Sugar, Flag]
                        options = [100, 1e-9, 1e-9, lb, ub, run, Wine[ii]]
                        # Error trapping
                        try: 
                            xsoln, status = NewtonsMethod(x0, options, param)
                        except:    
                            print('Error occurred, but other runs will\
                                  continue.')
                            pass  
                        if  status == True:
                            status2 = True
                            R_Convergence[jj] = 'Converged'
                            break
        else:
            if jj == 3:
                C_T_Array[S1_T_Index] = R_C_T_Array[1, S1_T_Index]
            C_T_solutes = np.append(C_T_Array, [C_Sugar])           #(µmol/L)
            C_T_solutes = C_T_solutes.reshape(len(C_T_solutes), 1)  #(µmol/L)
            C_H2O = WaterConc(C_EtOH, C_T_solutes)  #(µmol/L)
            x0 = InitialGuesses(C_T_Array)  #(µmol/L)
            Supersat_ratio_guess = np.ones([n_ppt, 1])
            x0 = np.append(x0, Supersat_ratio_guess)
            x0 = np.array([x0]).reshape(len(x0), 1)
            lb = np.ones([n_T+n_ppt, 1])*1e-20  # Lower bounds (µmol/L)
            ub = Bounds(C_T_Array)              # Upper bounds (µmol/L)
            run = str(jj+1)
            param = [C_T_Array, C_H2O, C_EtOH, C_Sugar, Flag]
            # Max iterations, x tolerance, f tolerance, lower bounds, upper
            # bounds, current run number, wine label
            options = [100, 1e-9, 1e-9, lb, ub, run, Wine[ii]]
            xsoln, status = NewtonsMethod(x0, options, param)  #(µmol/L)
            if status == True: 
                R_Convergence[jj] = 'Converged'
        x0 = xsoln[0:n_T+n_ppt, 0]
        if jj == 0: 
            x0_1 = x0
        Conc_Array = np.append(xsoln[0:n_T], [C_EtOH, C_Sugar]) #(µmol/L)
        Conc_Array = Conc_Array.reshape(len(Conc_Array), 1)     #(µmol/L)
        C_solutes = np.append(xsoln[0:n_T], [C_Sugar])          #(µmol/L)
        C_solutes = C_solutes.reshape(len(C_solutes), 1)        #(µmol/L)
        rho = Conc_Array/1e6*N_A*1000                           #(number/m^3)
        y = ActCoeff(z, sigma, rho, C_H2O, C_EtOH, C_solutes)
        R_y[jj, :] = y[0:n_T, 0]
        R_C[jj, :] = xsoln[0:n_T, 0]
        R_SS[jj, :] = xsoln[n_T:n_T+n_ppt, 0]
        A = R_y[jj, :] * R_C[jj, :]  
        A_H = (A[H_Index])/1e6  
        pH_jj = -np.log10(A_H)
        R_pH[jj, :] = pH_jj
        Final_Ksp = CaTar_Ksp(A.T, C_EtOH, pH_jj)
        pKsp = -np.log10(Final_Ksp)
        R_pKsp[jj, :] = pKsp
        R_C_T_Array[jj, :] = C_T_Array[:, 0]
    W_C[ii] = R_C
    W_y[ii] = R_y
    W_SS[ii] = R_SS
    W_pH[ii] = R_pH
    W_pKsp[ii] = R_pKsp
    W_C_T_Array[ii] = R_C_T_Array
    W_Convergence[ii] = R_Convergence
data = [0]*L
data2 = [0]*L
for ii in range(0, L, 1):
    pH = list(sublist[ii] for sublist in W_pH)
    pH = np.reshape(pH, [Wine_num,])
    SS = list(sublist[ii] for sublist in W_SS)
    SS = np.reshape(SS, [Wine_num,])
    Con = list(sublist[ii] for sublist in W_Convergence)
    Con = np.reshape(Con, [Wine_num,])
    data[ii] = {'Wine': Wine, 'pH': pH, 'Supersaturation Ratio': SS,
                'Convergence Status': Con}
    data2[ii] = {'Supersaturation Ratio': SS,'Convergence Status': Con}
Info = {'Wines': str(W_A.tolist()), 'Code Description': Description}
df_info = pd.DataFrame(Info, index=[0])

df_0 = pd.DataFrame(data[0])
# df_1 = pd.DataFrame(data[1])
df_1 = pd.DataFrame(data2[1])
if Crystals == True:
    df_2 = pd.DataFrame(data[2])
    df_3 = pd.DataFrame(data[3])
    if Excel_f != False:
        with pd.ExcelWriter(Output_Name, engine='openpyxl') as writer:
            df_info.to_excel(writer, sheet_name='Summary', index=False)
            df_0.to_excel(writer, sheet_name='No Crystals No pH Adjustment',
                          index=False)
            df_1.to_excel(writer, sheet_name='No Crystals pH Adjustment',
                          index=False)
            df_2.to_excel(writer, sheet_name='Crystals No pH Adjustment',
                          index=False)
            df_3.to_excel(writer, sheet_name='Crystals pH Adjustment',
                          index=False)
    print('\n''No Crystals, No pH Adjustment:')
    print(df_0)
    print('\n''No Crystals, pH Adjustment:')
    print(df_1)
    print('\n''Crystals, No pH Adjustment:')
    print(df_2)
    print('\n''Crystals, pH Adjustment:')
    print(df_3,'\n') 
else:
    if Excel_f != False:
        with pd.ExcelWriter(Output_Name, engine='openpyxl') as writer:
            df_0.to_excel(writer, sheet_name='No Crystals Model pH',
                          index=False)
            df_1.to_excel(writer, sheet_name='No Crystals pH Adjustment',
                          index=False)
    # print('\n''No Crystals, No pH Adjustment:')
    # print(df_0)
    # print('\n''No Crystals, pH Adjustment:')
    # print(df_1.to_string(index=False))
    
if df_1['Convergence Status'].iloc[0] == 'Converged':
    print('A solution was found. Supersaturation ratio: ' 
          + str(round(df_1['Supersaturation Ratio'].iloc[0], 1)))
    pass
else: 
    print('A solution could not be found and the results may be inaccurate. Supersaturation ratio: ' 
          + str(round(df_1['Supersaturation Ratio'].iloc[0], 1)))
# ##############################################################################
# # Precipitation over time was added in as an example of how this was
# # implemented in other codes. This has just been added to the end for a single
# # wine.
#                           ### KINETICS OPTIONS ###
# # Parameters for ODE Solver
# t0 = 0     # Initial time (min)
# tf = 30    # Final time (min)
# t_step = 1 # Time step (min)
# ##############################################################################
# t = np.arange(t0, tf+t_step, t_step)
# n_steps = len(t)-1
# t_values = [t0, tf, t_step, n_steps, t]
# Wine_num = len(W_A)
# Wine = [0]*Wine_num
# Wine[0] = 'W' + str(W_A[0])
# Filename = Wine[0] + ' Data.xlsx'
# Excel_data = pd.read_excel(Filename, sheet_name=['Ion_Data', 'Ka_Data',
#                                                   'Size_Data', 'Conc_Data',
#                                                   'Precipitate_Data',
#                                                   'Other_Data'])
# Ion_Data = Excel_data.get('Ion_Data')
# Ka_Data = Excel_data.get('Ka_Data')
# Size_Data = Excel_data.get('Size_Data')
# Conc_Data = Excel_data.get('Conc_Data')
# Precipitate_Data = Excel_data.get('Precipitate_Data')
# Other_Data = Excel_data.get('Other_Data')
# C_T_Array_o = Conc_Data['Wine Concentration (mol/L)'].values.tolist()
# Crystal_T_Array_o = Conc_Data['Crystal Concentration (mol/L)'].values.tolist()
# C_T_symbol = Conc_Data['Symbol'].values.tolist()
# Symbol_A = Ion_Data['Anion Symbol'].dropna()
# Symbol_A = Symbol_A.values.tolist()
# Symbol_C = Ion_Data['Cation Symbol'].dropna()
# Symbol_C = Symbol_C.values.tolist()
# Symbol_A_ion = Ka_Data.iloc[:, 1]
# Symbol_A_ion = Symbol_A_ion.values.tolist()
# Ion_Arrays = IonArrays()
# Names_A = Ion_Arrays[0]       
# n_A = Ion_Arrays[1]           
# Names_C = Ion_Arrays[2]       
# n_C = Ion_Arrays[3]        
# Names_ip = Ion_Arrays[4]      
# Names_ip_all = Ion_Arrays[5]  
# z = Ion_Arrays[6]             
# n_T = Ion_Arrays[7]
# Names_all = Ion_Arrays[8]      
# sigma = Ion_Arrays[9]
# Mw_Array = Ion_Arrays[10]      
# C_T_Array_set = Ion_Arrays[11] 
# Tot_Names = Ion_Arrays[12]    
# C_T_Array_Crystals = Ion_Arrays[13]
# n_ppt = len(Precipitate_Data)
# H_Index = Names_all.index('H')
# z = np.append(z, [z_EtOH, z_Sugar])
# z = np.array([z]).reshape(len(z), 1)
# sigma = 1e-9*(np.array([sigma + sigma_non_dissociating])).T  #(m)
# Mw_Array = np.array([Mw_Array]).reshape(len(Mw_Array), 1)/1e6   #(g/µmol)
# Mw_solutes = np.append(Mw_Array, [Mw_Sugar])                    #(g/µmol)
# Mw_solutes = np.array([Mw_solutes]).reshape(len(Mw_solutes), 1) #(g/µmol)
# C_T_Array_set = np.array(C_T_Array_set)*1e6  # Not changed each run (µmol/L)
# C_T_Array = np.array([C_T_Array_set]).T      # Changed each run (µmol/L)
# n_Tot = len(C_T_Array)
# if Crystals == True:
#     C_T_Array_set = (C_T_Array_set + np.array([C_T_Array_Crystals])*1e6).T
#     C_T_Array = C_T_Array + np.array([C_T_Array_Crystals]).T*1e6
# R_C_T_Array = W_C_T_Array[0] 
# C_T_Array[S1_T_Index] = R_C_T_Array[1, S1_T_Index]
# C_T_Array_set[S1_T_Index] = R_C_T_Array[1, S1_T_Index]
# ppt_anion_names = Precipitate_Data['Anion Names']
# ppt_cation_names = Precipitate_Data['Cation Names']
# ppt_names = ppt_cation_names + ppt_anion_names
# EtOH_ABV = Other_Data['Ethanol Content (vol%)']  
# EtOH_ABV = float(EtOH_ABV[0])
# C_EtOH = EtOH_Conversion(EtOH_ABV, T_C)*1e6 #(µmol/L)
# C_Sugar = Other_Data['Sugar Content (g/L)']
# C_Sugar = float(C_Sugar[0])
# C_Sugar = C_Sugar/Mw_Sugar #(µmol/L)
# if Sugar_Change == True:
#     C_Sugar = 0
# R_C = np.zeros([n_steps, n_T+n_ppt]) 
# R_y = np.zeros([n_steps, n_T])
# R_P = np.zeros([n_steps+1, n_ppt])
# R_CT = np.zeros([n_steps+1, n_Tot])  
# R_CT[0, :] = C_T_Array[0:n_Tot, 0]
# Seed = Precipitate_Data['Initial Seed (mol/L)'].T*1e6  #(µmol/L)
# Seed = np.array([Seed])
# R_P[0, :] = Seed
# x = InitialGuesses(C_T_Array) #(µmol/L)
# Supersat_ratio_guess = np.array([[1], [0.1]])
# Supersat_ratio_guess = np.ones([n_ppt, 1])
# x = np.append(x, Supersat_ratio_guess)
# x = np.array([x]).reshape(len(x), 1)
# # Choose the method to solve the ODE. 'Eulers' or 'Heuns'.
# ODE_Method = 'Heuns'
# # Options: ODE method, relative tolerance, absolute tolerance
# ODE_options = [ODE_Method, 0.001, 1e-6]
# Results = ODE_Solver(func, t_values, R_P[0, :], ODE_options)
# R_P = Results
# H_Index = Names_all.index('H')
# A_H = R_C[:, H_Index] * R_y[:, H_Index]
# pH = -np.log10(A_H/1e6)
# Ca_Index = Tot_Names.index('Ca')
# Ca_Tot = R_CT[:, Ca_Index] #(µmol/L)
# Ca_Tot = Ca_Tot * Mw_Array[Ca_Index]*1000 #(mg/L)
# plt.figure(1, dpi=300)
# plt.rc('font', family='Calibri')
# plt.plot(t[0:n_steps], pH, label='Model', linewidth=1, color='k')
# plt.legend(frameon=False, fontsize=9)
# plt.xlabel('Time (min)', fontsize=10)
# plt.ylabel('pH', fontsize=10)
# plt.xticks(fontsize=9)
# plt.yticks(fontsize=9)
# plt.ylim(3, 4)
# plt.figure(2, dpi=300)
# plt.rc('font', family='Calibri')
# plt.plot(t, Ca_Tot, linewidth=1, color='k', label='Model')
# plt.legend(frameon=False, fontsize=9)
# plt.xlabel('Time (min)', fontsize=10)
# plt.ylabel('Total Solution Calcium Concentration (mg/L)', fontsize=10)
# plt.xticks(fontsize=9)
# plt.yticks(fontsize=9)
# plt.ylim(0, 100)
