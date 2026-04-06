import math

def round_params(params: dict[str, float], ndigits: int = 4) -> dict[str, float]:
    return {k: round(v, ndigits) if isinstance(v, float) else v for k, v in params.items()}

def koning_delaroche_neutron_potential(
    E: float, zt: int, at: int, params: dict[str, float]
):
    """Koning-Delaroche neutron scattering optical model potential
    From  Koning, A.J., Delaroche, J.P., "Local and global nucleon optical models from 1 keV to 200 MeV", Nuclear Physics A, 713, 2003
    https://doi.org/10.1016/S0375-9474(02)01321-0

    Parameters
    ----------
    E: float
        The projectile energy in MeV
    zt: int
        The target Z
    at: int
        The target A
    """
    params = {}
    nt = at - zt
    v1 = 59.30 - 21.0 * (nt - zt) / at - 0.024 * at
    v2 = 0.007228 - 1.48e-6 * at    
    v3 = 1.994e-5 - 2.0e-8 * at
    v4 = 7.0e-9

    w1 = 12.195 + 0.0167 * at         
    w2 = 73.55 + 0.0795 * at

    d1 = 16.0 - 16.0 * (nt - zt) / at          
    d2 = 0.0180 + 0.003802 / (1.0 + math.exp((at - 156.0) / 8.0))
    d3 = 11.5

    vso1 = 5.922 + 0.0030 * at
    vso2 = 0.0040

    wso1 = -3.1
    wso2 = 160

    ef = -11.2814 + 0.02646 * at        
    delta_e = E - ef
    params["V"] = v1 * (
        1.0 - v2 * delta_e + v3 * delta_e**2.0 - v4 * delta_e**3.0
    )
    params["Vi"] = w1 * delta_e**2.0 / (delta_e**2.0 + w2**2.0)
    params["Vsi"] = (
        d1 * delta_e**2.0 / (delta_e**2.0 + d3**2.0) * math.exp(-1.0 * d2 * delta_e)
    )
    params["Vso"] = vso1 * math.exp(-1.0 * vso2 * delta_e)
    params["Vsoi"] = wso1 * delta_e**2.0 / (delta_e**2.0 + wso2**2.0)
    params["r0"] = 1.3039 - 0.4054 * at ** (-1/3)
    params["ri0"] = params["r0"]
    params["rsi0"] = 1.3424 - 0.01585 * at ** (1/3)
    params["rso0"] = 1.1854 - 0.647 * at ** (-1/3)
    params["rsoi0"] = params["rso0"]
    params["a"] = 0.6778 - 1.487e-4 * at
    params["ai"] = 0.6778 - 1.487e-4 * at
    params["asi"] = 0.5446 - 1.656e-4 * at
    params["aso"] = 0.59
    params["asoi"] = 0.59
    

    return params
    # rounded_params = round_params(params)
    # return rounded_params

def koning_delaroche_proton_potential(
    E: float, zt: int, at: int, params: dict[str, float]
):
    """Koning-Delaroche proton scattering optical model potential

    From  Koning, A.J., Delaroche, J.P., "Local and global nucleon optical models from 1 keV to 200 MeV", Nuclear Physics A, 713, 2003
    https://doi.org/10.1016/S0375-9474(02)01321-0

    Parameters
    ----------
    E: float
        The projectile energy in MeV
    zt: int
        The target Z
    at: int
        The target A
    """
    params = {}
    nt = at - zt

    
    v1 = 59.30 + 21.0 * (nt - zt) / at - 0.024 * at
    v2 = 0.007067 + 4.23e-6 * at
    v3 = 1.729e-5 + 1.136e-8 * at
    v4 = 7.0e-9

    w1 = 14.667 + 0.009629 * at
    w2 = 73.55 + 0.0795 * at

    d1 = 16.0 + 16.0 * (nt - zt) / at
    d2 = 0.0180 + 0.003802 / (1.0 + math.exp((at - 156.0) / 8.0))
    d3 = 11.5

    vso1 = 5.922 + 0.0030 * at
    vso2 = 0.0040

    wso1 = -3.1
    wso2 = 160

    ef = -8.4075 + 0.01378 * at
    rc = 1.198 + 0.697 * at ** (-2/3) + 12.994 * at ** (-5/3)
    vc = 1.73 / rc * zt * at ** (-1/3)


    delta_e = E - ef
    params["V"] = v1 * (
        1.0 - v2 * delta_e + v3 * delta_e**2.0 - v4 * delta_e**3.0
    ) + vc * v1 * (v2 - 2.0 * v3 * delta_e + 3.0 * v4 * delta_e**2.0)
    params["Vi"] = w1 * delta_e**2.0 / (delta_e**2.0 + w2**2.0)
    params["Vsi"] = (
        d1 * delta_e**2.0 / (delta_e**2.0 + d3**2.0) * math.exp(-1.0 * d2 * delta_e)
    )
    params["Vso"] = vso1 * math.exp(-1.0 * vso2 * delta_e)
    params["Vsoi"] = wso1 * delta_e**2.0 / (delta_e**2.0 + wso2**2.0)
    params["r0"] = 1.3039 - 0.4054 * at ** (-1/3)
    params["ri0"] = params["r0"]
    params["rsi0"] = 1.3424 - 0.01585 * at ** (1/3)
    params["rso0"] = 1.1854 - 0.647 * at ** (-1/3)
    params["rsoi0"] = params["rso0"]
    params["rc0"] = rc
    params["a"] = 0.6778 - 1.487e-4 * at
    params["ai"] = 0.6778 - 1.487e-4 * at
    params["asi"] = 0.5187 + 5.205e-4 * at
    params["aso"] = 0.59
    params["asoi"] = 0.59
    
    return params
    # rounded_params = round_params(params)
    # return rounded_params

def Wales_Johnson_deutron_AWDA(
    E: float, zt: int, at: int, params: dict[str, float]
):

    """
    Wales and Johnson deuteron scattering optical model potential
    From G.L. Wales and R.C. Johnson, "Deuteron break-up effects in (p, d) reactions at 65 MeV", Nuclear Physics A, 274, 1976
    https://doi.org/10.1016/0375-9474(76)90234-7    
    """
    E_deuteron = E / 2.0
    params_protons, params_neutrons = {}, {}
    params_protons = koning_delaroche_proton_potential(E_deuteron, zt, at, params)
    params_neutrons = koning_delaroche_neutron_potential(E_deuteron, zt, at, params)

    # print(params_protons)
    # print(params_neutrons)

    params = {}
    # For central potential terms
    r2 = 0.6 
    daBar = (5.0/(14.0* math.pi**2)) * (r2 / params_neutrons["a"])
    aBar = params_protons["a"] + daBar
    daBar_prime = (3.0/ (10.0 * math.pi**2)) * (r2 / params_neutrons["a"])
    aBar_prime = daBar_prime + params_neutrons["a"]
    dV = ( params_protons["V"] + params_neutrons["V"] ) * math.pi**2 / params_protons["r0"] ** 2 * (params_neutrons["a"] ** 2 - aBar ** 2) / (1 + (math.pi * aBar / params_protons["r0"])**2 )
    dVi =  (- daBar_prime / aBar_prime) * ( 1.0 + 2.0/3.0*(math.pi*params_protons["a"]/ params_protons["r0"])**2 )


    
    # For surface terms
    da_VSi = (3.0/(10.0 * math.pi**2)) * (r2 / params_neutrons["asi"])
    aBar_Vsi = da_VSi + params_neutrons["asi"]
    dVSi = ( -da_VSi / aBar_Vsi) * ( 1.0 + 2.0/3.0*(math.pi*params_neutrons["asi"]/ params_neutrons["rsi0"])**2 )
    

    # For spin-orbit terms
    daBar_Vso = (5.0/(14.0 * math.pi**2)) * (r2 / params_protons["aso"])
    aBar_Vso = params_protons["aso"] + daBar_Vso
    daBar_prime_Vso = (3.0/ (10.0 * math.pi**2)) * (r2 / params_protons["aso"])
    aBar_prime_Vso = daBar_prime_Vso + params_protons["aso"]
    dVso = ( params_protons["Vso"] + params_neutrons["Vso"] ) * math.pi**2 / params_protons["rso0"] ** 2 * (params_protons["aso"] ** 2 - aBar_Vso ** 2) / (1 + (math.pi * aBar_Vso / params_protons["rso0"])**2 )
    dVsoi = ( - daBar_prime_Vso / aBar_prime_Vso) * ( 1.0 + 2.0/3.0*(math.pi*params_protons["aso"]/ params_protons["rso0"])**2 )

    # Central potential
    params["V"] = (
        params_protons["V"] + params_neutrons["V"] + dV
    )
    params["Vi"] = (
        params_protons["Vi"] + params_neutrons["Vi"] + dVi
    )
    params["r0"] = params_protons["r0"]
    params["a"] = aBar
    params["ri0"] = params["r0"]
    params["ai"] = aBar_prime

    # Surface terms
    params["Vsi"] = ( 
        params_protons["Vsi"] + params_neutrons["Vsi"] + dVSi
    )

    params["rsi0"] = params_neutrons["rsi0"]
    params["asi"] = aBar_Vsi

    # Spin-orbit terms
    params["Vso"] = (
        params_protons["Vso"] + params_neutrons["Vso"] + dVso
    )

    params["Vsoi"] = (
        params_protons["Vsoi"] + params_neutrons["Vsoi"] + dVsoi
    )  
    params["rso0"] = params_neutrons["rso0"]

    params["aso"] = aBar_Vso
    params["rsoi0"] = params["rso0"]
    params["asoi"] = aBar_prime_Vso

    # Lastly, need the charge radius for target mass, depends on proton optical potential w/ mass, no energy dependence
    params["rc0"] = params_protons["rc0"]

    return params
    # rounded_params = round_params(params)
    # return rounded_params


def li_liang_cai_triton_potential(E: float, zt: int, at: int, params: dict[str, float]):
    """Global triton optical model potential

    From: Li, Liang, Cai, (2007) E < 40 | All masses | 48 < A < 232 | Tritons
    Parameters
    ----------
    E: float
        The projectile energy in MeV
    zt: int
        The target Z
    at: int
        The target A
    params: dict[str, float]
        The dictionary of optical model parameters to be filled out
    """
    nt = at - zt
    a3 = at ** (1.0 / 3.0)
    
    params["V"] = 137.6 - 0.1456 * E + 0.0436 * E**2 + 4.3751 * (nt - zt) / at + 1.0474 * zt / a3
    params["r0"] = 1.1201 - 0.1504 / a3
    params["a"] = 0.6833 + 0.0191 * a3

    params["Vi"] = 7.383 + 0.5025 * E - 0.0097 * E**2
    params["ri0"] = 1.3202 - 0.1776 / a3
    params["ai"] = 1.119 + 0.01913 * a3

    params["Vsi"] = 37.06 - 0.6451 * E - 47.19 * (nt - zt) / at
    params["rsi0"] = 1.251 - 0.4622 / a3
    params["asi"] = 0.8114 + 0.01159 * a3

    params["Vso"] = 1.9029
    params["rso0"] = 0.46991 + 0.1294 / a3
    params["aso"] = 0.3545 - 0.0522 * a3

    params["Vsoi"] = 0.0
    params["rsoi0"] = 0.0
    params["asoi"] = 0.0

    params["rc0"] = 1.422
    return params


def main():
    # Example usage
        E = 16.0
        zt = 92
        at = 238
        params = {}
        
        a1 = koning_delaroche_proton_potential(E, zt, at, params)
        print(a1)
        # a2 = koning_delaroche_neutron_potential(E, zt, at, params)
        # a3 = Wales_Johnson_deutron_AWDA(E, zt, at)
        # print(a3)

        # a4 = li_liang_cai_triton_potential(E, zt, at, params)
        # print(a4)
if __name__ == "__main__":
    main()
        

