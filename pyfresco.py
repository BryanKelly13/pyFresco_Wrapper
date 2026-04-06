import csv
import os
import subprocess
import re
import json
from pathlib import Path
import ADWA_potentials as adwa_pot
import sys




# Dictionaries used for formatting and parity assignment logic
l_dict = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6, 'j': 7}
frac_dict = {'0.5': 12, '1.5': 32, '2.5': 52, '3.5': 72, '4.5': 92, '5.5': 112, '6.5': 132, '7.5': 152} 

def write_sorted_from_fresco_output(infile: str, out_dir: str, out_basename: str) -> str:
    """
    Extract angular distribution block from a FRESCO text output (e.g. fort.16 renamed)
    and write a two-column .sorted file (angle, xsec) to out_dir.

    Returns the path to the written .sorted file.
    """
    angles = []
    cross_sections = []

    start_line = "@s1"
    end_line = "END"
    end_count = 0
    lines_after_start = 0
    between_lines = False

    with open(infile, "r") as f:
        for raw in f:
            line = raw.strip()

            if start_line in line:
                between_lines = True
                lines_after_start = 2
                continue

            if line == end_line:
                if end_count == 0:
                    end_count += 1
                    continue
                else:
                    between_lines = False
                    break

            if not between_lines:
                continue

            if lines_after_start > 0:
                lines_after_start -= 1
                continue

            # Expect 2 columns: angle xsec
            parts = line.split()
            if len(parts) < 2:
                continue
            ang, xsec = parts[0], parts[1]
            angles.append(ang)
            cross_sections.append(xsec)

    os.makedirs(out_dir, exist_ok=True)
    storefile = os.path.join(out_dir, f"{out_basename}.sorted")
    with open(storefile, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(zip(angles, cross_sections))

    print(f"Successfully written file {storefile}")
    return storefile
    
def print_optical_potentials(deuteron_pot, proton_pot, config):
    """
    Pretty-printer for ADWA optical model parameters.
    """

    def _print_block(label, particle, pot):
        print("=" * 60)
        print(f"{label} + {particle} Optical Model Parameters")
        print("=" * 60)

        fmt = "{:<8s} = {:>8.3f}"
        fmtV = "{:<8s} = {:>8.3f}"

        # Central real + imaginary
        print(fmtV.format("V",   pot["V"]))
        print(fmt.format("r0",  pot["r0"]))
        print(fmt.format("a",   pot["a"]))
        print(fmtV.format("Vi",  pot["Vi"]))
        print(fmt.format("ri0", pot["ri0"]))
        print(fmt.format("ai",  pot["ai"]))

        # Surface imaginary
        print(fmtV.format("Vsi",  pot["Vsi"]))
        print(fmt.format("rsi0", pot["rsi0"]))
        print(fmt.format("asi",  pot["asi"]))

        # Spin–orbit
        print(fmtV.format("Vso",  pot["Vso"]))
        print(fmt.format("rso0", pot["rso0"]))
        print(fmt.format("aso",  pot["aso"]))
        print(fmtV.format("Vsoi",  pot["Vsoi"]))
        print(fmt.format("rsoi0", pot["rsoi0"]))
        print(fmt.format("asoi",  pot["asoi"]))

        # Coulomb radius
        print(fmt.format("rc0", pot["rc0"]))
        print()

    _print_block(config["label_in"],  "d", deuteron_pot)
    _print_block(config["label_out"], "p", proton_pot)


def load_reaction_config(config_path):
    """
    Loads the reaction configuration from a JSON file.

    Parameters:
    - config_path (str or Path): Path to the JSON config file.

    Returns:
    - dict: A dictionary containing all reaction parameters.
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, "r") as file:
        config = json.load(file)

    return config
    
def getInput():
    """
    This scrapes the input_generator.inp file for the energies, jpi, and transfer configuration
    and returns them in a list. More than one state can be calculated at once, so however many
    states are given in input_generator.inp will be calculated, and fresco will be ran for each of them.
    """
    inputFileName = 'input_generator.inp'
    inputFile = open(inputFileName, 'r')
    energies = []
    ns       = []
    ls       = []
    js_transfer = []
    js_finalstate = []
    while True:
        line = inputFile.readline()
        if line.startswith('#'):
            continue
        if line == '':
            break
        args = line.split()
        energies.append(float(args[0]))
        js_finalstate.append(float(args[1]))
        ns.append(args[2])
        ls.append(args[3])
        js_transfer.append(args[4])
    
    return energies, ns, ls, js_transfer, js_finalstate

def runFresco(inFile, outfile, dir, string_suffix, transfer_info, evenMassFinal, config):
    """
    This function runs the FRESCO executable with the given input file and handles the output.
    It renames the output file to avoid overwriting in the next iteration, if there are multiple states in the input file.
    Make sure to adjust the path to the FRESCO executable in the reaction_config file.

    The fort.16 file is where the inelastic channel cross-section is stored. 
    A fresco .fri (input) and .fro (output) file, and {nucleus}_{energy}_{n}{l}{jpi}.txt (cross-section output) file is created for each state.
    The inelastic channel is the second partition, the elastic is always calculated and is the first set of cross-section data.
    """
    # Run the FRESCO executable with the input file
    executable_path = config['executable_path']
    command_string = f'{executable_path} < {inFile} > {outfile}'
    command = [command_string] 
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as error:
        print(f'Error running the command {error}')
    
    # take the cross-section output file & rename it so it does not get overwritten in the next iteration
    x_sec_outfile = f'{dir}/fort.16'
    renamed_outfile = f"{dir}/{config['label_out']}_{string_suffix}.txt"

    try:
        os.rename(x_sec_outfile, renamed_outfile)
        print(f"File '{x_sec_outfile}' has been renamed to '{renamed_outfile}'.")
        fresco_out_dir = config.get("fresco_output_dir")  
        if fresco_out_dir is None:
            print("WARNING: config['fresco_output_dir'] not set, skipping .sorted creation.")
        else:
            # Use same naming root as your pipeline expects
            # Example: 48Ti_7155_1f52_4+.sorted
            sorted_basename = f"{config['label_out']}_{string_suffix}"
            try:
                write_sorted_from_fresco_output(renamed_outfile, fresco_out_dir, sorted_basename)
            except Exception as e:
                print(f"WARNING: failed to write .sorted for {renamed_outfile}: {e}")




    except FileNotFoundError:
        print(f"Error: File '{x_sec_outfile}' not found.")
    except FileExistsError:
        print(f"Error: File '{renamed_outfile}' already exists.")
    except Exception as error:
        print(f"An error occurred: {error}")
    if evenMassFinal:
        print(f'Calculated the transfer to the Jpi={int(transfer_info["js_final"])}{transfer_info["parity"]}, {transfer_info["n"]}{transfer_info["l"]}{transfer_info["js_transfer"]} configuration')
    else:
        print(f'Calculated the transfer to the Jpi={frac_dict[str(transfer_info["js_final"])]}{transfer_info["parity"]}, {transfer_info["n"]}{transfer_info["l"]}{transfer_info["js_transfer"]} configuration')
    print(f"Successfully ran fresco & created cross-section file {renamed_outfile}")

def createInputFile(energies, ns, ls, js_transfer, js_finalstate, deuteron_pot, proton_pot, config):
    """
    This handles the creation of the input file for FRESCO.
    The input file is created based on the parameters in the reaction_config.json file, and uses the ADWA_potentials.py file to
    calculate the deuteron and proton potentials.

    Make sure that you updated the path to the FRESCO executable file in the reaction_config.json.

    DO NOT ADJUST SPACING, the f-string below looks crazy with all parameters filled but the formatting is very particular 
    """

    # Gather all necessary information to create the input file
    Z_fmt = f"{config['Z']:4.1f}"
    mass_in = f"{config['mass_in']:7.3f}"
    mass_out = f"{config['mass_out']:7.3f}"
    gs_spin_fmt = f"{config['gs_spin']:4.1f}"     # Ensures correct spacing and decimal format
    gs_parity_fmt = f"{config['gs_parity']:>2}"
    th_min = f"{config['angle_min']:4.1f}"
    th_max = f"{config['angle_max']:4.1f}"
    th_step = f"{config['angle_step']:3.1f}"
    at = float(config['AT'])
    at_fmt = f"{float(config['AT']):8.1f}"
    residual_mass = float(config['residual_mass'])
    residual_mass_fmt = f"{float(config['residual_mass']):8.1f}"
    match = re.match(r"(\d+)([A-Za-z]+)", config['label_out'])
    mass, element = match.groups()
    label_out_reversed = element.lower() + mass
    evenMassFinal = True if int(config['label_out'][0:2]) % 2 == 0 else False
    beam_energy = f"{float(config['beam_energy']):4.1f}"  
    
    # This is the loop that creates the input file for each state, which is based on the number of
    # states in the input_generator.inp file.
    for i in range(len(energies)):
        q = f"{config['Q_value']:6.4f}" # Q value
        e =  f"{energies[i]:5.3f}" # Energy of the state
        be = f"{(config['binding_energy'] - float(e)):6.4f}" # Effective binding energy
        energy = float(e) * 1000 # Convert MeV to keV
        l = l_dict[ls[i]]
        if l % 2 == 0:
            final_parity = config['gs_parity']
        else:
            final_parity = '+1' if config['gs_parity'] == '-1' else '-1'
        transfer_info = {
            'energy': int(energy),
            'n': ns[i],
            'l': ls[i],
            'js_transfer': frac_dict[js_transfer[i]],
            'js_final': js_finalstate[i],
            'parity': final_parity[0],  
        }
        
        fri=f'''{config['reaction']}, {js_finalstate[i]}{final_parity[0]} {e} MeV {ns[i]}{ls[i]}{js_transfer[i]}
0.10    55.0    0.20    0.20    30.0    -6.0
 00. 20.  +.00   F F
0  {th_min}     {th_max}  {th_step}  1
0.0    0 1   1 1  48          .000    0.   0.001
 1 1 0 0 2 3 0 0-3 1 0 0 1
2H      2.0141  1.0        1  {config['label_in']:<8}{mass_in} {Z_fmt}    0.0000
1.0   +1 0.0               1 {gs_spin_fmt}   {gs_parity_fmt} 0.000
1H      1.0078  1.0        1  {config['label_out']:<8}{mass_out} {Z_fmt}    {q}
0.5   +1 0.0               2  {js_finalstate[i]}   {final_parity} {e}

  1 0  0{at_fmt}     0.0   {deuteron_pot['rc0']:5.3f}
  1 1  0   {deuteron_pot['V']:5.2f}   {deuteron_pot['r0']:5.3f}   {deuteron_pot['a']:5.3f}   {deuteron_pot['Vi']:5.3f}   {deuteron_pot['ri0']:5.3f}   {deuteron_pot['ai']:5.3f}
  1 2  0                           {deuteron_pot['Vsi']:5.2f}   {deuteron_pot['rsi0']:5.3f}   {deuteron_pot['asi']:5.3f}
  1 3  0   {deuteron_pot['Vso']:5.2f}   {deuteron_pot['rso0']:5.3f}   {deuteron_pot['aso']:5.3f}   {deuteron_pot['Vsoi']:5.2f}   {deuteron_pot['rsoi0']:5.3f}   {deuteron_pot['asoi']:5.3f}
  2 0  0{residual_mass_fmt}     0.0   {proton_pot['rc0']:5.3f}
  2 1  0   {proton_pot['V']:5.2f}   {proton_pot['r0']:5.3f}   {proton_pot['a']:5.3f}   {proton_pot['Vi']:5.3f}   {proton_pot['ri0']:5.3f}   {proton_pot['ai']:5.3f}
  2 2  0                           {proton_pot['Vsi']:5.2f}   {proton_pot['rsi0']:5.3f}   {proton_pot['asi']:5.3f}
  2 3  0   {proton_pot['Vso']:5.2f}   {proton_pot['rso0']:5.3f}   {proton_pot['aso']:5.3f}   {proton_pot['Vsoi']:5.2f}   {proton_pot['rsoi0']:5.3f}   {proton_pot['asoi']:5.3f}
  3 0  0    1.00            1.25
  3 1  5    1.00            1.00
  3 3  5    1.00            1.00
  3 4  5    1.00            1.00
  3 7  5    1.00            1.00
  4 0  0{residual_mass_fmt}     0.0    1.25
  4 1  0    50.0    1.25    0.65
  4 3  0     6.0    1.25    0.65
0
  1  2 1 2-1 3   1 0 2 0.5 1 1.5 1  3  0  2.2260  0  3  0
  3    1 2-2 0   {ns[i]} {l}   0.5   {js_transfer[i]}    4  0  {be}  1  0  0

  -2   1   7 0-1 0
       1   1   1   1  1.0000
       2   1   1   3  1.0000
0
   0   1   1
{beam_energy}
EOF
'''
        # Handles the output logic, depending on whether you land on even/odd nucleus depends on interpretation of the js_finalstate
        # Whether your Jpi is integer or half-integer form
        if evenMassFinal: # landing on even nucleus
            string_suffix = f'{transfer_info["energy"]}_{transfer_info["n"]}{transfer_info["l"]}{transfer_info["js_transfer"]}_{int(transfer_info["js_final"])}{transfer_info["parity"]}'
            fileName = f'{label_out_reversed}dp_adwa_{string_suffix}.fri'
            inFile = fileName
            outfile = f'{label_out_reversed}dp_adwa_{string_suffix}.fro'
        else:
            string_suffix = f'{transfer_info["energy"]}_{transfer_info["n"]}{transfer_info["l"]}{transfer_info["js_transfer"]}_{frac_dict[str(transfer_info["js_final"])]}{transfer_info["parity"]}'
            fileName = f'{label_out_reversed}dp_adwa_{string_suffix}.fri'
            inFile = fileName
            outfile = f'{label_out_reversed}dp_adwa_{string_suffix}.fro'
        if not os.path.exists(inFile):
            with open(inFile, 'w') as file:
                file.write(fri)
        runFresco(inFile, outfile, os.getcwd(), string_suffix, transfer_info, evenMassFinal, config)

        
def main():
    config = load_reaction_config("reaction_config.json")
    beam_energy, zt, at, residual_mass = config["beam_energy"], config["Z"], config["AT"], config["residual_mass"]
    deuteron_pot = adwa_pot.Wales_Johnson_deutron_AWDA(beam_energy, zt, at, params={})
    proton_pot = adwa_pot.koning_delaroche_proton_potential(beam_energy, zt, residual_mass, params={})

    # Utility print function if you pass optional arg. 'print_params'
     # Optional print-only mode
    if len(sys.argv) > 1 and sys.argv[1] == "print_params":
        print_optical_potentials(deuteron_pot, proton_pot, config)
        return

    energies, ns, ls, js_transfer, js_finalstate = getInput()
    createInputFile(energies, ns, ls, js_transfer, js_finalstate, deuteron_pot, proton_pot, config)

if __name__ == "__main__":
    main()

