# pyFresco

pyFresco is a Python wrapper script that calculates optical model potential parameters and creates FRESCO input cards, which are then run using the FRESCO nuclear reaction code.
Features

    Calculates proton & neutron optical potentials using the Koning-Delaroche global optical model.

    Computes deuteron optical potentials using the Adiabatic Distorted Wave Approximation (ADWA).

    Automatically generates and runs FRESCO input files.

    Extracts and plots angular distributions.
    Python >= 3.8 (recommended 3.9+)

## Setup

Create and edit a reaction_config.json file. Example:
```
{
  "reaction": "49Ti(d,p)50Ti",
  "beam_energy": 16.0,
  "Z": 22,
  "AT": 49,
  "residual_mass": 50,
  "mass_in": 48.9478,
  "mass_out": 49.9447,
  "label_in": "49Ti",
  "label_out": "50Ti",
  "Q_value": 8.7146,
  "binding_energy": 10.939,
  "gs_spin": 3.5,
  "gs_parity": "-1",
  "angle_min": 15.0,
  "angle_max": 65.0,
  "angle_step": 0.1,
  "executable_path": "./bin/fresco"
}
```
This of course implies you have a working version of FRESCO on your local machine, update the executable_path to the directory where you have the fresco executable script located.

Also prepare an input_generator.inp file with information for each state of interest:
Excitation energy, spin (J), and transfer configuration (n, l, j).

```
# E	Jpi	n	l	j_transfer
4.172	3	2	p	1.5
```
**
Note any half-integer spins, should be input as a decimal value, not using the fraction value (i.e. 1.5 instead of 3/2)**
## Running the Code
Run the core script:
```
python pyfresco.py
```
For each state listed in the .inp file, this will:

    Create .fri (input) and .fro (output) files from running FRESCO

    Output angular distribution files with the naming convention:

{nucleus}{excitation_energy}{transfer_config}_{JPi}.txt

## Post-Processing & Plotting

To extract and sort FRESCO output:
```
python getFresco_output.py {nucleus}_{excitation_energy}_{transfer_config}_{JPi}.txt
```

This creates a .sorted version:

{nucleus}_{excitation_energy}_{transfer_config}_{JPi}.sorted

To plot the angular distribution:
```
python plot_fresco.py {nucleus}{excitation_energy}{transfer_config}_{JPi}.sorted
```

### Example
An example case for the reaction 49Ti(d,p)50Ti populating a 3⁺ state at 4.172 MeV is included in the default configuration files.
