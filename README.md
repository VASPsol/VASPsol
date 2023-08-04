# VASPsol++: A framework for implementing complex continuum fluid models in VASP density functional theory calculations

VASPsol++ is a framework for implementing complex continuum fluid models within density functional theory calculations performed using the Vienna Ab initio Simulation Package (VASP) ([https://www.vasp.at/](url)). It is being actively developed within the Plaisance group ([https://www.lsu.edu/eng/che/people/faculty/plaisance.php](url)) at Louisiana State University and has its origins in the VASPsol code ([https://github.com/henniggroup/VASPsol](url)) developed in the Hennig group ([https://hennig.mse.ufl.edu/](url)) at the University of Florida.

VASPsol++ adds a nonlocal and nonlinear implicit electrolyte model to the linear polarizable continuum model contained in the original VASPsol code. Additionally, it is written in a modular format that allows for easy addition of any new continuum solvation models that are developed in the future.

## Key features of the nonlinear+nonlocal model

* Nonlinear dielectric and ionic responses to model the high electric fields present at charged electrodes
* Uses a nonlocal cavity definition to prevent unphysical electrolyte leakage into small regions and to allow for separate dielectric and ionic cavities
* Efficient and robust solver for the nonlinear Poisson-Boltzmann equation based on Newton's method with a line search
* Option to performing calculations at constant potential rather than constant number of electrons (grand canonical DFT)
* Only slightly higher computational cost than the linear VASPsol model


## Installation

VASPsol++ is mostly implemented in a single fortran file, solvation.F, that can be found in the src/ directory of the repository. Additionally, a patch file is required to make modifications to some of the original VASP source files. These patch files are version specific and are located in the src/patches director of the repository. Currently, patch files are available only for VASP 5.4.4 although these could possibly work with other versions (no guarantees).

The instructions below assume you have downloaded the VASP source files in a directory `<VASP_SRC>` and configured the necessary makefiles specified in the VASP installation instructions.

1. Download the `src/solvation.F` file from the repository and copy it to `<VASP_SRC>/src/`, replacing the skeleton `solvation.F file` already there.

2. Download the appropriate patch file from src/patches and copy to `<VASP_SRC>/`. There are actually two patch files for each version of VASP. The first, named `vaspsol++-vtst-vasp_<version>.patch`, should be used if you are using the VTST add-on developed by the Henkelman group ([https://theory.cm.utexas.edu/vtsttools/](url)). Otherwise, use the patch file named `vaspsol++-vasp_<version>.patch`.

3. Apply the patch by running (where `<patch_file>` is the name of the patch file): `patch -p1 < <patch_file>`

4. Compile VASP as you normally would

## Input

### General

Parameters available for both <b>ISOL</b>=1 and <b>ISOL</b>=2:

* <b>LSOL</b>
    * <b>LSOL</b> = .TRUE. \
    turn on solvation
    * <b>LSOL</b> = .FALSE. (default) \
    turn off solvation
* <b>ISOL</b>
    * <b>ISOL</b> = 1 (default) \
    use the linear+local solvation model from the original VASPsol implementation
    * <b>ISOL</b> = 2 \
    use the nonlinear+nonlocal solvation model
* <b>LSOL_SCF</b>
    * <b>LSOL_SCF</b> = .TRUE. (default) \
    include solvation in the SCF cycle
    * <b>LSOL_SCF</b> = .FALSE. (not recommended) \
    do not include solvation in the SCF cycle, only as a correction at the end

Parameters only available for <b>ISOL</b>=2:

* <b>SOLTEMP</b> = 298 (default) \
Set the simulation temperature in K
* <b>A_K</b> = 0.125 (default) \
Set the smoothing length (&#197;) for eliminating FFT errors. The default value works well for the standard FFT grid used to represent the charge density in VASP.

### Cavity definition

Parameters available for both <b>ISOL</b>=1 and <b>ISOL</b>=2:

* <b>NC_K</b> = 0.0025 (default, <b>ISOL</b>=1) \
<b>NC_K</b> = 0.015 (default, <b>ISOL</b>=2) \
Cut-off charge density (e/&#197;<sup>3</sup>) for determining the vdW cavity
* <b>SIGMA</b> = 0.6 (default) \
Smoothness of the cavity transition. Increasing this value will make the cavity transition smoother, decreasing will make it sharper.
* <b>TAU</b> = 5.25e-4 (default, <b>ISOL</b>=1) \
<b>TAU</b> = 8.79e-4 (default, <b>ISOL</b>=2) \
Effective surface tension (eV/&#197;<sup>2</sup>) for computing the cavity formation free energy

Parameters only available for <b>ISOL</b>=2:

* <b>R_CAV</b> = 0. (default) \
Offset (&#197;) for determining the solute surface area used for calculating the cavity formation free energy

### Solvent specification

Parameters available for both <b>ISOL</b>=1 and <b>ISOL</b>=2:

* <b>EB_K</b> = 78.4 (default) \
Bulk dielectric constant of the solvent

Parameters only available for <b>ISOL</b>=2:

* <b>LNLDIEL</b>
    * <b>LNLDIEL</b> = .TRUE. (default)\
    use a nonlinear dielectric screening model
    * <b>LNLDIEL</b> = .FALSE. \
    use a linear dielectric screening model
* <b>EPSILON_INF</b> = 1.78 (default) \
Bulk optical dielectric constant of the solvent
* <b>N_MOL</b> = 0.0335 (default) \
Density of solvent molecules in the bulk (&#197;<sup>-3</sup>)
* <b>P_MOL</b> = 0.50 (default) \
Dipole moment of a solvent molecule (e-&#197;)
* <b>R_SOLV</b> = 1.40 (default) \
Solvent radius (&#197;) for constructing the solvent cavity
* <b>R_DIEL</b> = 1.00 (default) \
Dielectric radius (&#197;) for constructing the dielectric cavity
* <b>R_B</b> = A_K (default) \
Smearing length (&#197;) for the bound charge, used to reduce FFT truncation error

### Electrolyte specification

Parameters available for both <b>ISOL</b>=1 and <b>ISOL</b>=2:

* <b>LAMBDA_D_K</b> = 0 (default) \
Set the Debye screening length (&#197;) for <b>ISOL</b>=1. Setting <b>LAMBDA_D_K</b> $\le$ 0 will turn off ionic screening. Can also be used with <b>ISOL</b>=2 as an alternative to <b>C_MOLAR</b>.

Parameters only available for <b>ISOL</b>=2:

* <b>LNLION</b>
    * <b>LNLION</b> = .TRUE. (default)\
    use a nonlinear ionic screening model
    * <b>LNLION</b> = .FALSE. \
    use a linear ionic screening model
* <b>C_MOLAR</b> = 0.0 (default) \
Concentration of the electrolyte (mol/L)
* <b>ZION</b> = 1.0 (default) \
Electrolyte valency
* <b>D_ION</b> = 2<sup>5/6</sup>R_ION (default) \
Packing diameter of the ions (&#197;), defaults to the close-packed diameter corresponding to <b>R_ION</b>
* <b>R_ION</b> = R_SOLV (default) \
Ionic radius (&#197;) for constructing the ionic cavity

### Constant potential calculations (<b>ISOL</b>=2 only)

VASPsol++ can perform constant potential calculations where the number of electrons in varied. This is done by specifying the <b>EFERMI_ref</b> parameter in the <b>INCAR</b>. Note that this only works when ionic screening is present in the electrolyte (<b>C_MOLAR</b> $\gt$ 0); otherwise, the Fermi energy is undefined with respect to the vacuum level.

* <b>EFERMI_ref</b> = 0 (default) \
Electron chemical potential with respect to vacuum. Runs a constant potential calculation when <b>EFERMI_ref</b> $\lt$ 0
* <b>capacitance_init</b> = 1.0 (default) \
Initial guess for the capacitance of the unit cell (e/V), used for updating the number of electrons

### Recommended <b>INCAR</b> for an aqueous electrolyte

The default parameters correspond to pure water at 298 K. To model an electrolyte, it is necessary to specify the concentration (<b>C_MOLAR</b>) and the ionic radius (<b>R_ION</b>). If the latter is not specified, it defaults to the solvent radius which is likely too small. <b>EFERMI_ref</b> should be specified to run constant potential calculations. When using the default values for an aqueous electrolyte, we recommend setting <b>EFERMI_ref</b> $= -4.47 - U$, where $U$ is the electrode potential with respect to the standard hydrogen electrode.

```
LSOL = .TRUE.
ISOL = 2
C_MOLAR = 1.0        # set to the electrolyte concentration in mol/L
R_ION = 4.0          # set to the ionic radius in Angstrom
EFERMI_ref = -4.47   # set to the electron chemical potential in eV
```


## Output

The most important output quantities from VASPsol++ are the free energy and the Fermi level. The free energy printed by the main VASP program contains all necessary solvation corrections and can be used directly for computing free energy differences between states. Importantly, the free energy and Fermi level are properly referenced even for charged systems as long as ionic screening is present in the electrolyte. Unlike the original VASPsol implementation, there are no additional corrections that need to be applied these quantities.

In the case of constant potential calculations (<b>EFERMI_ref</b> $\lt$ 0), VASP prints the grand canonical potential $\Omega$ rather than the free energy $F$. This is defined as,

$$ \Omega = F - q_\mathrm{sol} \mu_\mathrm{e} $$

where $q_\mathrm{sol}$ is the solute charge using an <b><i> electron is positive </i></b> convention and $\mu_\mathrm{e}$ is the electron chemical potential with respect to vacuum (equal to the Fermi level printed in VASP).

All solvation corrections, including the additional $q_\mathrm{sol} \mu_\mathrm{e}$ term subtracted in constant potential calculations, are also applied to the $E_0$ value printed by VASP. This is the electronic energy extrapolated to 0 K and is the suggested quantity to use when computing free energy differences between states, even with solvation. The entropic contribution removed from $F$ to obtain $E_0$ arises from electronic entropy in the solute at high temperatures assossiated with the smearing procedure; these have nothing to do with entropy in the electrolye or any other non-electronic entropy.

Specifying <b>LVHAR = .TRUE.</b> or <b>LVTOT = .TRUE.</b> in the <b>INCAR</b> will cause the following additional files to be written at the end of the calculation. These files have the same format as <b>LOCPOT</b>.

Files written for both <b>ISOL</b>=1 and <b>ISOL</b>=2:

* <b>PHI</b> \
electrostatic potential, should be same as <b>LOCPOT</b>
* <b>PHI_SOLV</b> \
electrostatic potential from the solvent
* <b>VSOLV</b> \
cavity correction to the KS potential
* <b>RHOB</b> \
bound charge density
* <b>RHOION</b> \
electrolyte ionic charge density

Files written only for <b>ISOL</b>=1:

* <b>S</b> \
solvent cavity

Files written only for <b>ISOL</b>=2:

* <b>ELOC</b> \
electrostatic field in the z direction
* <b>P</b> \
solvent polarization density in the z direction
* <b>SVDW</b> \
vdW cavity
* <b>SSOLV</b> \
solvent cavity
* <b>SION</b> \
ionic cavity
* <b>SDIEL</b> \
dielectric cavity
* <b>SCAV</b> \
cavity used for calculating the cavity formation free energy, equal to <b>SSOLV</b> by default unless <b>R_CAV</b> is specified

