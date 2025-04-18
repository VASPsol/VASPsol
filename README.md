# VASPsol++: A framework for implementing complex continuum fluid models in VASP density functional theory calculations

VASPsol++ is a framework for implementing complex continuum fluid models within density functional theory calculations performed using the [Vienna Ab initio Simulation Package](https://www.vasp.at/) (VASP). It is being actively developed within the [Plaisance group](https://www.lsu.edu/eng/che/people/faculty/plaisance.php) at Louisiana State University and has its origins in the [VASPsol](https://github.com/henniggroup/VASPsol) code developed in the [Hennig group](https://hennig.mse.ufl.edu/) at the University of Florida.

VASPsol++ adds a [nonlocal and nonlinear implicit electrolyte model](#references) to the [linear polarizable continuum model](#references) contained in the original VASPsol code. Additionally, it is written in a modular format that allows for easy addition of any new continuum solvation models that are developed in the future.

Although there are some significant differences in the implementation, the general ideas for the [nonlocal and nonlinear models]((#references)) were developed in the group of [Tomás Arias](https://www.lassp.cornell.edu/people/tomas-arias) at Cornell University.

## Key features of the nonlinear+nonlocal model

* Nonlinear dielectric and ionic responses to model the high electric fields present at charged electrodes
* Uses a nonlocal cavity definition to prevent unphysical electrolyte leakage into small regions and to allow for separate dielectric and ionic cavities
* Efficient and robust solver for the nonlinear Poisson-Boltzmann equation based on Newton's method with a line search
* Option to performing calculations at constant potential rather than constant number of electrons (grand canonical DFT)
* Only slightly higher computational cost than the linear VASPsol model


## Installation

VASPsol++ is mostly implemented in a single fortran file, solvation.F, that can be found in the src/ directory of the repository. Additionally, a patch file is required to make modifications to some of the original VASP source files. These patch files are version specific and must be obtained by emailing [Craig Plaisance](mailto:plaisance@lsu.edu). Currently, the following patch files are available:
* VASP 5.4.4 (with and without VTST)
* VASP 6.3.2 (with and without VTST)

The instructions below assume you have downloaded the VASP source files in a directory `<VASP_SRC>` and configured the necessary makefiles specified in the VASP installation instructions.

1. Download (using `wget`) one of the source file bundles from underneath 'Assets' on the page for the most recent release. Then extract it. For the `tar.gz` file, this is done with the command `tar -xzf vaspsol-pp-<version>.tar.gz`. We refer to the resulting source file directory as `<VASPSOL_SRC>`.

2. Copy `<VASPSOL_SRC>/src/solvation.F` to `<VASP_SRC>/src/`, replacing the skeleton `solvation.F` file already there.

3. Copy the appropriate patch file from `<VASPSOL_SRC>/src/patches` to `<VASP_SRC>/`. There are actually two patch files for each version of VASP. The first, named `vaspsol++-vtst-vasp_<version>.patch`, should be used if you are using the [VASP Transition State Tools](https://theory.cm.utexas.edu/vtsttools/) add-on developed by the [Henkelman group](http://henkelmanlab.org/) at UT Austin. Otherwise, use the patch file named `vaspsol++-vasp_<version>.patch`.

4. Apply the patch by running `patch -p1 < <patch_file>` in `<VASP_SRC>`.

5. Compile VASP as you normally would. We recommend using the Intel Fortran compiler since this is the only compiler we have tested.

## Input

### General

Parameters available for all solvation models:

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
* <b>EFERMI_tol</b> = 0 (default) \
Convergence criteria for the Fermi level. If <b>EFERMI_tol</b> $\le$ 0, it is set to 10&times;<b>EDIFF</b>
* <b>capacitance_init</b> = 1.0 (default) \
Initial guess for the capacitance of the unit cell (e/V), used for updating the number of electrons

### Recommended <b>INCAR</b> for an aqueous electrolyte

The default parameters correspond to pure water at 298 K. To model an electrolyte, it is necessary to specify the concentration (<b>C_MOLAR</b>) and the ionic radius (<b>R_ION</b>). If the latter is not specified, it defaults to the solvent radius which is likely too small. <b>EFERMI_ref</b> should be specified to run constant potential calculations. In principle, the value of the Fermi level corresponding to the standard hydrogen electrode ($\varepsilon_{\rm F, SHE}$) should be determined self consistently for the DFT functional and solvation parameters you are using. The procedure for this is given in our 2023 JCP manuscript. One would then set <b>EFERMI_ref</b> to $\varepsilon_{\rm F} = \varepsilon_{\rm F, SHE} - U_{\rm SHE}$, where $U_{\rm SHE}$ is the electrode potential with respect to the standard hydrogen electrode. When using the default values for an aqueous electrolyte with the BEEF-vdW functional, we obtain a value of $\varepsilon_{\rm F, SHE} = -4.57 \rm eV$. (<b>Note:</b> the value of $-4.47 \rm eV$ reported in our manuscript is incorrect)

```
LSOL = .TRUE.
ISOL = 2
C_MOLAR = 1.0        # set to the electrolyte concentration in mol/L
R_ION = 4.0          # set to the ionic radius in Angstrom
EFERMI_ref = -4.57   # set to the electron chemical potential in eV (read how to determine this value above)
```


## Output

The most important output quantities from VASPsol++ are the free energy and the Fermi level. The free energy printed by the main VASP program contains all necessary solvation corrections and can be used directly for computing free energy differences between states. Importantly, the free energy and Fermi level are properly referenced even for charged systems as long as ionic screening is present in the electrolyte. Unlike the original VASPsol implementation, there are no additional corrections that need to be applied these quantities.

In the case of constant potential calculations (<b>EFERMI_ref</b> $\lt$ 0), VASP prints the grand canonical potential $\Omega$ rather than the free energy $F$. This is defined as,

$$ \Omega = F - q_\mathrm{sol} \mu_\mathrm{e} $$

where $q_\mathrm{sol}$ is the solute charge using an <b><i> electron is positive </i></b> convention and $\mu_\mathrm{e}$ is the electron chemical potential with respect to vacuum (equal to the Fermi level printed in VASP).

All solvation corrections, including the additional $q_\mathrm{sol} \mu_\mathrm{e}$ term subtracted in constant potential calculations, are also applied to the $E_0$ value printed by VASP. This is the electronic energy extrapolated to 0 K and is the suggested quantity to use when computing free energy differences between states, even with solvation. The entropic contribution removed from $F$ to obtain $E_0$ arises from electronic entropy in the solute at high temperatures assossiated with the smearing procedure; these have nothing to do with entropy in the electrolye or any other non-electronic entropy.

Specifying <b>LVHAR = .TRUE.</b> or <b>LVTOT = .TRUE.</b> in the <b>INCAR</b> will cause the following additional files to be written at the end of the calculation. These files have the same format as <b>LOCPOT</b>.

Files written for both <b>ISOL</b>=1 and <b>ISOL</b>=2:

* <b>PHI</b> : electrostatic potential, should be same as <b>LOCPOT</b>
* <b>PHI_SOLV</b> : electrostatic potential from the solvent
* <b>VSOLV</b> : cavity correction to the KS potential
* <b>RHOB</b> : bound charge density
* <b>RHOION</b> : electrolyte ionic charge density

Files written only for <b>ISOL</b>=1:

* <b>S</b> : solvent cavity

Files written only for <b>ISOL</b>=2:

* <b>ELOC</b> : electrostatic field in the z direction
* <b>P</b> : solvent polarization density in the z direction
* <b>SVDW</b> : vdW cavity
* <b>SSOLV</b> : solvent cavity
* <b>SION</b> : ionic cavity
* <b>SDIEL</b> : dielectric cavity
* <b>SCAV</b> : cavity used for calculating the cavity formation free energy, equal to <b>SSOLV</b> by default unless <b>R_CAV</b> is specified


## References

<b>Nonlinear+nonlocal model</b>
* An implicit electrolyte model for plane wave density functional theory exhibiting nonlinear response and a nonlocal cavity definition. S.M.R. Islam, F. Khezeli, S. Ringe, and C. Plaisance, J. Chem. Phys. 159, 234117 (2023), (https://doi.org/10.1063/5.0176308).

<b>Original linear+local VASPsol method</b>
* Implicit solvation model for density-functional study of nanocrystal surfaces and reaction pathways. K. Mathew, R. Sundararaman, K. Letchworth-Weaver, T.A. Arias, and R.G. Hennig, J. Chem. Phys. 140, 084106 (2014), (https://doi.org/10.1063/1.4865107).
* Implicit self-consistent electrolyte model in plane-wave density-functional theory. K. Mathew, V.S. C. Kolluru, S. Mula, S.N. Steinmann, and R.G. Hennig, J. Chem. Phys. 151, 234101 (2019), (https://doi.org/10.1063/1.5132354).

<b>Original nonlinear+nonlocal models</b>
* The importance of nonlinear fluid response in joint density-functional theory studies of battery systems. D. Gunceler, K. Letchworth-Weaver, R. Sundararaman, K.A. Schwarz, T.A. Arias, Model. Simul. Mater. Sc. 21, 074005 (2013), (https://doi.org/10.1088/0965-0393/21/7/074005).
* Spicing up continuum solvation models with SaLSA: The spherically averaged liquid susceptibility ansatz. R. Sundararaman, K.A. Schwarz, K. Letchworth-Weaver, T.A. Arias, J. Chem. Phys. 142, 054102 (2015), (https://doi.org/10.1063/1.4906828).
