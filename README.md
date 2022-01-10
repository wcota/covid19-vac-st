
# Code for *Effects of infection fatality ratio and social contact matrices on vaccine prioritization strategies*

The codes for simulations were written in Fortran and compiled with the Intel Fortran Compiler. Data analysis and figures were done Python 3.10 and the following open source libraries: [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), and [seaborn](https://seaborn.pydata.org/).

In this repository we show only the codes for the simulations.

## Folder structure

All codes and data are available in the folder `codes`. Each `*.f90` inside `codes/main` and `codes/supp` folders correspond to the respective figure in the main paper and supplementary material, respectively. The data (contact matrices and vaccine efficacies) are available at a common folder `codes/data`, with files:

- `vaccines-data.dat`: file with four columns - reduction in deaths with Vaxzevria, CoronaVac, and of infections with Vaxzevria and CoronaVac, respectively.
- `<country>/contact_<contact>.dat`: contact matrix $C_{ij}$ (16x16) for each `<country>` and contact scenarios `all`  (unmitigated) and `all-sd` (social distancing).

## How to run

The codes can be compiled with the [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) or with [GNU Fortran compiler](https://gcc.gnu.org/wiki/GFortran). Examples:

- Intel: `ifort program.f90 -o program`
- GFortran: `gfortran program.f90 -o program`

To execute, use `./program > output.dat` or simply `./program`.

Some of the `*.f90` codes generates the time evolution data, and others the heatmaps. In each `*.f90` file, enter the corresponding parameters after the `!#Initial conditions` comment.

## Description of main parameters

1. *R0*: R0 value
2. *tx*: value of vaccination rate x 1e5
3. *dia*: delay in starting vaccination
4. *ig*: Strategies used
	- *ig2*: Last age group (75+) => (70-74) => (65-69) => ... => first age group (0-4)				!DAP
	- *ig3*: Priority only for elders, (75+) => ... => (60-64) => all adults (20-59) => all children (0-19)		!HVP
	- *ig7*: All at same time											!NP
	- *ig8*: Descending order of the total number of contacts for each age group					!DCP
5. *fxi* or *T*: results for age group "i" or for the total population, respectively
6. *cenario*: contact patterns scenario
	- *cenario* = 0: Unmitigated, 100% of all contacts 						
	- *cenario* = 2: Social distancing, considers the reduction of contact in the social distancing scenario 
7. *vacina*: which efficacies will be used (this variable is used for when we have more than one type of vaccine, eg vaccine = 0 => coronavac data)

## Output files

For the temporal evolution program, the results will have the value of each parameter in their name:

1. *R0*: R0 value
2. *tx*: value of vaccination rate x 1e5
3. *fxi* or *T*: results for age group "i" or for the total population, respectively
4. *ig*: Strategie used
5. *dia*: delay in starting vaccination
6. *cenario*: contact patterns scenario
7. *CV*: indicates which vaccine was used for efficacy (in this case "CV" refers for coronavac), you can change the acronym to indicate any other data you have available

Examples:

### `R0-13-tx-15-fx4-ig2-dia60-CV-cenario-2.dat`

Here we have the result for the temporal evolution with the parameters:

*R0* = 1.3
*tx-15* = Total vaccination rate = 0.0015
*fx4* = the results are for age group 4
*ig2* = vaccination strategy 2 (DAP)
*dia60* = 60 days delay
*cenario-2* = Scenario 2 contact patterns (Social distancing)
*CV*: indicates which vaccine was used for efficacy (in this case "CV" refers for coronavac), you can change the acronym to indicate any other data you have available

For the programs that generate the heatmap the name of the results files will have the following meanings:

*R0-tx*: heatmap for vaccination rate vs R0
*dia-tx*: heatmap for vaccination rate vs delay in starting vaccination
*fxi* or *T*: resultados para o grupo etário "i" ou para a população total, respectivamente
*ig*: Strategie used
*dia*: delay in starting vaccination

An exemple:

### `R0-tx-T-ig8-dia30-CV-cenario-0.dat`

Here we have the result for the heatmap of vaccination rate vs R0 (represented by "R0-tx") with the parameters:

*T* = The results are for total population
*ig8* = vaccination strategy 8 (DCP)
*dia30* = 30 days delay
*CV*: indicates which vaccine was used for efficacy (in this case "CV" refers for coronavac), you can change the acronym to indicate any other data you have available
*cenario-0* = Scenario 0 contact patterns (Unmitigated)


OBS: Within each result (.dat) is indicated in the first line what each column represents.
	






