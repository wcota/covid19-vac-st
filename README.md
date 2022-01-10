
# Code for *Effects of infection fatality ratio and social contact matrices on vaccine prioritization strategies*

The codes for simulations were written in Fortran and compiled with the Intel Fortran Compiler. Data analysis and figures were done Python 3.10 and the following open source libraries: [pandas](https://pandas.pydata.org/), [matplotlib](https://matplotlib.org/), and [seaborn](https://seaborn.pydata.org/).

In this repository we show only the codes for the simulations.

## Folder structure

All codes and data are available in the folder `codes`. Each `*.f90` inside `codes/main` and `codes/supp` folders correspond to the respective figure in the main paper and supplementary material, respectively. The data (contact matrices and vaccine efficacies) are available at a common folder `codes/data`, with files:

- `vaccines-data.dat`: file with four columns - reduction in deaths with Vaxzevria, CoronaVac, and of infections with Vaxzevria and CoronaVac, respectively.
- `<country>/contact_<contact>.dat`: contact matrix C<sub>ij</sub> (16x16) for each `<country>` (Brazil, Germany or Uganda) and contact scenarios `<cenario>` - `all`  (unmitigated) or `all-sd` (social distancing).

## How to run

The codes can be compiled with the [Intel Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) or with [GNU Fortran compiler](https://gcc.gnu.org/wiki/GFortran). Examples:

- Intel: `ifort program.f90 -o program`
- GFortran: `gfortran program.f90 -o program`

To execute, use `./program > output.dat` or simply `./program`.

Some of the `*.f90` codes generates the time evolution data, and others the heatmaps. In each `*.f90` file, enter the corresponding parameters after the `!#Initial conditions` comment.

## Description of main parameters

The parameters used in the codes are the following:

1. `R0`: infection rate parameter ϖ
2. `csipop`: value of vaccination rate ξ
3. `dia`: days to start the vaccination (t<sub>v</sub>)
4. `ig`: One of the following prioritization strategies:
	- `2`: DAP strategy
	- `3`: HVP strategy
	- `7`: NP strategy
	- `8`: DCP strategy
5. `cenario`: contact patterns scenario
	- `0`: Unmitigated, 100% of all contacts 						
	- `2`: Social distancing, considers the reduction of contact in the social distancing scenario 
6. `vacina`: which vaccine data to use (see file `codes/data/vaccines-data.dat`)
	- `0`: CoronaVac
	- `1`: Vaxzevria 

## Output files

For the temporal evolution program, the results will have the value of each parameter in their name, as mentioned in the description of the main parameters. For the vaccine efficacies, `CV` refers to CoronaVac, but can be changed in the code for the corresponding name to the data used.

For the programs that generate the heatmap the name can have `R0-tx`, corresponding to the heatmap of vaccination rate ξ versus ϖ, or `dia-tx`, corresponding to the heatmap of delay t<sub>v</sub> versus ϖ.

The first line of each `*.dat` file indicates what the values of the columns mean.

### Examples:

#### `R0-13-tx-15-fx4-ig2-dia60-CV-cenario-2.dat`

Here we have the result for the temporal evolution with the parameters:

- `R0-13`: `R0 = 1.3`
- `tx-15`: `csipop = 0.0015`
- `fx4`: results are for age group 4
- `ig2`: `ig = 2` (DAP strategy)
- `dia60`: `dia = 60` (t<sub>v</sub> = 60)
- `cenario-2`: `cenario = 2` (social distancing scenario)
- `CV`: `vacina = 0` (CoronaVac)

#### `R0-tx-T-ig8-dia30-CV-cenario-0.dat`

Here we have the results for the heatmap of vaccination rate versus ϖ (represented by "R0-tx") with the parameters:

- `T`: the results are for total population
- `ig8`: `ig = 8` (DCP strategy)
- `dia30`: `dia = 30` (t<sub>v</sub> = 30)
- `CV`: `vacina = 0` (CoronaVac)
- `cenario-0`: `cenario = 0` (unmitigated scenario)

