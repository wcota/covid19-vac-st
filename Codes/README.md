
**Description of main parameters:**

	1. *R0*: R0 value
	2. *tx*: value of vaccination rate x 1e5
	3. *dia*: delay in starting vaccination
	4. *ig*: Strategies used
		-*ig2*: Last age group (75+) => (70-74) => (65-69) => ... => first age group (0-4)				!DAP
		-*ig3*: Priority only for elders, (75+) => ... => (60-64) => all adults (20-59) => all children (0-19)		!HVP
		-*ig7*: All at same time											!NP
		-*ig8*: Descending order of the total number of contacts for each age group					!DCP

	5. *fxi* or *T*: results for age group "i" or for the total population, respectively
	6. *cenario*: contact patterns scenario
		-*cenario* = 0: Unmitigated, 100% of all contacts 						
		-*cenario* = 2: Social distancing, considers the reduction of contact in the social distancing scenario 

	7. *vacina*: which efficacies will be used (this variable is used for when we have more than one type of vaccine, eg vaccine = 0 => coronavac data)

**Description of files:**

	1. *contact_all.dat*: Unmitigated, 100% of all contacts 	
	2. *contact_all-sd.dat*: Social distancing, considers the reduction of contact in the social distancing scenario 
	3. *vaccines-data.dat*: vaccines efficacies

**Parameter adjustment:**

	There are two main types of .f90 programs, one generates the time evolution data, the other generates the data for the heatmap.
	For the .f90 program you need to enter the analysis parameters in the "initial conditions" within the program.

**Results:**

	For the temporal evolution program, the results will have the value of each parameter in their name:

		1. *R0*: R0 value
		2. *tx*: value of vaccination rate x 1e5
		5. *fxi* or *T*: results for age group "i" or for the total population, respectively
		4. *ig*: Strategie used
		3. *dia*: delay in starting vaccination
		6. *cenario*: contact patterns scenario
		7. *CV*: indicates which vaccine was used for efficacy (in this case "CV" refers for coronavac), you can change the acronym to indicate any other data you have available

	An exemple:

		*R0-13-tx-15-fx4-ig2-dia60-CV-cenario-2.dat*

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

		*R0-tx-T-ig8-dia30-CV-cenario-0.dat*

		Here we have the result for the heatmap of vaccination rate vs R0 (represented by "R0-tx") with the parameters:
	
		*T* = The results are for total population
		*ig8* = vaccination strategy 8 (DCP)
		*dia30* = 30 days delay
		*CV*: indicates which vaccine was used for efficacy (in this case "CV" refers for coronavac), you can change the acronym to indicate any other data you have available
		*cenario-0* = Scenario 0 contact patterns (Unmitigated)


		OBS: Within each result (.dat) is indicated in the first line what each column represents.
		
**Compilation**

		It is necessary to have in the same program directory (.f90) the "Data" folder, which contains the contact matrices of the country of interest and the data for the vaccine's efficacy.

		Fortran(Gfortran):

			*gfortran program.f90 -o exec
			*./exec

		Fortran(Intel):

			*ifort program.f90 -o exec
			*./exec





