program sirs
implicit none

real*8 :: t0, tN, h, ti, dt, poptotal, R0, soma, kmed, csipop, gama, ptime,imtime, auxreal
integer :: HH, ii, iii, j, jj, amost, namost, ig, dia, cenario, l, faixa, ifaixa, vacina

real*8, allocatable :: lambda(:,:), csi(:,:), mu(:,:), beta(:,:), alfa(:,:), nu(:,:), diafx(:),kmedioaux(:,:),alfadp(:),ef(:),&
&tau(:,:), N(:,:), Naux(:,:), k(:,:,:), kk(:,:), kmedio(:), Nmedio(:,:,:), peso(:), populacao(:), theta(:), kaux(:),reducao(:)

character*199 :: string, string1, string2, string3, string4, string5, string6

!#######################################
!parameters

t0 = 0d0		!Beginning of the epidemic
tN = 500d0		!Maximum time considered for analyzing the results
HH = 1000		!Number of intervals for numerical integration
h = (tN - t0)/(HH*1d0)
poptotal = 100000d0	!Total population considered (what interferes the results is the proportion of the population for each age group)
namost = 16		!Number of age groups considered

!#######################################
!Initial conditions

allocate(lambda(2,16), csi(5,16), mu(2,16), beta(4,16), alfa(4,16), nu(3,16), tau(2,16), N(17,16),diafx(16), kaux(2), ef(16), &
&Naux(17,16), kk(16,16), k(17,16,4), kmedio(16), Nmedio(HH,17,16),peso(16), populacao(16), theta(16),kmedioaux(2,16), alfadp(16),&
&reducao(16))

    !Indices for each age group
    !1:  0-4
    !2:  5-9
    !3:  10-14
    !4:  15-19
    !5:  20-24
    !6:  25-29
    !7:  30-34
    !8:  35-39
    !9:  40-44
    !10: 45-49
    !11: 50-54
    !12: 55-59
    !13: 60-64
    !14: 65-69
    !15: 70-74
    !16: 75+

 dia = 30
 ig = 2
 cenario = 2
 R0 = 1.30d0
 csipop = 0.00150d0 	
 vacina = 0

!dia = delay in starting vaccination
!ig = vaccination strategy
!R0 = R0 value calculated from contacts for social distancing
!cenario = contact patterns scenario
!vacina = which efficacies will be used (this variable is used for when we have more than one type of vaccine, eg vaccine = 0 => coronavac data)
!csipop = total vaccination rate

! cenario = 0 all			!100% of all contacts 							!used in results
! cenario = 1 all without school	!excludes contacts only from schools 					!not used in results
! cenario = 2 all with sd		!considers the reduction of contact in the social distancing scenario 	!used in the results
! cenario = 3 lockdown			!considers contact reduction in lockdown scenario			!not used in results

!vacina = 0 (coronavac), vacina = 1 (pfizer/astrazeneca)

!ig's:
!1: All elders, (60+) => all adults (20-59) => all children (0-19)
!2: Last age group (75+) => (70-74) => (65-69) => ... => first age group (0-4)						!DAP
!3: Priority only for elders, (75+) => ... => (60-64) => all adults (20-59) => all children (0-19)			!HVP
!4: All adults (20-59) => all elders (60+)
!5: First age group of adults and rising (20-24) => (25-29) => ... (75+)
!6: Same as strategy 5, but starting with children (0-4)
!7: All at same time													!NP
!8: Descending order of the total number of contacts for each age group							!DCP

ef = 0			!where the data will be stored for effectiveness against infection
reducao = 0		!where data will be stored for effectiveness against deaths

open(1,file="../data/vaccines-data.dat")			!reads the efficacy data of the vaccine considered

read(1,*)

if (vacina .eq. 0) then   !auxreal has no meaning, it's just to read the column we don't want and jump to the next value of interest

    ptime = 21d0 !1/nu_p		
    imtime = 7d0 			!CORONAVAC

    do ii = 1,16

        read(1,*) auxreal, reducao(ii), auxreal, ef(ii)

    end do

else if (vacina .eq. 1) then

    ptime = 21d0 !tempo de proteção
    imtime = 35d0 !tempo de imunização			!AZ/PFIZER

    do ii = 1,16

        read(1,*) reducao(ii), auxreal, ef(ii), auxreal

    end do

end if

reducao = reducao/100d0
ef = ef/100d0

close(1)


faixa = 0

if (ig .eq. 2 .or. ig .eq. 3) then

    faixa = 15

else if (ig .eq. 5) then

    faixa = 6

else if (ig .eq. 6 .or. ig .eq. 8) then

    faixa = 2

end if

 Nmedio = 0d0

if (cenario .eq. 0) then  			!these if's will read the contact patterns equivalent to the scenario (and country) used in the simulation

 kmedio = [9.06104972622234,&
&18.53416619469897,&
&24.944777393540704,&
&28.83814524492348,&
&18.44064843675028,&
&17.7443525026541,&
&16.764214270277833,&
&15.833199279359286,&
&15.038997939348693,&
&12.749506122497436,&
&11.548589867882228,&
&8.908087617136696,&
&5.766233304795698,&
&4.488684584905327,&
&3.6386329582863897,&
&2.926288313680333]

    !com escola


 open(33,file="../data/Brazil/contact_all.dat")

    do ii = 1,16

        read(33,*) kk(ii,:)

    end do

 close(33)

else if (cenario .eq. 1) then

 kmedio = [8.49786590562505,&
&11.775275350653459,&
&16.102568319432287,&
&18.66140733649442,&
&17.197672679076398,&
&16.909074320839306,&
&15.858722276247798,&
&15.142699108411534,&
&14.315104974852478,&
&11.854082354191457,&
&10.609107377404492,&
&8.213165399250961,&
&5.504699508670181,&
&4.392636558980293,&
&3.5519870345904874,&
&2.892123631908877]

 kk = 0d0 

 open(33,file="../data/Brazil/contact_all-w_school.dat")

    do ii = 1,16

        read(33,*) kk(ii,:)

    end do

 close(33)

else if (cenario .eq. 2) then

 kmedio = [6.066509213452985,&
&10.930551348598055,&
&13.616615299414576,&
&14.627455323011965,&
&9.133033662101811,&
&9.051618251406843,&
&8.853469601672689,&
&8.70091453873895,&
&8.370970709656403,&
&7.300762923669087,&
&6.581118313232487,&
&5.417379902686264,&
&3.7341551974408596,&
&3.126298556926259,&
&2.44983480267834,&
&2.418587539545057]

 kk = 0d0 

 open(33,file="../data/Brazil/contact_all-sd.dat")

    do ii = 1,16

        read(33,*) kk(ii,:)

    end do

 close(33)

else if (cenario .eq. 3) then

 kmedio = [5.784915703083305,&
&7.535414253794091,&
&9.083154884909021,&
&8.940805526782992,&
&7.635891366785771,&
&7.439945416740162,&
&7.207761595740647,&
&7.118569562043913,&
&6.731608109797762,&
&5.790894946576102,&
&5.206235849688046,&
&4.559399532909547,&
&3.433405274413712,&
&3.0781809674286924,&
&2.40644339764307,&
&2.4014630350224184]

 kk = 0d0 

 open(33,file="../data/Brazil/contact_lockdown.dat")

    do ii = 1,16

        read(33,*) kk(ii,:)

    end do

 close(33)

end if



do ii = 1,16				!the next two loops will be used in strategy 8

    kmedioaux(1,ii) = ii
    kmedioaux(2,ii) = kmedio(ii)

end do

do ii = 1,15

    do jj = ii+1, 16

        if (kmedioaux(2,ii) .lt. kmedioaux(2,jj)) then

            kaux(:) = kmedioaux(:,ii)

            kmedioaux(:,ii) = kmedioaux(:,jj)

            kmedioaux(:,jj) = kaux(:)

        end if

    end do

end do

 gama = 0.5d0 			!factor reduction in the transmissibility of infected vaccinated (IP)

 theta = 0d0			!theta(i) = Infection fatality rate for the age group i

	theta(1) = 0.0000161d0 !0-4
	theta(2) = 0.0000161d0 !5-9
	theta(3) = 0.0000695d0 !10-14
	theta(4) = 0.0000695d0 !15-19
	theta(5) = 0.000309d0  !20-24
	theta(6) = 0.000309d0  !25-29
	theta(7) = 0.000844d0  !30-34 
	theta(8) = 0.000844d0  !35-39
	theta(9) = 0.00161d0   !40-44
	theta(10) = 0.00161d0  !45-49
	theta(11) = 0.00595d0  !50-54
	theta(12) = 0.00595d0  !55-59
	theta(13) = 0.01930d0  !60-64
	theta(14) = 0.01930d0  !65-69
	theta(15) = 0.04280d0  !70-74
	theta(16) = 0.07800d0  !75+

  k = 0d0

 mu = 0d0

	mu = 1d0/5.2d0

	!mu(1,i) = mu A
	!mu(2,i) = mu Av

 beta = 0d0

	beta(1,:) = 1d0/5.8d0
	beta(2,:) = 1d0/2.6d0
	beta(3,:) = 1d0/2.6d0
	beta(4,:) = 1d0/5.8d0

	!beta(1,i) = beta R
	!beta(2,i) = beta I
	!beta(3,i) = beta Iv
	!beta(4,i) = beta Rv

 alfa = 0d0

	alfa(1,:) = 1d0/3.2d0

	do ii = 1,16

		alfa(2,ii) = theta(ii)/(0.69d0-theta(ii))*alfa(1,ii)

	end do

	alfa(4,:) = alfa(1,:)

	alfa(3,:) = alfa(2,:)

	!alfa(1,i) = alfa R
	!alfa(2,i) = alfa D
	!alfa(3,i) = alfa Dv
	!alfa(4,i) = alfa Rv

 alfadp = 0d0  		!alfa_Dp

 do ii = 1,16

    if (reducao(ii) .eq. 1d0) then
    
        alfadp(ii) = 0d0

    else

        alfadp(ii) =  ((1d0-reducao(ii))*theta(ii)*beta(4,ii))/(1d0-(1d0-reducao(ii))*theta(ii))  

    end if

 end do

 nu = 0d0

    do ii = 1,16

	nu(1,ii) = 1d0/ptime
	nu(2,ii) = (1d0-ef(ii))*1d0/imtime
    	nu(3,ii) = ef(ii)*1d0/imtime

    end do

    	!nu(1,i) = nu P
	!nu(2,i) = nu Sp
	!nu(3,i) = nu Rp

 N = 0d0
 populacao = 0d0

    populacao = [0.0689267 , 0.0688616 , 0.06866363, 0.07279936, 0.07996429,&
    &0.07974595, 0.08032925, 0.08027336, 0.0746512 , 0.06530365,&
    &0.05977749, 0.05383225, 0.04534166, 0.03573888, 0.02645304,&
    &0.0393377 ]

    populacao = populacao*poptotal

	!N(1,i) = S
	!N(2,i) = E
	!N(3,i) = A
	!N(4,i) = I
	!N(5,i) = D
	!N(6,i) = R
	!N(7,i) = Sv
	!N(8,i) = Sp
	!N(9,i) = Ip
	!N(10,i) = Rp*		!!! Before, this compartment was for the fully immunized (coming from P)
	!N(11,i) = Ev
	!N(12,i) = Av
	!N(13,i) = Iv
	!N(14,i) = Rv
	!N(15,i) = Dv
	!N(16,i) = Rp**		!!! Before, this compartment was for protected vaccinated who recovered (coming from IP)
   	!N(17,i) = P	

				!!! before we divided the immunized and recovered protected into two compartments (N(10,i) and N(16,i)), the Rp is the sum of both

	!The compartment for the dead of the group of protected (Dp) was not separated in this program, as it was counted together with the Dv compartment

 peso = 0d0

 do ii = 1,16

    peso(ii) = kmedio(ii)*populacao(ii)

 end do

 kmed = sum(peso)/poptotal

 peso = peso/sum(peso)

	 lambda = 0d0

		soma = 0d0

        if (cenario .eq. 2) go to 1

         kmedio = [6.066509213452985,&
        &10.930551348598055,&
        &13.616615299414576,&
        &14.627455323011965,&
        &9.133033662101811,&
        &9.051618251406843,&
        &8.853469601672689,&
        &8.70091453873895,&
        &8.370970709656403,&
        &7.300762923669087,&
        &6.581118313232487,&
        &5.417379902686264,&
        &3.7341551974408596,&
        &3.126298556926259,&
        &2.44983480267834,&
        &2.418587539545057]

         kk = 0d0

         open(33,file="../data/Brazil/contact_all-sd.dat")

         do ii = 1,16

                read(33,*) kk(ii,:)

         end do

         close(33)

         peso = 0d0

         do ii = 1,16

            peso(ii) = kmedio(ii)*populacao(ii)

         end do

         kmed = sum(peso)/poptotal

         peso = peso/sum(peso)

         1 continue

	do j = 1,16

		do jj = 1,16

			soma = soma + (populacao(j)/poptotal)*(populacao(jj)/poptotal)*kmedio(j)*kmedio(jj)&
			&*kk(j,jj)/(kmed*(beta(1,j)+beta(2,j)))*(1d0+beta(2,j)/(alfa(1,j)+alfa(2,j)))

		end do

	end do	

	!lambda = transmission rate per contact of infectious compartments (we consider the same for all)
	!lambda(1,:) = lambda(2,:), inheritance from an old model, was not omitted as it does not change the results as all lambdas are equal
	!note that the transmission rate per contact is the same for all compartments, with the exception of IP in which we include a reduction factor

	lambda = R0/soma	!the transmission rate is calculated from the omega parameter (R0 of social distance)

        tau = 0d0

	tau(1,:) = lambda(1,:)
	tau(2,:) = beta(1,:)


    	!tau(1,i) = lambda, the reduction factor for IP transmission is included in the equations
    	!tau(2,i) = alfa Rp

        if (cenario .eq. 0) then		!here we go back to correcting the contact patterns for the considered scenario

   	      kmedio = [9.06104972622234,&
        &18.53416619469897,&
        &24.944777393540704,&
        &28.83814524492348,&
        &18.44064843675028,&
        &17.7443525026541,&
        &16.764214270277833,&
        &15.833199279359286,&
        &15.038997939348693,&
        &12.749506122497436,&
        &11.548589867882228,&
        &8.908087617136696,&
        &5.766233304795698,&
        &4.488684584905327,&
        &3.6386329582863897,&
        &2.926288313680333]

         	open(33,file="../data/Brazil/contact_all.dat")

            	do ii = 1,16

                	read(33,*) kk(ii,:)

            	end do

         	close(33)

        else if (cenario .eq. 1) then

         	kmedio = [8.49786590562505,&
        &11.775275350653459,&
        &16.102568319432287,&
        &18.66140733649442,&
        &17.197672679076398,&
        &16.909074320839306,&
        &15.858722276247798,&
        &15.142699108411534,&
        &14.315104974852478,&
        &11.854082354191457,&
        &10.609107377404492,&
        &8.213165399250961,&
        &5.504699508670181,&
        &4.392636558980293,&
        &3.5519870345904874,&
        &2.892123631908877]

         	kk = 0d0

         	open(33,file="../data/Brazil/contact_all-w_school.dat")

            	do ii = 1,16

                	read(33,*) kk(ii,:)

            	end do

         	close(33)

        else if (cenario .eq. 2) then

         	kmedio = [6.066509213452985,&
        &10.930551348598055,&
        &13.616615299414576,&
        &14.627455323011965,&
        &9.133033662101811,&
        &9.051618251406843,&
        &8.853469601672689,&
        &8.70091453873895,&
        &8.370970709656403,&
        &7.300762923669087,&
        &6.581118313232487,&
        &5.417379902686264,&
        &3.7341551974408596,&
        &3.126298556926259,&
        &2.44983480267834,&
        &2.418587539545057]

         	kk = 0d0 

         	open(33,file="../data/Brazil/contact_all-sd.dat")

            	do ii = 1,16

                	read(33,*) kk(ii,:)

            	end do

         	close(33)

        else if (cenario .eq. 3) then

         	kmedio = [5.784915703083305,&
        &7.535414253794091,&
        &9.083154884909021,&
        &8.940805526782992,&
        &7.635891366785771,&
        &7.439945416740162,&
        &7.207761595740647,&
        &7.118569562043913,&
        &6.731608109797762,&
        &5.790894946576102,&
        &5.206235849688046,&
        &4.559399532909547,&
        &3.433405274413712,&
        &3.0781809674286924,&
        &2.40644339764307,&
        &2.4014630350224184]

         	kk = 0d0 

         	open(33,file="../data/Brazil/contact_lockdown.dat")

            	do ii = 1,16

                	read(33,*) kk(ii,:)

            	end do

         	close(33)

        end if

     	peso = 0d0

     	do ii = 1,16

        	peso(ii) = kmedio(ii)*populacao(ii)

     	end do

     	kmed = sum(peso)/poptotal

     	peso = peso/sum(peso)

	write(string,*) int(ig)
	write(string1,*) dia
	write(string3,*) cenario
	write(string4,*) int(R0*10)
	write(string5,*) int(csipop*10000)
	if (vacina .eq. 0) then
	    string2 =  'CV'
	else
	    string2 =  'Az-Pf'
	end if

	! fx i = age group i
	! T = total population

    	open(10,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx1-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(20,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx2-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(30,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx3-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(40,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx4-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(50,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx5-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(60,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx6-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(70,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx7-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(80,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx8-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
   	 open(90,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx9-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(100,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx10-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(110,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx11-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
   	 open(120,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx12-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(130,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx13-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(140,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx14-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(150,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx15-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')
    	open(160,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-fx16-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')

    	open(170,file = 'R0-'//trim(adjustl(string4))//'-tx-&
&'//trim(adjustl(string5))//'-T-ig'//trim(adjustl(string))//'-dia'&
&//trim(adjustl(string1))//'-'//trim(adjustl(string2))//'-cenario-'&
&//trim(adjustl(string3))//'.dat') !write(i,'(5E20.8)')

   	!ordem t, 1, 2, 3, 4, 6, 5, 7, 11, 12, 13, 14, 15, 17, 8, 9, 16, 10

        do ifaixa = 1,16
  
    		write(ifaixa*10,*) "#t,S,E,A,I,R,D,S_V,E_V,A_V,I_V,R_V,D_V,P,S_P,I_P,R_P**,R_P* ## absolute values"

        end do

    	write(170,*) "#t,S,E,A,I,R,D,S_V,E_V,A_V,I_V,R_V,D_V,P,S_P,I_P,R_P**,R_P* ## absolute values"


	!#############
	!Dynamics

	Nmedio = 0d0

	do amost = 1, namost

            faixa = 0

            if (ig .eq. 2 .or. ig .eq. 3) then

                faixa = 15

            else if (ig .eq. 5) then

                faixa = 6

            else if (ig .eq. 6 .or. ig .eq. 8) then

                faixa = 2

            end if

			N = 0d0

            do iii = 1,16
    
                N(1,iii) = populacao(iii)

            end do

            N(2,amost) = 1d0				!initial condition, one exposed in the "amost" age group population
            N(1,amost) = N(1,amost) - N(2,amost)

	    ti = t0

            csi = 0d0

	    do ii = 1, HH

	          !csi(j,i) = vaccination rate j for each group i
	          !j = 1 = S
	          !j = 2 = E
	          !j = 3 = A
	          !j = 4 = I
	          !j = 5 = R

                  if (ii*h .ge. dia*1d0) then		!if time t has already reached the necessary delay to start vaccination

                    	if (ig .eq. 1) then

                        	if (faixa .eq. 2) go to 10
                        	if (faixa .eq. 3) go to 11

		                if ((sum(N(1:6,13:16)))/sum(populacao(13:16)) .ge. 0.2d0) then

!	                        	csi(:,1:12) = 0d0 

	                        	csi(:,13:16) = csipop*poptotal/(sum(N(1:3,13:16)) + sum(N(6,13:16)))

                            		go to 666

                        	end if

                        	10 continue

		                if (sum(N(1:6,5:12))/sum(populacao(5:12)) .ge. 0.2d0) then

!	                        	csi(:,1:4) = 0d0
                            		csi(:,5:16) = csipop*poptotal/(sum(N(1:3,5:16)) + sum(N(6,5:16)))
                            		faixa = 2

                            		go to 666

                        	end if
                
                        	11 continue

                        	if (sum(N(1:6,:))/poptotal .ge. 0.005) then

    	                    		csi = csipop*poptotal/(sum(N(1:3,:)) + sum(N(6,:)))
                            		faixa = 3   

                       		 else
                
                            		csi = 0            

		          	end if

                    else if (ig .eq. 2) then

                        	if (faixa .eq. 1) go to 201            

		                if ((sum(N(1:6,16)))/populacao(16) .ge. 0.2d0) then

!	                        	csi(:,1:15) = 0d0 

	                        	csi(:,16) = csipop*poptotal/(sum(N(1:3,16)) + (N(6,16)))

		                else if ((sum(N(1:6,faixa)))/(populacao(faixa)) .ge. 0.2d0 .and. faixa .gt. 1) then

!	                        	csi(:,1:faixa-1) = 0d0 

	                        	csi(:,faixa:16) = csipop*poptotal/(sum(N(1:3,faixa:16)) + sum(N(6,faixa:16)))

                        	else if (faixa .gt. 1) then

  !                          		csi(:,1:faixa-1) = 0d0

                            		faixa = faixa - 1

	                        	csi(:,faixa:16) = csipop*poptotal/(sum(N(1:3,faixa:16)) + sum(N(6,faixa:16)))

                        	end if

                        	go to 666

                        	201 continue
                
                        	if (sum(N(1:6,:))/poptotal .ge. 0.005) then

    	                    		csi = csipop*poptotal/(sum(N(1:3,:)) + sum(N(6,:)))


                        	else
                
                            		csi = 0           

		        	end if

                    else if (ig .eq. 3) then   

                        	if (faixa .eq. 12) go to 301
                        	if (faixa .eq. 5) go to 302                                   

		                if ((sum(N(1:6,16)))/populacao(16) .ge. 0.2d0) then

	!                        	csi(:,1:15) = 0d0 

	                        	csi(:,16) = csipop*poptotal/(sum(N(1:3,16)) + (N(6,16)))

		                else if ((sum(N(1:6,faixa)))/(populacao(faixa)) .ge. 0.2d0 .and. faixa .gt. 12) then

!	                        	csi(:,1:faixa-1) = 0d0 

	                        	csi(:,faixa:16) = csipop*poptotal/(sum(N(1:3,faixa:16)) + sum(N(6,faixa:16)))

                        	else if (faixa .gt. 13) then

   !                         		csi(:,1:faixa-1) = 0d0

                            		faixa = faixa - 1

	                        	csi(:,faixa:16) = csipop*poptotal/(sum(N(1:3,faixa:16)) + sum(N(6,faixa:16)))

                        	else if (faixa .eq. 13) then 

                            		faixa = 12
                            		go to 301

                        	end if

                        	go to 666

                        	301 continue

		                if (sum(N(1:6,5:12))/sum(populacao(5:12)) .ge. 0.2d0) then

!	                        	csi(:,1:4) = 0d0
                            		csi(:,5:16) = csipop*poptotal/(sum(N(1:3,5:16)) + sum(N(6,5:16)))

                            		go to 666

                        	end if

                        	302 continue
                
                        	if (sum(N(1:6,:))/poptotal .ge. 0.005) then

    	                    		csi = csipop*poptotal/(sum(N(1:3,:)) + sum(N(6,:)))
                            		faixa = 5  


                        	else

                            		csi = 0         

		                end if                    

                     	else if (ig .eq. 4) then

                        	if (faixa .eq. 2) go to 20
                        	if (faixa .eq. 3) go to 21

		                if (sum(N(1:6,5:12))/sum(populacao(5:12)) .ge. 0.2d0) then

	  !                      	csi(:,1:4) = 0d0

	                        	csi(:,5:12) = csipop*poptotal/(sum(N(1:3,5:12)) + sum(N(6,5:12)))

!	                        	csi(:,13:16) = 0d0

                            		go to 666

                        	end if

                        	20 continue

		                	if (sum(N(1:6,13:16))/sum(populacao(13:16)) .ge. 0.2d0) then

 !                           		csi(:,1:4) = 0d0
	                        	csi(:,5:16) = csipop*poptotal/(sum(N(1:3,5:16)) + sum(N(6,5:16)))
                    
                            		faixa = 2
                            		go to 666

                        	end if

                        	21 continue

                        	if (sum(N(1:6,5:16))/sum(populacao(5:16)) .ge. 0.005) then

    	                    		csi(:,5:16) = csipop*poptotal/(sum(N(1:3,5:16)) + sum(N(6,5:16)))

                            		faixa = 3 

                        	else 

                            		csi = 0                

		                end if 

                    else if (ig .eq. 5) then 

                        	if (faixa .eq. 16) go to 501            

		                if ((sum(N(1:6,5)))/populacao(5) .ge. 0.2d0) then

!	                        	csi(:,1:15) = 0d0 

	                        	csi(:,5) = csipop*poptotal/(sum(N(1:3,5)) + (N(6,5)))

		                else if ((sum(N(1:6,faixa)))/(populacao(faixa)) .ge. 0.2d0 .and. faixa .lt. 16) then

!	                        	csi(:,1:faixa-1) = 0d0 

	                        	csi(:,5:faixa) = csipop*poptotal/(sum(N(1:3,5:faixa)) + sum(N(6,5:faixa)))

                        	else if (faixa .lt. 16) then

  !                          		csi(:,1:faixa-1) = 0d0

                            		faixa = faixa + 1

	                        	csi(:,5:faixa) = csipop*poptotal/(sum(N(1:3,5:faixa)) + sum(N(6,5:faixa)))

                        	end if

                        	go to 666

                        	501 continue
                
                        	if (sum(N(1:6,5:16))/sum(populacao(5:16)) .ge. 0.005) then

    	                    		csi(:,5:16) = csipop*poptotal/(sum(N(1:3,5:16)) + sum(N(6,5:16))) 

                        	else 

                            		csi = 0                

		                end if 

                    else if (ig .eq. 6) then 

                        	if (faixa .eq. 16) go to 601            

		                if ((sum(N(1:6,1)))/populacao(1) .ge. 0.2d0) then

!	                        	csi(:,1:15) = 0d0 

	                        	csi(:,1) = csipop*poptotal/(sum(N(1:3,1)) + (N(6,1)))

		                	else if ((sum(N(1:6,faixa)))/(populacao(faixa)) .ge. 0.2d0 .and. faixa .lt. 16) then

!	                        	csi(:,1:faixa-1) = 0d0 

	                        	csi(:,1:faixa) = csipop*poptotal/(sum(N(1:3,1:faixa)) + sum(N(6,1:faixa)))

                        	else if (faixa .lt. 16) then

  !                          		csi(:,1:faixa-1) = 0d0

                            		faixa = faixa + 1

	                        	csi(:,1:faixa) = csipop*poptotal/(sum(N(1:3,1:faixa)) + sum(N(6,1:faixa)))

                        	end if

                        	go to 666

                        	601 continue
                
                        	if (sum(N(1:6,:))/poptotal .ge. 0.005) then

    	                    		csi = csipop*poptotal/(sum(N(1:3,:)) + sum(N(6,:)))

                        	else
                
                            		csi = 0           

		                end if  

                    else if (ig .eq. 7) then 

                        	if (sum(N(1:6,5:16))/poptotal .ge. 0.005) then

    	                    		csi = csipop*poptotal/(sum(N(1:3,:)) + sum(N(6,:)))

                        	else
                
                            		csi = 0           

		                end if  

                    else if (ig .eq. 8) then 

                       		if (faixa .eq. 16) go to 801            

		                if ((sum(N(1:6,int(kmedioaux(1,1)))))/populacao(int(kmedioaux(1,1))) .ge. 0.2d0) then

!	                        	csi(:,1:15) = 0d0 

	                       		csi(:,int(kmedioaux(1,1))) = csipop*poptotal/(sum(N(1:3,int(kmedioaux(1,1))))&
& + (N(6,int(kmedioaux(1,1)))))

		                else if ((sum(N(1:6,int(kmedioaux(1,faixa)))))/(populacao(int(kmedioaux(1,faixa))))&
& .ge. 0.2d0 .and. faixa .lt. 16) then

!	                        	csi(:,1:faixa-1) = 0d0 

                            		soma = 0d0

                            		do jj = 1,faixa

                                		soma = soma + (sum(N(1:3,int(kmedioaux(1,jj)))) + (N(6,int(kmedioaux(1,jj)))))

                           		end do

                            		do jj = 1,faixa

	                             		csi(:,int(kmedioaux(1,jj))) = csipop*poptotal/soma

                            		end do

                        	else if (faixa .lt. 16) then

  !                          		csi(:,1:faixa-1) = 0d0

                            		faixa = faixa + 1

                            		soma = 0d0

                            		do jj = 1,faixa

                                		soma = soma + (sum(N(1:3,int(kmedioaux(1,jj)))) + (N(6,int(kmedioaux(1,jj)))))

                            		end do

                            		do jj = 1,faixa

	                             		csi(:,int(kmedioaux(1,jj))) = csipop*poptotal/soma

                            		end do

                        	end if

                        	go to 666

                        	801 continue
                
                        	if (sum(N(1:6,:))/poptotal .ge. 0.005) then

    	                    		csi = csipop*poptotal/(sum(N(1:3,:)) + sum(N(6,:)))

                        	else
                
                            		csi = 0           

		                end if 

                    end if

                    666 continue

                    csi(4,:) = 0		  !the vaccination rate for the infected is always 0 (they don't get vaccinated)
    
                else
    
                    csi = 0d0

                end if

		Naux = N

		ti = ti + h

		do jj = 1,16

	        	soma = 0d0

			do l = 1,16		 !Runge-Kutta fourth order numerical integration

				soma = soma + kmedio(jj)*kk(jj,l)*(Naux(3,l) + Naux(4,l) + Naux(12,l) + Naux(13,l) + gama*Naux(9,l))/poptotal

			end do

			do j = 1,17

				k(j,jj,1) = h*dt(j,jj,Naux,lambda,csi,mu,beta,alfa,gama,nu,&
				&tau,kk,kmedio,poptotal,soma,alfadp)

			end do

		end do

		do jj = 1,16

			do j = 1,17

				Naux(j,jj) = N(j,jj) + k(j,jj,1)*0.5d0
					
			end do

		end do

		do jj = 1,16

	       		soma = 0d0

	 		do l = 1,16

		                soma = soma + kmedio(jj)*kk(jj,l)*(Naux(3,l) + Naux(4,l) + Naux(12,l) + Naux(13,l) + gama*Naux(9,l))/poptotal

	                end do

			do j = 1,17

				k(j,jj,2) = h*dt(j,jj,Naux,lambda,csi,mu,beta,alfa,gama,nu,&
				&tau,kk,kmedio,poptotal,soma,alfadp)
					
			end do

		end do

		do jj = 1,16

			do j = 1,17

				Naux(j,jj) = N(j,jj) + k(j,jj,2)*0.5d0
					
			end do

		end do

		do jj = 1,16

	        	soma = 0d0

	                do l = 1,16

		                soma = soma + kmedio(jj)*kk(jj,l)*(Naux(3,l) + Naux(4,l) + Naux(12,l) + Naux(13,l) + gama*Naux(9,l))/poptotal

	                end do

			do j = 1,17

				k(j,jj,3) = h*dt(j,jj,Naux,lambda,csi,mu,beta,alfa,gama,nu,&
				&tau,kk,kmedio,poptotal, soma,alfadp)
					
			end do

		end do

		do jj = 1,16

			do j = 1,17

				Naux(j,jj) = N(j,jj) + k(j,jj,3)
					
			end do

		end do

		do jj = 1,16

	        	soma = 0d0

	                do l = 1,16

		                soma = soma + kmedio(jj)*kk(jj,l)*(Naux(3,l) + Naux(4,l) + Naux(12,l) + Naux(13,l) + gama*Naux(9,l))/poptotal

	                end do

			do j = 1,17

				k(j,jj,4) = h*dt(j,jj,Naux,lambda,csi,mu,beta,alfa,gama,nu&
				&,tau,kk,kmedio,poptotal,soma,alfadp)
	
			end do

		end do

		do jj = 1,16

			do j = 1,17

				N(j,jj) = N(j,jj) + 1d0/6d0*(k(j,jj,1) + 2*k(j,jj,2) + 2*k(j,jj,3) + k(j,jj,4))
					
			end do

		end do

		do jj = 1,16

			do j = 1,17

				Nmedio(ii,j,jj) = Nmedio(ii,j,jj) + N(j,jj)*peso(amost)
					
			end do

		end do

	end do

     end do

   ti = t0

   !ordem t, 1, 2, 3, 4, 6, 5, 7, 11, 12, 13, 14, 15, 17, 8, 9, 16, 10

do ii = 1,HH

        do ifaixa = 1,16
    
    		write(ifaixa*10,'(18E20.8)') ti, (Nmedio(ii,1,ifaixa)), (Nmedio(ii,2,ifaixa)), (Nmedio(ii,3,ifaixa)),&
& (Nmedio(ii,4,ifaixa)), (Nmedio(ii,6,ifaixa)), (Nmedio(ii,5,ifaixa)), (Nmedio(ii,7,ifaixa)), (Nmedio(ii,11,ifaixa)), &
&(Nmedio(ii,12,ifaixa)), (Nmedio(ii,13,ifaixa)), (Nmedio(ii,14,ifaixa)), (Nmedio(ii,15,ifaixa)), (Nmedio(ii,17,ifaixa))&
&, (Nmedio(ii,8,ifaixa)), (Nmedio(ii,9,ifaixa)), (Nmedio(ii,16,ifaixa)), (Nmedio(ii,10,ifaixa))

        end do

    	write(170,'(18E20.8)') ti, sum(Nmedio(ii,1,:)), sum(Nmedio(ii,2,:)), sum(Nmedio(ii,3,:)),&
& sum(Nmedio(ii,4,:)), sum(Nmedio(ii,6,:)), sum(Nmedio(ii,5,:)), sum(Nmedio(ii,7,:)), sum(Nmedio(ii,11,:)), &
&sum(Nmedio(ii,12,:)), sum(Nmedio(ii,13,:)), sum(Nmedio(ii,14,:)), sum(Nmedio(ii,15,:)), sum(Nmedio(ii,17,:))&
&, sum(Nmedio(ii,8,:)), sum(Nmedio(ii,9,:)), sum(Nmedio(ii,16,:)), sum(Nmedio(ii,10,:))

	ti = ti + h

end do

end program sirs

real*8 function dt(j,jj,N,lambda,csi,mu,beta,alfa,gama,nu,tau,kk,kmedio,poptotal, soma,alfadp)
integer :: j, jj
real*8 :: lambda(2,16), csi(5,16), mu(2,16), beta(4,16), alfa(4,16), nu(3,16), tau(2,16), N(17,16), &
&Naux(17,16),kk(16,16),kmedio(16), alfadp(16)
real*8 :: soma, gama, poptotal

if (j .eq. 1) then			!S

	dt = -lambda(1,jj)*N(1,jj)*soma - csi(1,jj)*N(1,jj)

else if (j .eq. 2) then 		!E

	dt = lambda(1,jj)*N(1,jj)*soma - (csi(2,jj) + mu(1,jj))*N(2,jj)

else if (j .eq. 3) then			!A

	dt = mu(1,jj)*N(2,jj) - (csi(3,jj) + beta(1,jj) + beta(2,jj))*N(3,jj)

else if (j .eq. 4) then			!I


	dt = beta(2,jj)*N(3,jj) - (alfa(1,jj) + alfa(2,jj))*N(4,jj)

else if (j .eq. 5) then			!D

	dt = alfa(2,jj)*N(4,jj)

else if (j .eq. 6) then			!R

	dt = beta(1,jj)*N(3,jj) + alfa(1,jj)*N(4,jj) - csi(5,jj)*N(6,jj)

else if (j .eq. 7) then			!Sv

	dt = -lambda(2,jj)*N(7,jj)*soma + csi(1,jj)*N(1,jj) - nu(1,jj)*N(7,jj)

else if (j .eq. 8) then			!Sp

	dt = nu(2,jj)*N(17,jj) - tau(1,jj)*N(8,jj)*soma

else if (j .eq. 9) then			!Ip

	dt = tau(1,jj)*(N(8,jj)+N(17,jj))*soma - tau(2,jj)*N(9,jj) - alfadp(jj)*N(9,jj)

else if (j .eq. 10) then		!"Rp*" (Completamente imunizados vindos de P) 

	dt = nu(3,jj)*N(17,jj)

else if (j .eq. 11) then		!Ev

	dt = lambda(2,jj)*N(7,jj)*soma + csi(2,jj)*N(2,jj) - mu(2,jj)*N(11,jj)

else if (j .eq. 12) then		!Av

	dt = mu(2,jj)*N(11,jj) + csi(3,jj)*N(3,jj) - (beta(3,jj) + beta(4,jj))*N(12,jj)

else if (j .eq. 13) then		!Iv

	dt = beta(3,jj)*N(12,jj) - (alfa(3,jj)+alfa(4,jj))*N(13,jj)

else if (j .eq. 14) then

	dt = alfa(4,jj)*N(13,jj) + beta(4,jj)*N(12,jj) + csi(5,jj)*N(6,jj)

else if (j .eq. 15) then		!Rv

	dt = alfa(3,jj)*N(13,jj) + alfadp(jj)*N(9,jj)

else if (j .eq. 16) then		!"Rp**" (Vacinados protegidos recuperados vindos a partir do Ip)

	dt = tau(2,jj)*N(9,jj)

else if (j .eq. 17) then		!P

    dt = nu(1,jj)*N(7,jj) - (nu(2,jj) + nu(3,jj))*N(17,jj) - tau(1,jj)*N(17,jj)*soma

end if

return
end
