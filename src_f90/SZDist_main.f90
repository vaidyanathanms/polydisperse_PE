program SZdist
implicit none

real::PDI ! PDI (polydispersive index)
real::k ! Related to PDI [k=1/(PDI-1)]
real::range = 5 ! Max value of Sigma
real::area = 0 ! Area under the curve 
real::step ! Step size = range/end  
real::randnum ! random number from 0 to 1
real::average = 0 ! Sum of all chain MW, to be divided by num
real::tol = 10 ! tolerance value for PDI of generated chains (%)
real::PDIgen ! PDI of generated list
real::Mi = 0 ! used to calculate PDI of list
real::Mi2 = 0 ! used to calculate PDI of list

integer::i
integer::j
integer::l
integer::M ! Number avg MW of chains
integer,parameter::num = 100 ! Number of chains being created
integer,parameter::end = 100000 ! Number of steps 
integer,parameter::numlist = 10 ! Number of polymer lists to be created
integer::subtract = 1 ! Conditional: equal to 1 when the program is
! generating polymer, 0 when it has decided the M of the chain
integer::loop = 1 ! Similarly, will be equal to 1 until generated
! PDI is within the tolerance, then set to 0
integer::maxiteration = 100 ! sets the maximum amount of times the program will
! try to get a PDI within the tolerance before exiting
    
real,dimension(1:end)::S ! Ratios of length of chain to mean length   
real,dimension(1:end)::P ! Probability of SZ distribution 
real,dimension(1:end)::Normal ! Normalized distributation 
real,dimension(1:end)::IntNormal ! Integrated Areas for Normalized Distribution 
real,dimension(1:num)::MolWt ! List of MW for N polymers

Character(len = 20),dimension(numlist)::Filename
Character(len = 1)::FileNumberSmall
Character(len = 2)::FileNumberLarge
! Assign filenames to each list 
do i=1,numlist
   if (i .lt. 10) then
      write(FileNumberSmall,'(I1)') i
      Filename(i) = "Polymer"//FileNumberSmall//".dat"
   else if (i .lt. 100) then 
      write(FileNumberLarge,'(I2)') i 
      Filename(i) = "Polymer"//FileNumberLarge//".dat"
   end if
end do 

step = range/end 

! Prompt for values of PDW, Mn, and N
print *, "What is the PDI?" 
read *, PDI
print *, "What is Mn?"
read *, M

k = 1.0/real(PDI - 1.0)

! Sigma (S) can only range from 0 to infinity but not more than maybe 5?            
do i=1,end
   S(i) = (range/real(end))*i
end do

do i=1,end
   P(i) = (k**k)*(GAMMA(S(i))**-1)*(S(i)**(k-1))*(EXP(-1*k*S(i)))
end do

! Calculate total area for normalization
area = 0.5*(range/real(end))*P(1)

do i=2,end
   area = area + 0.5*(range/real(end))*(P(i)+P(i-1))
end do

! Normalize the probability function
do i=1,end
   Normal(i) = P(i)/area 
end do

! Calculate the normalized area for each slice such that the sum of normalized areas is one

IntNormal(1) = 0.5*Normal(1)*(range/real(end)) ! definte IntNormal(1) so that it doesnt index out of bounds for first calculation

do i=2,end
   IntNormal(i)=0.5*(Normal(i)+Normal(i-1))*(range/real(end))
end do

! Calculate random MW for each polymer using subtraction method

call random_seed()

! ====================== BEGIN GENERATING LISTS HERE ======================

do l=1,numList

   ! Resetting variables
   loop = 1
   Mi = 0
   Mi2 = 0

   print *, "Generating list", l
   do while (loop .eq. 1)

      ! Generates polymer list using subtraction method
      do i=1,num
         call random_number(randnum)
         subtract = 1
         j = 1
   
         do while (subtract == 1) 
            randnum = randnum - IntNormal(j)
      
            if (randnum .le.  0) then
               subtract = 0
               MolWt(i) = int(S(j-1)*M) 
            else      
               j = j + 1
               if(j == end+1) then
                  print *, 'array out of bounds error imminent'
                  exit  
               endif
            endif
         enddo
      end do
   
      ! Calculate PDI of list
      do i=1,num
         Mi2 = Mi2 + (MolWt(i)**2)
         Mi = Mi + Molwt(i)
      end do

      PDIgen = (Mi2*num)/(Mi**2)

      ! Checks if generated PDI is within tolerance of desired PDI and does
      ! not have an Mi smaller than 3

      if (ABS(PDIgen - PDI) .le. (PDI*tol)) then
         if (minval(MolWt) .ge. 3) then 
            loop = 0
         end if
      endif

   end do

   ! Calculates the average MW of the 100 chains to check if we got Mn back
   do i=1,num
      average = average + MolWt(i)
   end do

   average = average/num
   print *, 'the number average molecular wt is', average 
   print *, 'the PDI of the generated list is', PDIgen 

   ! End of loop

   open(unit = l, file = Filename(l))

   write(l,*),num
   write(l,*),Mi
   
   do i=1,num
      write(l,*),i, '      ', MolWt(i)
   end do

end do

! Prints the names of the generated lists
print *, "The generated files are:"

do i=1,numlist
   print *, Filename(i)
end do

end 
