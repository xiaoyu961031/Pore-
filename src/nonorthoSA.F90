Program nonorthoSA
!
! This program calculates the accessible surface area by using a simple Monte
! Carlo integration for both orthorhombic and non-orthorhombic unit cells. 
! To run the program the following files are needed:
! 1) an input file containing the name of the file containing the diameters of 
!    the atoms, the name of the file containing the coordinates of the 
!    framework atoms, the diameter of the probe, the number of random trials
!    around each framework atom, the lenghts of the unit cells (a, b, c), the
!    cell angles (alpha, beta and gamma) and the crystal density.
! 2) A file containing the diameters of the framework atoms. The program is
!    written to read in an xyz file. Please note that the coordinates of all
!    framework atoms need to be specified (P1) as the program does not handle
!    symmetry operations.
!
! Things that you might have to change:
! - The maximum number of framework atoms 5000 at the moment. Change the value
!   of max_no if your unit cells contains more atoms.
!
! Changes to the code
!
! 2020-2021  Xiaoyu Wu
!            Append atom speices across periodic table for nanoporous materials informatic
!            Atoms radius collected from CCDC
!            Atoms types followed the order of Universal Forcefield
!
! 2009-2013  Linjiang Chen
!            Decompose the total SA into contributions from different chemical spieces, 
!            currently defined as different partial charges
!
! 9/3/2009   Tina Duren
!            fixed a bug that wrote the error message to a file rather than
!            the screen if an atom type was not specified in the diameter list
!            Thanks to Kenji Sumida from the University of Berkeley for 
!            pointing this out.
!  	
Use matrix
Use defaults
        
Implicit None

!
! Maximum number of atoms in your unit cell. Increase if your unit cell 
! contains more atoms.
!
Integer, Parameter   :: max_no = 50000

Character(len = 3)   :: char
Character(len = 10)  :: chargetype(1:100)
Character(len = 3)   :: atom
Character(len = 10)  :: symbol
Character(len = 10)  :: atomname(1:max_no), atomtype(1:100)
Character(len = 100) :: charge_file, atom_file, coord_file, output_file

Type(MatrixType)     :: lij, lijinv

Real(kind=RDbl)      :: chargevalue(1:max_no)
Real(kind=RDbl)      :: sigmatype(1:100), atomsigma(1:5000)
Real(kind=RDbl)      :: x(1:max_no), y(1:max_no), z(1:max_no), charge(1:max_no)
Real(kind=RDbl)      :: xnew(1:max_no), ynew(1:max_no), znew(1:max_no)
Real(kind=RDbl)      :: alpha, beta, gamma, alpha_r, beta_r, gamma_r
Real(kind=RDbl)      :: a, b, c, ai, bi, ci, aj, bj, cj, ak, bk, ck
Real(kind=RDbl)      :: xpoint, ypoint, zpoint, xpointnew, ypointnew, zpointnew
Real(kind=RDbl)      :: phi, costheta, theta
Real(kind=RDbl)      :: rho_crys, dprobe, Nsample
Real(kind=RDbl)      :: dx, dy, dz, dx_trans, dy_trans, dz_trans, dist2
Real(kind=RDbl)      :: sjreal, stotal, sfrac, uc_volume, stotalreduced
Real(kind=RDbl)      :: stotal_char(1:max_no), stotalreduced_char(1:max_no)
Real(kind=RDbl)      :: ran_num1, ran_num2, trash
Real(kind=RDBL)      :: ran0

Integer              :: seed, N, ntypes, ntypes_char, i, j, k, ncount
Integer              :: ierror
INTEGER              :: i1, i2, i3
INTEGER              :: w_char

Logical              :: match, deny
!
! Seed for the random number generator, change if you want to start from a
! different random number
!
seed = -20170601
!
! to properly initialise the random number generator, a negative seed is
! needed
!
If (seed> 0) seed = -1 * seed
!
! Read the data from the input file
!
Read(*,('(A)')) charge_file      ! file containing the different chargetypes and chargevalues
Read(*,('(A)')) atom_file        ! file containing the diameters of atoms
Read(*,('(A)')) coord_file       ! file containing the cartesian coordinates
Read(*,*) dprobe                 ! diameter of probe in A
Read(*,*) Nsample                ! Number of trials per framework atom
Read(*,*) a, b, c                ! Cell parameters in A
Read(*,*) alpha, beta, gamma     ! cell angles
Read(*,*) rho_crys               ! density of crystal in g / cm3


alpha_r = alpha*degtorad
beta_r = beta*degtorad
gamma_r = gamma*degtorad

ai = a
aj = 0
ak = 0
bi = b*cos(gamma_r)
bj = b*sin(gamma_r)
bk = 0
ci = c*cos(beta_r)
cj = (b*c*cos(alpha_r)-bi*ci)/bj
ck = sqrt(c**2-ci**2-cj**2)

lij%comp = 0
lij%comp(1,1) = ai
lij%comp(1,2) = bi
lij%comp(1,3) = ci
lij%comp(2,1) = aj
lij%comp(2,2) = bj
lij%comp(2,3) = cj
lij%comp(3,1) = ak
lij%comp(3,2) = bk
lij%comp(3,3) = ck

lijinv=matrix_inverse(lij)

!
! Open the file containing the charge types
!
Open(30, file=charge_file, status='old', IOSTAT = ierror)
if (ierror /= 0) Then
    Write(*,*) 'Error opening file: ',charge_file
    Write(*,*) 'You sucked and try again!'
END IF

ntypes_char = 0

Do
   Read(30,*) char
   IF(Trim(char) == 'EOF') EXIT
   ntypes_char = ntypes_char + 1
END DO

Rewind(30)

Do i = 1, ntypes_char
   Read(30,*) chargetype(i), chargevalue(i)
End do

    Write(*,*) 'ntypes_char = ', ntypes_char
!    Write(*,*) 'chargetype(3) = ', chargetype(3)
!    Write(*,*) 'chargevalue(3) = ', chargevalue(3)

!
! Open the files and read the data
!
Open(10, file=atom_file, status='old', IOSTAT = ierror)
If (ierror /= 0) Then
    Write(*,*) 'Error opening file: ',atom_file
    Write(*,*) 'Make sure it exists or check the spelling of the filename'
END IF


ntypes = 0

Do
   Read(10,*) atom
   IF(Trim(atom) == 'EOF') EXIT
   ntypes = ntypes + 1
END DO

Rewind(10)

Do i = 1, ntypes
   Read(10,*) atomtype(i), sigmatype(i)
End do

!
! Openeing and reading the coord file
!
Open(20, file=coord_file, status='old', IOSTAT = ierror)
If (ierror /= 0) Then
    Write(*,*) 'Error opening file: ',coord_file
    Write(*,*) 'Make sure it exists or check the spelling of the filename'
END IF

!
! Skip first 3 lines
!
Read(20,*)
Read(20,*)
Read(20,*)

!
! The 4th line contains the number of framework atoms
!
Read(20,*) N

!
! Skip the blank line of the xyz file
!
!Read(20,*)
Do i=1, N

! The format of the coordinate file corresponds to an xyz file as exported
! from Diamond.
! x(i), y(i), z(i): x,y,z coordinates in Angstrom
! symbol: chemical symbol of framework atom.
!
! Change the following read statement if you want to use a different
! file format. But ensure that you end up with the coordinates in A and
! the name of the framework atoms!
!
   Read(20,*) i1, x(i), y(i), z(i), symbol, charge(i), i2, i3

! The diameters of the framework atoms are stored together with the atom
! name. So assign atom name according to chemical symbol.
!
  symbol = Trim(symbol)
  SELECT CASE(symbol)
    CASE('Ac')
      atomname(i) = 'Actinium'
    CASE('Actinium')
      atomname(i) = 'Actinium'
    CASE('Ag')
      atomname(i) = 'Silver'
    CASE('Silver')
      atomname(i) = 'Silver'
    CASE('Al')
      atomname(i) = 'Aluminum'
    CASE('Aluminum')
      atomname(i) = 'Aluminum'
    CASE('Am')
      atomname(i) = 'Americium'
    CASE('Americium')
      atomname(i) = 'Americium'
    CASE('As')
      atomname(i) = 'Arsenic'
    CASE('Arsenic')
      atomname(i) = 'Arsenic'
    CASE('At')
      atomname(i) = 'Astatine'
    CASE('Astatine')
      atomname(i) = 'Astatine'
    CASE('Au')
      atomname(i) = 'Gold'
    CASE('Gold')
      atomname(i) = 'Gold'
    CASE('B')
      atomname(i) = 'Boron'
    CASE('Boron')
      atomname(i) = 'Boron'
    CASE('Ba')
      atomname(i) = 'Barium'
    CASE('Barium')
      atomname(i) = 'Barium'
    CASE('Be')
      atomname(i) = 'Beryllium'
    CASE('Beryllium')
      atomname(i) = 'Beryllium'
    CASE('Bi')
      atomname(i) = 'Bismuth'
    CASE('Bismuth')
      atomname(i) = 'Bismuth'
    CASE('Bk')
      atomname(i) = 'Berkelium'
    CASE('Berkelium')
      atomname(i) = 'Berkelium'
    CASE('Br')
      atomname(i) = 'Bromine'
    CASE('Bromine')
      atomname(i) = 'Bromine'
    CASE('C')
      atomname(i) = 'Carbon'
    CASE('Carbon')
      atomname(i) = 'Carbon'
    CASE('Ca')
      atomname(i) = 'Calcium'
    CASE('Calcium')
      atomname(i) = 'Calcium'
    CASE('Cd')
      atomname(i) = 'Cadmium'
    CASE('Cadmium')
      atomname(i) = 'Cadmium'
    CASE('Ce')
      atomname(i) = 'Cerium'
    CASE('Cerium')
      atomname(i) = 'Cerium'
    CASE('Cf')
      atomname(i) = 'Californium'
    CASE('Californium')
      atomname(i) = 'Californium'
    CASE('Cl')
      atomname(i) = 'Chlorine'
    CASE('Chlorine')
      atomname(i) = 'Chlorine'
    CASE('Cm')
      atomname(i) = 'Curium'
    CASE('Curium')
      atomname(i) = 'Curium'
    CASE('Co')
      atomname(i) = 'Cobalt'
    CASE('Cobalt')
      atomname(i) = 'Cobalt'
    CASE('Cr')
      atomname(i) = 'Chromium'
    CASE('Chromium')
      atomname(i) = 'Chromium'
    CASE('Cu')
      atomname(i) = 'Cesium'
    CASE('Cesium')
      atomname(i) = 'Cesium'
    CASE('Cs')
      atomname(i) = 'Copper'
    CASE('Copper')
      atomname(i) = 'Copper'
    CASE('Dy')
      atomname(i) = 'Dysprosium'
    CASE('Dysprosium')
      atomname(i) = 'Dysprosium'
    CASE('Eu')
      atomname(i) = 'Erbium'
    CASE('Erbium')
      atomname(i) = 'Erbium'
    CASE('Er')
      atomname(i) = 'Einsteinium'
    CASE('Einsteinium')
      atomname(i) = 'Einsteinium'
    CASE('Es')
      atomname(i) = 'Europium'
    CASE('Europium')
      atomname(i) = 'Europium'
    CASE('F')
      atomname(i) = 'Fluorine'
    CASE('Fluorine')
      atomname(i) = 'Fluorine'
    CASE('Fe')
      atomname(i) = 'Iron'
    CASE('Iron')
      atomname(i) = 'Iron'
    CASE('Fm')
      atomname(i) = 'Fermium'
    CASE('Fermium')
      atomname(i) = 'Fermium'
    CASE('Fr')
      atomname(i) = 'Francium'
    CASE('Francium')
      atomname(i) = 'Francium'
    CASE('Ga')
      atomname(i) = 'Gallium'
    CASE('Gallium')
      atomname(i) = 'Gallium'
    CASE('Ge')
      atomname(i) = 'Gadolinium'
    CASE('Gadolinium')
      atomname(i) = 'Gadolinium'
    CASE('Gd')
      atomname(i) = 'Germanium'
    CASE('Germanium')
      atomname(i) = 'Germanium'
    CASE('H')
      atomname(i) = 'Hydrogen'
    CASE('Hydrogen')
      atomname(i) = 'Hydrogen'
    CASE('Hf')
      atomname(i) = 'Hafnium'
    CASE('Hafnium')
      atomname(i) = 'Hafnium'
    CASE('Hg')
      atomname(i) = 'Mercury'
    CASE('Mercury')
      atomname(i) = 'Mercury'
    CASE('Ho')
      atomname(i) = 'Holmium'
    CASE('Holmium')
      atomname(i) = 'Holmium'
    CASE('I')
      atomname(i) = 'Iodine'
    CASE('Iodine')
      atomname(i) = 'Iodine'
    CASE('In')
      atomname(i) = 'Indium'
    CASE('Indium')
      atomname(i) = 'Indium'
    CASE('Ir')
      atomname(i) = 'Iridium'
    CASE('Iridium')
      atomname(i) = 'Iridium'
    CASE('K')
      atomname(i) = 'Potassium'
    CASE('Potassium')
      atomname(i) = 'Potassium'
    CASE('Kr')
      atomname(i) = 'Krypton'
    CASE('Krypton')
      atomname(i) = 'Krypton'
    CASE('La')
      atomname(i) = 'Lanthanum'
    CASE('Lanthanum')
      atomname(i) = 'Lanthanum'
    CASE('Li')
      atomname(i) = 'Lithium'
    CASE('Lithium')
      atomname(i) = 'Lithium'
    CASE('Lu')
      atomname(i) = 'Lawrencium'
    CASE('Lawrencium')
      atomname(i) = 'Lawrencium'
    CASE('Lr')
      atomname(i) = 'Lutetium'
    CASE('Lutetium')
      atomname(i) = 'Lutetium'
    CASE('Md')
      atomname(i) = 'Mendelevium'
    CASE('Mendelevium')
      atomname(i) = 'Mendelevium'
    CASE('Mg')
      atomname(i) = 'Magnesium'
    CASE('Magnesium')
      atomname(i) = 'Magnesium'
    CASE('Mn')
      atomname(i) = 'Manganese'
    CASE('Manganese')
      atomname(i) = 'Manganese'
    CASE('Mo')
      atomname(i) = 'Molybdenum'
    CASE('Molybdenum')
      atomname(i) = 'Molybdenum'
    CASE('N')
      atomname(i) = 'Nitrogen'
    CASE('Nitrogen')
      atomname(i) = 'Nitrogen'
    CASE('Na')
      atomname(i) = 'Sodium'
    CASE('Sodium')
      atomname(i) = 'Sodium'
    CASE('Ne')
      atomname(i) = 'Niobium'
    CASE('Niobium')
      atomname(i) = 'Niobium'
    CASE('Nb')
      atomname(i) = 'Neodymium'
    CASE('Neodymium')
      atomname(i) = 'Neodymium'
    CASE('Nd')
      atomname(i) = 'Neon'
    CASE('Neon')
      atomname(i) = 'Neon'
    CASE('Ni')
      atomname(i) = 'Nickel'
    CASE('Nickel')
      atomname(i) = 'Nickel'
    CASE('No')
      atomname(i) = 'Nobelium'
    CASE('Nobelium')
      atomname(i) = 'Nobelium'
    CASE('Np')
      atomname(i) = 'Neptunium'
    CASE('Neptunium')
      atomname(i) = 'Neptunium'
    CASE('O')
      atomname(i) = 'Oxygen'
    CASE('Oxygen')
      atomname(i) = 'Oxygen'
    CASE('Os')
      atomname(i) = 'Osmium'
    CASE('Osmium')
      atomname(i) = 'Osmium'
    CASE('P')
      atomname(i) = 'Phosphorus'
    CASE('Phosphorus')
      atomname(i) = 'Phosphorus'
    CASE('Pa')
      atomname(i) = 'Protactinium'
    CASE('Protactinium')
      atomname(i) = 'Protactinium'
    CASE('Pb')
      atomname(i) = 'Lead'
    CASE('Lead')
      atomname(i) = 'Lead'
    CASE('Pd')
      atomname(i) = 'Palladium'
    CASE('Palladium')
      atomname(i) = 'Palladium'
    CASE('Pm')
      atomname(i) = 'Promethium'
    CASE('Promethium')
      atomname(i) = 'Promethium'
    CASE('Po')
      atomname(i) = 'Polonium'
    CASE('Polonium')
      atomname(i) = 'Polonium'
    CASE('Pr')
      atomname(i) = 'Praseodymium'
    CASE('Praseodymium')
      atomname(i) = 'Praseodymium'
    CASE('Pt')
      atomname(i) = 'Platinum'
    CASE('Platinum')
      atomname(i) = 'Platinum'
    CASE('Pu')
      atomname(i) = 'Plutonium'
    CASE('Plutonium')
      atomname(i) = 'Plutonium'
    CASE('Ra')
      atomname(i) = 'Radium'
    CASE('Radium')
      atomname(i) = 'Radium'
    CASE('Rb')
      atomname(i) = 'Rubidium'
    CASE('Rubidium')
      atomname(i) = 'Rubidium'
    CASE('Re')
      atomname(i) = 'Rhenium'
    CASE('Rhenium')
      atomname(i) = 'Rhenium'
    CASE('Rh')
      atomname(i) = 'Rhodium'
    CASE('Rhodium')
      atomname(i) = 'Rhodium'
    CASE('Rn')
      atomname(i) = 'Radon'
    CASE('Radon')
      atomname(i) = 'Radon'
    CASE('Ru')
      atomname(i) = 'Ruthenium'
    CASE('Ruthenium')
      atomname(i) = 'Ruthenium'
    CASE('S')
      atomname(i) = 'Sulfur'
    CASE('Sulfur')
      atomname(i) = 'Sulfur'
    CASE('Sb')
      atomname(i) = 'Antimony'
    CASE('Antimony')
      atomname(i) = 'Antimony'
    CASE('Sc')
      atomname(i) = 'Scandium'
    CASE('Scandium')
      atomname(i) = 'Scandium'
    CASE('Se')
      atomname(i) = 'Selenium'
    CASE('Selenium')
      atomname(i) = 'Selenium'
    CASE('Si')
      atomname(i) = 'Silicon'
    CASE('Silicon')
      atomname(i) = 'Silicon'
    CASE('Sm')
      atomname(i) = 'Samarium'
    CASE('Samarium')
      atomname(i) = 'Samarium'
    CASE('Sn')
      atomname(i) = 'Tin'
    CASE('Tin')
      atomname(i) = 'Tin'
    CASE('Sr')
      atomname(i) = 'Strontium'
    CASE('Strontium')
      atomname(i) = 'Strontium'
    CASE('Ta')
      atomname(i) = 'Tantalum'
    CASE('Tantalum')
      atomname(i) = 'Tantalum'
    CASE('Tb')
      atomname(i) = 'Terbium'
    CASE('Terbium')
      atomname(i) = 'Terbium'
    CASE('Tc')
      atomname(i) = 'Technetium'
    CASE('Technetium')
      atomname(i) = 'Technetium'
    CASE('Te')
      atomname(i) = 'Tellurium'
    CASE('Tellurium')
      atomname(i) = 'Tellurium'
    CASE('Th')
      atomname(i) = 'Thorium'
    CASE('Thorium')
      atomname(i) = 'Thorium'
    CASE('Ti')
      atomname(i) = 'Titanium'
    CASE('Titanium')
      atomname(i) = 'Titanium'
    CASE('TI')
      atomname(i) = 'Thallium'
    CASE('Thallium')
      atomname(i) = 'Thallium'
    CASE('Tm')
      atomname(i) = 'Thulium'
    CASE('Thulium')
      atomname(i) = 'Thulium'
    CASE('U')
      atomname(i) = 'Uranium'
    CASE('Uranium')
      atomname(i) = 'Uranium'
    CASE('V')
      atomname(i) = 'Vanadium'
    CASE('Vanadium')
      atomname(i) = 'Vanadium'
    CASE('W')
      atomname(i) = 'Tungsten'
    CASE('Tungsten')
      atomname(i) = 'Tungsten'
    CASE('Y')
      atomname(i) = 'Yttrium'
    CASE('Yttrium')
      atomname(i) = 'Yttrium'
    CASE('Yb')
      atomname(i) = 'Ytterbium'
    CASE('Ytterbium')
      atomname(i) = 'Ytterbium'
    CASE('Zn')
      atomname(i) = 'Zinc'
    CASE('Zinc')
      atomname(i) = 'Zinc'
    CASE('Zr')
      atomname(i) = 'Zirconium'
    CASE('Zirconium')
      atomname(i) = 'Zirconium'
    CASE default
       Write(*,*) symbol ,' not in list'
       Write(*,*) 'List can be found in nonorthoSA.F90'
  END SELECT

! Transform molecules coordinates
   xnew(i) = x(i)*lijinv%comp(1,1) + y(i)*lijinv%comp(1,2) + &
             z(i)*lijinv%comp(1,3)
   ynew(i) = x(i)*lijinv%comp(2,1) + y(i)*lijinv%comp(2,2) + &
             z(i)*lijinv%comp(2,3)
   znew(i) = x(i)*lijinv%comp(3,1) + y(i)*lijinv%comp(3,2) + &
             z(i)*lijinv%comp(3,3)
   xnew(i) = a*xnew(i)
   ynew(i) = b*ynew(i)
   znew(i) = c*znew(i)

! Translate the transformed coordinates so they lie in the box a,b,c

    If(xnew(i)<0.0) xnew(i) = xnew(i) + a
    If(xnew(i)>=a) xnew(i) = xnew(i) - a
    If(ynew(i)<0.0) ynew(i) = ynew(i) + b
    If(ynew(i)>=b) ynew(i) = ynew(i) - b
    If(znew(i)<0.0) znew(i) = znew(i) + c
    If(znew(i)>=c) znew(i) = znew(i) - c

! Match sigmas with coordinates

    atomname(i)=Trim(atomname(i))

    match=.False.

    Do j=1, ntypes
       If(atomname(i)==atomtype(j)) Then
          atomsigma(i)=sigmatype(j)+dprobe
          match=.True.
          Exit
       End If
    End Do

    If(.Not.match) Then
      Write(*,*) 'Could not find match for atom: ', i, ' ', atomname(i)
      Write(*,*) 'The name is either read in incorrectly or does not exist'
      Write(*,*) 'in the list of available atoms in ',atom_file
      Stop
    End If

 End Do

Write(*,*)
Write(*,*) 'Calculating the accessible surface area for the following input parameters'
Write(*,*) '-------------------------------------------------------------------------'
Write(*,*) 'File with framework coordinates: ',TRIM(coord_file)
Write(*,*) 'File with atom diameters: ', TRIM(atom_file)
Write(*,'(A,F12.3)') ' Probe diameter in A: ',dprobe

! Main sampling cycle

stotal=0.0

 DO w_char=1, ntypes_char
    stotal_char(w_char) = 0.0
 END DO

Do i=1, N ! Loop over all framework atoms

    ncount=0

    Do j=1, Nsample ! Number of trial positions around each framework atom
!
! Generate random vector of length 1
! First generate phi 0:+2pi
!
       phi = ran0(seed)*twopi

! Generate cosTheta -1:1 to allow for even distribution of random vectors
! on the unit sphere.See http://mathworld.wolfram.com/SpherePointPicking.html
! for further explanations
!
       costheta = 1 - ran0(seed) * 2.0
       theta = Acos(costheta)
       xpoint = sin(theta)*cos(phi)
       ypoint = sin(theta)*sin(phi)
       zpoint = costheta


! Make this vector of (sigma+probe_diameter)/2.0 length

       xpoint=xpoint*atomsigma(i)/2.0
       ypoint=ypoint*atomsigma(i)/2.0
       zpoint=zpoint*atomsigma(i)/2.0

! Transform random vector

       xpointnew = xpoint*lijinv%comp(1,1) + ypoint*lijinv%comp(1,2) &
                   + zpoint*lijinv%comp(1,3)
       ypointnew = xpoint*lijinv%comp(2,1) + ypoint*lijinv%comp(2,2) &
                   + zpoint*lijinv%comp(2,3)
       zpointnew = xpoint*lijinv%comp(3,1) + ypoint*lijinv%comp(3,2) &
                   + zpoint*lijinv%comp(3,3)
       xpointnew = a*xpointnew
       ypointnew = b*ypointnew
       zpointnew = c*zpointnew

! Translate the center of coordinate to the particle i center and apply PBC

       xpointnew = xpointnew + xnew(i)
       ypointnew = ypointnew + ynew(i)
       zpointnew = zpointnew + znew(i)

       If(xpointnew < 0.0) xpointnew = xpointnew + a
       If(xpointnew >= a) xpointnew = xpointnew - a
       If(ypointnew < 0.0) ypointnew = ypointnew + b
       If(ypointnew >= b) ypointnew = ypointnew - b
       If(zpointnew < 0.0) zpointnew = zpointnew + c
       If(zpointnew >= c) zpointnew = zpointnew - c

! Now we check for overlap

       deny=.False.

       Do k=1,N
          if(k==i) cycle
          dx = xpointnew - xnew(k)
          dx = dx - a*int(2.0*dx/a)
          dy = ypointnew - ynew(k)
          dy = dy - b*int(2.0*dy/b)
          dz = zpointnew-znew(k)
          dz = dz-c*int(2.0*dz/c)

          dx = dx / a
          dy = dy / b
          dz = dz / c

          dx_trans = dx*lij%comp(1,1) + dy*lij%comp(1,2) + dz*lij%comp(1,3)
          dy_trans = dx*lij%comp(2,1) + dy*lij%comp(2,2) + dz*lij%comp(2,3)
          dz_trans = dx*lij%comp(3,1) + dy*lij%comp(3,2) + dz*lij%comp(3,3)

          dist2=dx_trans*dx_trans+dy_trans*dy_trans+dz_trans*dz_trans

          If(sqrt(dist2)<0.999*atomsigma(k)/2.0) then
             deny=.True.
             Exit
          End If
       End Do

       If(deny) Cycle

       ncount=ncount+1

    End Do

! Fraction of the accessible surface area for sphere i

    sfrac=Real(ncount)/Real(Nsample)

! Surface area for sphere i in real units (A^2)

    sjreal=pi*atomsigma(i)*atomsigma(i)*sfrac
    stotal=stotal+sjreal

    DO w_char=1, ntypes_char
       IF (charge(i) == chargevalue(w_char)) THEN
          stotal_char(w_char) = stotal_char(w_char) + sjreal
       END IF
    END DO

End Do

! Converting stotal on Surface per Volume
! Unit volume calculated from the absolute value of the determinate of the 3x3
! matrix containing the vectors defining the unit cell.
! See e.g. Marsden, et al. "Basic Multivariable Calculus" 1993 pg. 53

uc_volume=abs(ai*(bj*ck-bk*cj)-aj*(bi*ck-bk*ci)+ak*(bi*cj-bj*ci))
stotalreduced=stotal/uc_volume*1.E4

! Report results

 DO w_char=1, ntypes_char
    stotalreduced_char(w_char) = stotal_char(w_char)/uc_volume*1.E4
    Write(*,'(A,I5,A5,F12.2)') ' Total surface area (m^2/cm^3) of chargetype: ', &
                              w_char, 'is', stotalreduced_char(w_char)
    Write(*,'(A,I5,A5,F12.2)') ' Total surface area (m^2/g) of chargetype: ', &
                              w_char, 'is', stotalreduced_char(w_char)/rho_crys
 END DO

Write(*,'(A,F12.2)') ' Total surface area in Angstroms^2: ', stotal
Write(*,'(A,F12.2)') ' Total surface area per volume in m^2/cm^3: ', &
                      stotalreduced
Write(*,'(A,F12.2)') ' Total surface area per volume in m^2/g: ', &
                            stotalreduced / rho_crys


End Program NonorthoSA

!----------------FUNCTIONS-------------------------------------

FUNCTION ran0(idum)
!
! Random number generator from W.H. Press et al, Numerical Recipes in
! FORTRAN, Cambridge University Press, 1992
!
 INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
 REAL(kind = 8) ran0,AM,EPS,RNMX
 PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
   NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
 INTEGER j,k,iv(NTAB),iy
 
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
 if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do 11 j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11   continue
     iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran0=min(AM*iy,RNMX)
return

END FUNCTION ran0
