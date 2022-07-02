from pyne import data, particle, nucname
import math
import matplotlib as plt
import csv
from matplotlib import pyplot as plt

# Subatomic Particle Masses in AMU
# 
electronMass = 0.00054858
neutronMass = 1.008665
protonMass = 1.007825

# Speed of Light in m/s
lightSpeed = 299792500

# Z: number of protons in the nucleus of a given element
def z(isotope):
	return isotope//10000000

# A: number of nucleons (neutrons + protons) in the nucleus of a given isotope
def a(isotope):
	return (isotope//10000)%1000
	
def state(isotope):
	return isotope%10000

def neutCount(isotope):
	return a(isotope) - z(isotope)


def massDefect(isotope):
	return z(isotope)*protonMass + neutCount(isotope)*neutronMass - data.atomic_mass(isotope)

def bindingEnergy(isotope):
	return (massDefect(isotope))*931.494893	
	
def bindingEnergyPerNucleon(isotope):
	return bindingEnergy(isotope)/a(isotope)
	

def readIsotopes():
	with open('input.txt') as file:
		for line in file:
			i = line.rstrip()
			print(i)
			be = bindingEnergy(nucname.id(i))
			print('BE: %.6f MeV' % be)
			bepn = bindingEnergyPerNucleon(nucname.id(i))
			print('BE Per Nucleon: %.6f MeV\n' % bepn)
	file.close()


# Empirical Mass Formula Coefficients (von Weizsacker, 1935)
#

#a_volume = 15.56 # 16 MeV
#a_surface = 17.23 # 18 MeV
#a_coulomb = 0.7 # 0.76 # MeV
#a_asymmetry = 23.285 # 23.5 # MeV
#a_pairing = 11 # MeV

a_volume = 16 # MeV
a_surface = 18 # MeV
a_coulomb = 0.72 # 0.76 # MeV
a_asymmetry = 23.5 # 23.5 # MeV
a_pairing = 11 # MeV


def semiEmpiricalBindingEnergy(isotope):

	p = z(isotope)
	nuc = a(isotope)
	neu = nuc - p
	
	term1 = a_volume*nuc
	term2 = a_surface*(pow(nuc,(2/3)))
	term3 = a_coulomb*(p*(p - 1))/(pow(nuc,(1/3)))
	term4 = a_asymmetry*(pow((neu - p),2)/nuc)
	#term5 = 0
	
	p_pairing = (p % 2)
	n_pairing = (neu % 2)
	
	if (n_pairing + p_pairing) == 0:
		term5 = a_pairing/pow(nuc, 0.5)
		
	if (n_pairing == 0) != (p_pairing == 0):
		term5 = 0
	
	if (n_pairing > 0) and (p_pairing > 0):
		term5 = -a_pairing/pow(nuc, 0.5)
		
	mass = term1 - term2 - term3 - term4 + term5
	
	return mass

def semf(isotope):
	
	p = z(isotope)
	neu = a(isotope) - p
	
	term1 = p*data.atomic_mass(nucname.id('H1'))
	term2 = neu*neutronMass
	term3 = semiEmpiricalBindingEnergy(isotope)/pow(lightSpeed, 2)
	
	semi = term1 + term2 - term3
	
	return semi
	

def BindingEnergyFromMassExcess(a, z, me):
	term1 = z*(protonMass - neutronMass)
	term2 = a*(neutronMass - 1)
	
	be = (term1 + term2)*931.494893 - me
	return be
	

def radCalc(filepath):	
	
	print("Isotope,A,Z,Excess Mass,BE,BE/a")
	
	with open(filepath, 'r') as f:
	
		reader = csv.reader(f)
		next(reader)
		for column in reader:
		
			i = column[0]
			em = float(column[1])
			be = BindingEnergyFromMassExcess(a(nucname.id(i)), z(nucname.id(i)), em) 
			be_per_n = be/a(nucname.id(i))
			#print("z: ", z(nucname.id(i)), "a: ", a(nucname.id(i)))
			print (i,",",a(nucname.id(i)),",",z(nucname.id(i)),",",em, ",",be,",",be_per_n)
			
	f.close()
	
def radCalcarray(filepath):	
	#print("Isotope,A,Z,Excess Mass,BE,BE/a")
	
	filepath2="/home/pyne-user/PycharmProjects/ne_calcs/output/output.csv"
	out = open(filepath2, 'w')

	
	with open(filepath, 'r') as f:
		
		
		writer = csv.writer(out)
		header = ['Isotope', 'A', 'Z', 'Excess Mass', 'BE', 'BE/a']
		writer.writerow(header)
		
		reader = csv.reader(f)
		next(reader)
		
		for column in reader:
		
			i = column[0]
			em = float(column[1])
			be = BindingEnergyFromMassExcess(a(nucname.id(i)), z(nucname.id(i)), em) 
			be_per_n = be/a(nucname.id(i))
			
			datarow = [i, a(nucname.id(i)), z(nucname.id(i)), em, be, be_per_n]
			writer.writerow(datarow)

			#print("z: ", z(nucname.id(i)), "a: ", a(nucname.id(i)))
			#print (i,",",a(nucname.id(i)),",",z(nucname.id(i)),",",em, ",",be,",",be_per_n)
		
			
	f.close()
	
	out.close()
	
def plotData(filepath):
	
	a_number = []
	excess_mass = []
	be_per_a = []
	errorlist = []
	
	with open(filepath, 'r') as f:
		plots = csv.reader(f, delimiter=',')
		next(f)
		for row in plots:
			
			an = float(row[1])
			emn = float(row[3])
			errorn = float(row[4]) - semiEmpiricalBindingEnergy(nucname.id(row[0]))
			beperan = float(row[5])
			
			a_number.append(an)
			excess_mass.append(emn)
			be_per_a.append(beperan)
			errorlist.append(errorn)
						
		#plt.plot(a_number,excess_mass,marker='o')
		#plt.show()
		
		#plt.plot(a_number,be_per_a,marker='o')
		#plt.show()
		
		plt.plot(a_number,errorlist, marker='o')
		plt.show()
		
radCalcarray('input.csv')
plotData("/home/pyne-user/PycharmProjects/ne_calcs/output/output.csv")
