#def calc_tm_all(sequences, oligo_conc, salt_conc):
#	oligo_concentration =  
#	salt_concentration = 
#	import os

def variance(x):
        x_sum = 0
        for entry in x:
                x_sum += entry
        mean = float(x_sum)/len(x)
        vari = 0
        for entry in x:
                vari += (entry-mean)*(entry-mean)
        vari = vari/float(len(x))
        return vari


def calc_tm_with_parameters(params, sequence, oligo_conc, salt_conc):
	from math import log
	from math import log10
	dHsum = 0
	dSsum = 0
	init_params = params['init']
	dHi = init_params[0]
	dSi = init_params[1]
	dSself = params['self_comp'][1]
	R = 1.987
	b = 4
	for i in range(len(sequence)-1):
		dinuc = sequence[i:(i+2)]
		dinuc_params = params[dinuc]
		dH = dinuc_params[0]
		dS = dinuc_params[1]
		dG = dinuc_params[2]
		dHsum += dH
		dSsum += dS
	Tm = (dHsum + dHi)/(dSsum + dSi + dSself + (R*log(oligo_conc/b))) + (16.6 * log10(salt_conc))
	return Tm

def calc_tm_for_long(sequence, salt_conc):
	from math import log10
	sequence = sequence.upper()
	A = sequence.count('A')
	T = sequence.count('T')
	G = sequence.count('G')
	C = sequence.count('C')

	return 81.5 + (41*float(G + C)/float(A+T+G+C)) - (500/float(A+T+G+C)) + 16.6*log10(salt_conc) + 273.15
	
def calc_tm_sugi(sequences, oligo_conc, salt_conc):

	output = []

	sequences = [sequence.upper() for sequence in sequences]
	parameters = {}
	parameters['AA'] = (-8000.0, -21.9, -1.2)
	parameters['AT'] = (-5600, -15.2, -0.9)
	parameters['AC'] = (-9400, -25.5, -1.5)
	parameters['AG'] = (-6600, -16.4, -1.5)
	
	parameters['TA'] = (-6600, -18.4, -0.9)
	parameters['TT'] = (-8000, -21.9, -1.2)
	parameters['TC'] = (-8800, -23.5, -1.5)
	parameters['TG'] = (-8200, -21.0, -1.7)
	
	parameters['CA'] = (-8200, -21.0, -1.7)
	parameters['CT'] = (-6600, -16.4, -1.5)
	parameters['CC'] = (-10900, -28.4, -2.1)
	parameters['CG'] = (-11800, -29.0, -2.8)

	parameters['GA'] = (-8800, -23.5, -1.5)
	parameters['GT'] = (-9400, -25.5, -1.5)
	parameters['GC'] = (-10500, -26.4, -2.3)
	parameters['GG'] = (-10900, -28.4, -2.1)

	parameters['init'] = (600, -9.0, 3.4)
	parameters['self_comp'] = (0.0, -1.4, 0.4)

	for sequence in sequences:
		if len(sequence) > 30 or len(sequence)<16:
			print 'WARNING: sequence is not optimal size for Tm prediction with NN method!'
		output.append(calc_tm_with_parameters(parameters, sequence, oligo_conc, salt_conc))

	return output

def calc_tm_long(sequences, salt_conc):
	output = []
	for sequence in sequences:
		if len(sequence) < 50:
			print 'WARNING: this sequence may be too short for the long salt-adjusted approximation'
		output.append(calc_tm_for_long(sequence, salt_conc))

	return output

def find_gaussian_parameters(values, expected):
#find the weighted average
	expected_sum = sum(expected)
	weighted_average = 0
	for i in range(len(values)):
		weighted_average += (float(expected[i])/expected_sum)*values[i]

	starting_omegas = [1, 2, 3, 4, 5, 6, 7, 8, 9]
	resids = find_residuals(starting_omegas, weighted_average, values, expected)
	best_primary_omega =  starting_omegas[resids.index(min(resids))]
	new_omegas = []
	for i in range(-5,6):
		new_omegas.append(best_primary_omega + float(i)/10.0)
	
	resids = find_residuals(new_omegas, weighted_average, values, expected)
	best_omega = new_omegas[resids.index(min(resids))]
	return [weighted_average, best_omega]
	

def test_omega(values, expected, omega, u):
	from pylab import lstsq
	from pylab import polyfit

	calculated_values = [calc_density_function(entry, u, omega) for entry in values]
	fit_params = polyfit(calculated_values, values, 1)
	normalization_factor = fit_params[0]
	

def calc_density_function(x, u, omega):
	from math import sqrt
	from math import pi
	from math import exp

	buffer = 1.0/(omega * sqrt(2*pi))
	buffer2 = -((x-u)*(x-u))/(2*omega*omega)
	return buffer * exp(buffer2)

def find_residuals(omegas, u, values, expected):
	from pylab import lstsq
	from pylab import polyfit

	output =[] 
	for omega in omegas:
		calculated_values = [calc_density_function(entry, u, omega) for entry in values]
		fit_params = polyfit(calculated_values, values, 1)
		normalization_factor = fit_params[0]
		residuals_total = 0
		for i in range(len(values)):
			value = values[i]
			dens_value = normalization_factor * calculated_values[i]
			residual = (expected[i] - dens_value)*(expected[i] - dens_value)
			residuals_total +=residual
		output.append(residuals_total)
	return output

def calc_gc_content(sequence):
	sequence = sequence.upper()
	num_g = sequence.count('G')
	num_c = sequence.count('C')
	gc = num_g + num_c
	return float(gc)/len(sequence)

def calc_expected_twist(sequence):
	seq_len = len(sequence)
	total_twist = (float(seq_len)*35.9)%360
	return total_twist

def make_binary_matrix(sequence):
	sequence = sequence.upper()
	matrix_dict = {}
	matrix_dict['A'] = [1,0,0,0]
	matrix_dict['T'] = [0,1,0,0]
	matrix_dict['C'] = [0,0,1,0]
	matrix_dict['G'] = [0,0,0,1]
	output = []
	for entry in sequence:
		for val in matrix_dict[entry]:
			output.append(val)
	return output

# Target Arm Sequence
# -fit calculated Tm to a gaussian, deviation from optimum

# Targetted Sequence
# find %GC content, fit that to a gaussian?

# Sequence Length
# - closeness to proper frequency plus offset
# - discrete function from length vs. reads



# Junction scores
