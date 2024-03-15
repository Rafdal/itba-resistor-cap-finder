import math

class ResistorFinder:
	def __init__(self):
		self.serieE12 = [1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2]
		self.serieE24 = [1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1]
		self.serieE96 = [1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30, 1.33, 1.37, 1.40, 1.43, 1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74, 1.78, 1.82, 1.87, 1.91, 1.96, 2.00, 2.05, 2.10, 2.15, 2.21, 2.26, 2.32, 2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09, 3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12, 4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49, 5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81, 6.98, 7.15, 7.32, 7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76]
		  
		self.Rbase = self.serieE12	# HARD CODED: E12 series [!][!][!]
		  
		self.R = []	# resistance combinations array
		self.G = []	# conductance combinations array
		self.n_max = 0	# maximum valid index
		self.out_r1 = [None]*20	# output array for R1
		self.out_r2 = [None]*20	# output array for R2
		self.out_op = [""]*20	# output array for operation
		self.out_rres = [None]*20	# output array for resulting resistance
		self.out_tol = [None]*20	# output array for resulting tolerance
		self.r1 = 0	# first resistor
		self.r2 = 0	# second resistor
		self.r1_idx = 0
		self.rres = 0
		self.rres_tol = 0
		self.best_tol = 0
		self.out_idx = 0
		self.op = ""
  
		self.out_prres, self.out_vrres = 0, 0
		self.select_series()
  
	def select_series(self, exp = [0, 7]):
     
		# Create the resistance array
		self.R = []
		for mult in range(exp[0], exp[1]):
			for idx in range(len(self.Rbase)):
				# Need to round to compensate for pow() errors; allow max two decimals, needed for E96
				self.R.append(round(self.Rbase[idx] * math.pow(10, mult) * 100, 2) / 100)
		
		# Calculate the maximum valid index
		self.n_max = len(self.R) - 1
		
		# Compute the conductances array, lowest conductance first to have an array sorted in ascending order
		self.G = [1.0 / self.R[self.n_max - idx] for idx in range(0, self.n_max + 1)]

	def find_index(self, vect, value):
		index_min = 0
		index_max = self.n_max + 1
		index = math.floor((index_min + index_max) / 2.0)
		i = 0

		while ((index_max - index_min) > 1) and (i < 2000):
			if vect[index] == value:
				break
			elif vect[index] > value:
				index_max = index
			elif vect[index] < value:
				index_min = index

			index = math.floor((index_min + index_max) / 2.0)
			i += 1

		if index < self.n_max:
			tol1 = abs(vect[index] / (value + 1e-30) - 1.0)
			tol2 = abs(vect[index + 1] / (value + 1e-30) - 1.0)
			if tol1 < tol2:
				return index
			else:
				return index + 1
		else:
			return index

	def roundWithMultiplier(self, value, offset = 4):
		#			  0   1   2   3   4  5   6   7   8   9   10  11  12
		#			-12	 -9	 -6	 -3	  0	 3	 6	 9   12  15  18  21  24
		highMults = ['p','n','u','m','','K','M','G','T','P','E','Z','Y']
		multCount = 0
		ogVal = value
		if(value == 0):
			return ["0", "0"]
		elif(value >= 0.999):
			while (value // 999 > 0):
				value = value / 1000.0
				multCount += 1
		elif(value < 0.0999):
			while (value < 0.0999):
				value = value * 1000.0
				multCount -= 1

		if multCount < (-12) or multCount > 7:
			print("OUT OF RANGE:", multCount)
			multCount = 0
			
		ogVal = str(round(ogVal,2)).rstrip('0').rstrip('.')
		return [str(round(value,2)).rstrip('0').rstrip('.') + highMults[offset + multCount], ogVal]


	def compute(self, rd, options = "both", scale = 1.0):
		if options not in ["both", "series", "parallel"]:
			raise ValueError("Invalid options value")

		r1, r2, r1_idx, rres, rres_tol, out_idx = 0, 0, 0, 0, 0, 0

		out_prres, out_vrres = 0, 0
		iter = 0  # number of iterations

		rd = rd * scale  # scale the input value

		# compute assuming resistors in series
		# locate nearest approximation with standard resistor values
		r1_idx = self.find_index(self.R, rd)
		r1 = self.R[r1_idx]
		# other resistor
		r2 = 0
		rres = r1
		rres_tol = (rres - rd) / rd  # relative tolerance
		best_tol = rres_tol

		if options == "both" or options == "series":
			out_idx = 0
			self.out_r1[out_idx] = r1
			self.out_r2[out_idx] = r2
			self.out_op[out_idx] = "+"
			self.out_rres[out_idx] = rres
			self.out_tol[out_idx] = rres_tol
			out_idx += 1

			while self.R[r1_idx] >= rd / 2.0:
				iter += 1
				r1 = self.R[r1_idx]

				r2d = abs(rd - r1)  # this is the value needed
				if r2d < 0:  # might happen...
					continue

				r2_idx = self.find_index(self.R, r2d)
				r2 = self.R[r2_idx]  # get the nearest standard value
				rres = r1 + r2  # compute the resulting composition
				rres_tol = rres / rd - 1.0  # and its tolerance

				if abs(rres_tol) < abs(best_tol):
					self.out_r1[out_idx] = r1
					self.out_r2[out_idx] = r2
					self.out_op[out_idx] = "+"
					self.out_rres[out_idx] = rres
					self.out_tol[out_idx] = rres_tol
					out_idx += 1
		
				r1_idx -= 1

		rd = 1.0 / rd
  
		if options == "both" or options == "parallel":
      
			# compute assuming resistors in parallel
			r1_idx = self.find_index(self.G, rd)
			while self.G[r1_idx] >= rd / 2.001 and iter < 10000:
				iter += 1
				r1 = self.G[r1_idx]

				r2d = abs(rd - r1)  # this is the value needed
				if r2d < 0:  # might happen...
					continue

				r2_idx = self.find_index(self.G, r2d)
				r2 = self.G[r2_idx]  # get the nearest standard value
				rres = r1 + r2  # compute the resulting composition
				rres_tol = rd / rres - 1.0  # and its tolerance

				if abs(rres_tol) < abs(best_tol):
					self.out_r1[out_idx] = self.R[self.n_max - r1_idx]  # 1.0 / r1;
					self.out_r2[out_idx] = self.R[self.n_max - r2_idx]  # 1.0 / r2;
					self.out_op[out_idx] = "||"
					self.out_rres[out_idx] = 1.0 / rres
					self.out_tol[out_idx] = rres_tol
					out_idx += 1
				r1_idx -= 1

			# sort the results
			for i in range(1, out_idx):
				r1 = self.out_r1[i]
				r2 = self.out_r2[i]
				op = self.out_op[i]
				rres = self.out_rres[i]
				rres_tol = self.out_tol[i]
				j = i - 1
				while j >= 0 and abs(self.out_tol[j]) > abs(rres_tol):
					self.out_r1[j + 1] = self.out_r1[j]
					self.out_r2[j + 1] = self.out_r2[j]
					self.out_op[j + 1] = self.out_op[j]
					self.out_rres[j + 1] = self.out_rres[j]
					self.out_tol[j + 1] = self.out_tol[j]
					j -= 1
				self.out_r1[j + 1] = r1
				self.out_r2[j + 1] = r2
				self.out_op[j + 1] = op
				self.out_rres[j + 1] = rres
				self.out_tol[j + 1] = rres_tol

		results = []
		for r1_idx in range(out_idx):
			out_prres = round((self.out_rres[r1_idx]) * 1000.0) / 1000.0
			
			out_prres = out_prres / scale
			out_vrres = self.roundWithMultiplier(out_prres, offset=4)[0]

			r1 = self.out_r1[r1_idx]
			r2 = self.out_r2[r1_idx]
   
			# offset = 4 - math.floor(math.log10(scale))
			offset = 4
			r1 = r1 / scale
			r2 = r2 / scale

			r1_text = self.roundWithMultiplier(r1, offset=offset)[0]
			r2_text = self.roundWithMultiplier(r2, offset=offset)[0]
			result = {
				"r1": r1_text,
				"op": self.out_op[r1_idx],
				"r2": r2_text,
				"rres": out_vrres,
				"tol": round(self.out_tol[r1_idx] * 100000.0) / 1000.0,
				"r1Val": r1,
				"r2Val": r2,
			}
			
			results.append(result)

			# order the results by tolerance (from best to worst)
			results = sorted(results, key=lambda k: abs(k['tol']))

		return results

	def prettyFormat(self, result):
		return f"{result['rres']} \t=\t{result['r1']} \t{result['op']}\t{result['r2']}\t(tol: {result['tol']}%)"




class PañoleroVirtual:
	def __init__(self):
		self.rf = ResistorFinder()
		self.res_bom = {}
		self.cap_bom = {}

		self.designator_components = {}

	def clear(self):
		self.res_bom = {}
		self.cap_bom = {}
		self.designator_components = {}
  
	def parse_params(self, input_string):
		# Define the multipliers
		multipliers = {
			'g': 1e9,
			'meg': 1e6,
			'k': 1e3,
			'': 1,
			'm': 1e-3,
			'u': 1e-6,
			'n': 1e-9,
			'p': 1e-12,
		}
		
		# remove all '=' from the input string
		input_string = input_string.replace('=', ' ')

		# Split the input string into parameters
		params = input_string.split()

		# Initialize the output list
		output = []

		# Initialize the last key
		last_key = None

		# Process each parameter
		for param in params:
			# Check if the parameter starts with a number
			if param[0].isdigit():
				# If it does, treat the last key and the current parameter as a key-value pair
				value = param

				# Check if the value has a suffix
				suffix = ''.join(filter(str.isalpha, value))
				value = value.replace(suffix, '')

				# If it does, remove the suffix and multiply the value by the corresponding multiplier
				if suffix.lower() in multipliers:
					value = float(value) * multipliers[suffix.lower()]
				else:
					# If it doesn't, just convert the value to a float
					value = float(value)

				# Add the key-value pair to the output list
				output.append((last_key, value))

				# Reset the last key
				last_key = None
			else:
				# If the parameter does not start with a number, save it as the last key
				last_key = param

		# Return the output list
		return output
  
	def add_resistor_to_bom(self, value, textValue):
		if self.res_bom.get(value) == None:
			self.res_bom[value] = {'count': 1, 'name': textValue}
		else:
			self.res_bom[value]['count'] += 1	# increment the value of the key

	def add_capacitor_to_bom(self, value, textValue):
		if self.cap_bom.get(value) == None:
			self.cap_bom[value] = {'count': 1, 'name': textValue}
		else:
			self.cap_bom[value]['count'] += 1	# increment the value of the key
  
	def add_components(self, text):
		# Parse the input string
		params = self.parse_params(text)
  
		components_added = False

		# Add the components to the designator_components dictionary
		for key, value in params:
      
			if key in self.designator_components:	# the key should not be in the dictionary
				print("The key", key, "is already in the dictionary")
				continue

			comboText = ""
			if key[0] == 'R':
				# combo = [r1, op, r2, rres, tol]
				self.rf.select_series(exp=[0, 7])
				combo = self.rf.compute(value, "parallel")[0]
				r1 = combo['r1Val']
				r2 = combo['r2Val']
				if r1 != 0:
					self.add_resistor_to_bom(r1, combo['r1'])
				if r2 != 0:
					self.add_resistor_to_bom(r2, combo['r2'])

				comboText = self.rf.prettyFormat(combo)

			elif key[0] == 'C':
				self.rf.select_series(exp=[-3, 12])

				combosLocos = self.rf.compute(value, "series", scale=1e12)
				combo = combosLocos[0]
				c1 = combo['r1Val']
				c2 = combo['r2Val']
				if c1 != 0:
					self.add_capacitor_to_bom(c1, combo['r1'])
				if c2 != 0:
					self.add_capacitor_to_bom(c2, combo['r2'])

				comboText = self.rf.prettyFormat(combo)

			self.designator_components[key] = comboText
			components_added = True
		return components_added

	def print_all(self):
		print("Resistor BOM:")
		print('Key'.ljust(15), 'Count'.ljust(15))

		self.res_bom = dict(sorted(self.res_bom.items()))

		for key, value in self.res_bom.items():
			count = value['count']
			name = value['name']
			print('R ' + str(name).ljust(15), str(count).ljust(5))
		print("\nCapacitor BOM:")
		print('Key'.ljust(15), 'Count'.ljust(15))

		self.cap_bom = dict(sorted(self.cap_bom.items()))

		for key, value in self.cap_bom.items():
			count = value['count']
			name = value['name']
			print('C ' + str(name).ljust(15), str(count).ljust(5))
		print("\nDesignator - Component:")
		print('Key'.ljust(15), 'Value'.ljust(15))

		for key, value in self.designator_components.items():
			print(str(key).ljust(20), str(value).ljust(20))

# Continuous program to compute resistor combinations

# Create the object
pibe = PañoleroVirtual()

while True:
	# Get the input
	text = input("Enter the components: ")

	# Check if the input is empty
	if text == "":
		# Print the BOMs
		pibe.print_all()
	elif text == "exit":
		break
	elif text == "clear":
		pibe.clear()
		print("\nBOMs cleared\n")
		continue
	elif text == "help":
		print("\nCommands:")
		print("  - clear: Clear the BOMs")
		print("  - exit: Exit the program")
		print("  - help: Show this help message")
		continue

	# Add the components to the BOM
	if pibe.add_components(text):
		pibe.print_all()
