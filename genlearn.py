# Daniel Browne
# Lab 4 Genetic Algorithm Framework Used
# Most functions were modified for the purpose of this project

import random, string, numpy
import matplotlib.pyplot as plt

# Initialize the population with random individuals
def init(): 
   
	# Create an empty population
	pop = []
   
	for i in range(0, pop_size):
		genotype = bin(random.randint(0, 2**(genotype_size)))[2:]
		if len(genotype) < genotype_size:
			genotype = (genotype_size - len(genotype))*"0" + genotype
		pop.append(genotype)
	return(pop)

# Go through the population list, pick each genotype, evaluate it, and record its fitness score
def evaluate(pop):
	# Create an empty list of fitness values
	fitness = []
	for i in range(pop_size):
		genotype = pop[i]
		temp_fit = 0
		# Fitness is described by how long it takes to get to a stable correct output 
		# for each set of inputs and initial states
		for j in input_list:
			exp = str(int(j.count('1') >= j.count('0')))
			total_count = 0
			true_count = 0
			for k in range(2**(internal_nodes+1)):
				rand_part = bin(k)[2:]
				if len(rand_part) < internal_nodes+1:
					rand_part = (internal_nodes + 1 - len(rand_part))*"0" + rand_part
				network = str(j) + str(rand_part)
				while total_count < total_length and true_count != true_length:
					network = step(network,genotype)
					if network[-1] == exp:
						true_count = true_count + 1
					else:
						true_count = 0
					total_count = total_count + 1
				if true_count < true_length:
					true_count = total_count
				temp_fit = temp_fit + true_count

		fitness.append(temp_fit)
	
	# Normalize fitness. IMPORTANT: If you don't multiple by 1.0, you will get all zeros
	fitness = [f * 1.0/sum(fitness) for f in fitness]
	
	return(fitness)

# Using one genotype, evaluate it, and record its fitness score
def single_eval(genotype):
	# Create an empty list of fitness values
	fitness = 0
	# Fitness is described by how long it takes to get to a stable correct output 
	# for each set of inputs and initial states
	for j in input_list:
		exp = str(int(j.count('1') >= j.count('0')))
		total_count = 0
		true_count = 0
		for k in range(2**(internal_nodes+1)):
			rand_part = bin(k)[2:]
			if len(rand_part) < internal_nodes+1:
				rand_part = (internal_nodes + 1 - len(rand_part))*"0" + rand_part
			network = str(j) + str(rand_part)
			while total_count < total_length and true_count != true_length:
				network = step(network,genotype)
				if network[-1] == exp:
					true_count = true_count + 1
				else:
					true_count = 0
				total_count = total_count + 1
			if true_count < true_length:
				true_count = total_count
			fitness = fitness + true_count
	
	return(fitness)

# Divide each genotype into its individual genes that represent each internal/output node
def subgene(genotype):
	genotype = str(genotype)
	gene_size = 3 + gene_nodes
	result = []
	for i in range(1,internal_nodes+2):
		result.append(str(genotype[gene_size*(i-1):gene_size*i]))
	return result

# Generate the binary network at the next time step based on what the genotype dictates
def step(network, genotype):
	new_network = network
	split = subgene(genotype)
	for i in range(input_size,gene_nodes+1):
		gene = split[i-input_size]
		view = []
		for j in range(3,len(gene)):
			if gene[j] == '1':
				view.append(network[j-3])
		match gene[0:3]:
			case '000': #0
				new_num = '0'
			case '001': #or
				if '1' in view:
					new_num = '1'
				else:
					new_num = '0'
			case '010': #xor
				if view.count('1') == 1:
					new_num = '1'
				else:
					new_num = '0'
			case '011': #and
				if not '0' in view:
					new_num = '1'
				else:
					new_num = '0'
			case '100': #1
				new_num = '1'
			case '101': #nor
				if not '1' in view:
					new_num = '1'
				else:
					new_num = '0'
			case '110': #xnor
				if view.count('1') != 1:
					new_num = '1'
				else:
					new_num = '0'
			case _: #nand
				if view.count('0') >= 1:
					new_num = '1'
				else:
					new_num = '0'
		new_network = new_network[0 : i] + new_num +  new_network[i+1:]
	return new_network

# Select two parents using roulette wheel method, also known as 'fitness proportionate' selection
# Reference: http://en.wikipedia.org/wiki/Fitness_proportionate_selection
def select(fitness, pop):

	# Favor lower fitness values
	fitness = [1.0/f for f in fitness]
	fitness = [f * 1.0/sum(fitness) for f in fitness]

	# Select first parent
	parent_1_index = -9999  #indicates that parent_1 is yet to be chosen
	cumulated_fitness = 0
	roulette_marker = random.random()  #indicates that the 'roulette wheel' has settled on a random number
	for i in range(pop_size):
		cumulated_fitness += fitness[i]
		if cumulated_fitness >= roulette_marker:
			parent_1_index = i
			break
	
	# Select second parent different from the first parent
	parent_2_index = parent_1_index  #indicates that parent_2 is yet to be chosen
	while parent_2_index == parent_1_index:  #this ensures that the two parents chosen are distinct
		cumulated_fitness = 0
		roulette_marker = random.random()  #indicates that the 'roulette wheel' has settled on a random number
		for i in range(pop_size):
			cumulated_fitness += fitness[i]
			if cumulated_fitness >= roulette_marker:
				parent_2_index = i
				break

	return([pop[parent_1_index], pop[parent_2_index]])


# Recombine two parents to produce two offsprings, with a certain probability
# specified by recombination_rate
def recombine(parent_1, parent_2, recombination_rate):
   
	r = random.random()
   
	if r <= recombination_rate:  #recombination_rate is a value between 0 and 1
		#recombine
		slice_num = random.randint(1, 2*(internal_nodes+1))
		if slice_num % 2 == 1:
			slice_point = 3 + int(slice_num/2)*(3+gene_nodes)
		else:
			slice_point = int(slice_num/2)*(3+gene_nodes)
		offspring_1 = parent_1[0 : slice_point] + parent_2[slice_point : genotype_size]
		offspring_2 = parent_2[0 : slice_point] + parent_1[slice_point : genotype_size]
		return([offspring_1, offspring_2])
   
	else:  #don't recombine
		return([parent_1, parent_2])

# Mutate a genotype with a certain probability of mutation per-bit specified by 'mutation_rate'.
def mutate(genotype, mutation_rate):
   
	mutated_genotype = genotype  #indicates that the genotype is yet to be mutated
   
	for i in range(genotype_size):
   
		r = random.random()
   
		if r <= mutation_rate:
			if (genotype[i]) == "1": new_bit = "0"
			else: new_bit = "1"

			mutated_genotype = mutated_genotype[0 : i] + new_bit + mutated_genotype[(i+1) : genotype_size]

	return(mutated_genotype)

# GA parameers
pop_size = 100  #population size
num_generations = 100
num_experiments = 10
input_size = 3
internal_nodes = 5
recombination_rate = 0.2
mutation_rate = 0.1
elite_rate = 0.05

true_length = 5
total_length = 10

gene_nodes = input_size + internal_nodes
genotype_size = (internal_nodes+1)*(gene_nodes+3)

min_score = 5 * 2**(gene_nodes+1)
graph_x = []

final_y = []
input_list = []
for i in range(num_generations):
	graph_x.append(i)
	final_y.append(0)
for i in range(2**input_size):
	new_in = bin(i)[2:]
	if len(new_in) < input_size:
		new_in = (input_size - len(new_in))*"0" + new_in
	input_list.append(new_in)

for n in range(num_experiments):
	# Main GA loop.

	pop = init()

	print("Initial population generated.")
	graph_y = []
	for gen in range(num_generations):
		fitness = evaluate(pop)
		best = numpy.argmin(fitness)
		percent = ((min_score/true_length)-(single_eval(pop[best])-min_score)/(total_length-true_length))/(min_score/true_length)
		print("Generation",gen,"best genotype:",pop[best], ", with a success rate of:",percent)

		#graph_x.append(gen)
		graph_y.append(percent)

		if percent == 1.0:
			print("Found a genotype that has a 100% successs rate in generation" , gen,)
			for i in range(gen+1,num_generations):
				#graph_x.append(i)
				graph_y.append(1.0)
			break
		
		# create a new population that will hold the next generation of genotypes
		# if you are using 'elite' selection, then simply copy the top x% of individuals from pop to new_pop
		temp_fitness = fitness.copy()
		new_pop = []
		for i in range(int(elite_rate*pop_size)):
			min_i = numpy.argmin(temp_fitness)
			new_pop.append(pop[min_i])
			temp_fitness[min_i] = -1

		while (len(new_pop) < pop_size):  # continue loop until new_pop has pop_size number of individuals

			[parent_1, parent_2] = select(fitness, pop)

			[offspring_1, offspring_2] = recombine(parent_1, parent_2, recombination_rate)

			mutated_genotype_1 = mutate(offspring_1, mutation_rate)
			mutated_genotype_2 = mutate(offspring_2, mutation_rate)
			
			new_pop.append(mutated_genotype_1)
			new_pop.append(mutated_genotype_2)
			
			# continue loop until new_pop has pop_size number of individuals

		pop = new_pop  #replace current population with new population


	final = evaluate(pop)
	best = numpy.argmin(final)

	best_fit = single_eval(pop[best])
	print("Best genotype:", subgene(pop[best]))
	print("This genotype scored a",best_fit, "out of a minimum of", min_score,". Or a", percent, "correct detection rate.")

	for i in range(len(graph_y)):
		final_y[i] = final_y[i] + graph_y[i]

final_y = [y * 1.0 / num_experiments for y in final_y]

plt.plot(graph_x,final_y)
plt.xlabel('# of Generations')
plt.ylabel('Correct Detection Rate')
title = "Input Nodes = " + str(input_size) + " Internal Nodes = " + str(internal_nodes) + " Averaged Accross " +  str(num_experiments) +  " Experiments"
plt.title(title)
plt.ylim(0,1)
plt.show()