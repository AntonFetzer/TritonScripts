import numpy as np
import matplotlib.pyplot as plt
from deap import base, creator, tools, algorithms


# Your Weibull function
def weibull_func(x, A, x0, W, s):
    result = np.zeros_like(x)
    # Remove all values of LET that are less than x0
    mask = x > x0
    # A * (1 - np.exp(-((x - x0) / W) ** s))
    result[mask] = A * (1 - np.exp(-((x[mask] - x0) / W) ** s))

    return result

# Your data
x_data = np.array([0.963, 2.246, 3.479, 6.123])   # LET [MeV cm2 mg-1]
y_data = np.array([9.031E-08, 2.632E-07, 1.447E-06, 2.254E-06])   # cross section [cm-2]
y_lower = np.array([9.1610E-09, 5.2084E-08, 4.7947E-07, 8.2306E-07])   # lower error [cm-2]
y_upper = np.array([3.2799E-07, 7.6510E-07, 2.8480E-06, 4.8890E-06])   # upper error [cm-2]

# x_data = np.array([0.00984463519313305, 0.00583690987124464, 0.00438369098712446, 0.00362845493562232])   # LET [MeV cm2 mg-1]
# y_data = np.array([2.37388724035608E-09, 1.84119677790564E-09, 2.14564369310793E-09, 2.68987341772152E-09])   # cross section [cm-2]
# y_lower = y_data * 0.9   # lower error [cm-2]
# y_upper = y_data * 1.1   # upper error [cm-2]

def chi_squared(individual):
    A, x0, W, s = individual
    
    # Ensure non-negative values for the Weibull function calculation
    if A <= 0 or W <= 0 or s <= 0 or (x_data <= x0).any():
        return np.inf,  # Return a large chi-squared value to penalize invalid parameters
    
    expected = weibull_func(x_data, A, x0, W, s)

    # Calculate the distances (errors) to the bounds
    errors_upper = y_upper - y_data  # Distances to the upper bound
    errors_lower = y_data - y_lower  # Distances to the lower bound
    
    # Calculate chi-squared using asymmetric errors
    chi2 = 0
    for i in range(len(y_data)):
        if expected[i] > y_data[i]:
            # If the model is above the data, use the distance to the upper bound
            error = errors_upper[i]
        else:
            # If the model is below the data, use the distance to the lower bound
            error = errors_lower[i]
        chi2 += ((y_data[i] - expected[i]) ** 2) / error**2

    return chi2,

# Genetic algorithm constants
POPULATION_SIZE = 500
P_CROSSOVER = 0.5  # probability for crossover
P_MUTATION = 0.01   # probability for mutating an individual
MAX_GENERATIONS = 25
HALL_OF_FAME_SIZE = 10

# Set the random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Create the Weibull fitting problem class
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

# Define a function to generate log-uniform distributed values
def loguniform(low=0, high=1, size=None):
    return np.exp(np.random.uniform(np.log(low), np.log(high), size))

toolbox = base.Toolbox()
toolbox.register("A", loguniform, max(y_lower), max(y_upper))
toolbox.register("x0", loguniform, min(x_data)/2, min(x_data))
toolbox.register("W", np.random.uniform, np.std(x_data)/2, np.std(x_data)*2)
toolbox.register("s", np.random.uniform, 0.5, 2)

# Structure initializers
toolbox.register("individualCreator", tools.initCycle, creator.Individual,
                 (toolbox.A, toolbox.x0, toolbox.W, toolbox.s), n=1)

toolbox.register("populationCreator", tools.initRepeat, list, toolbox.individualCreator)
toolbox.register("evaluate", chi_squared) # Fitness function
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("mate", tools.cxBlend, alpha=0.5)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.1)

# Genetic algorithm
def main():
    # Create initial population
    population = toolbox.populationCreator(n=POPULATION_SIZE)
    
    # Prepare the statistics object
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("min", np.min)
    
    # Define the hall-of-fame object
    hof = tools.HallOfFame(HALL_OF_FAME_SIZE)
    
    # Perform the Genetic Algorithm flow with hof feature added
    population, logbook = algorithms.eaSimple(population, toolbox,
                                              cxpb=P_CROSSOVER,
                                              mutpb=P_MUTATION,
                                              ngen=MAX_GENERATIONS,
                                              stats=stats,
                                              halloffame=hof,
                                              verbose=False)
    
    # Print info about the final solution
    best = hof.items[0]
    bestFitness = best.fitness.values[0]
    #print("-- Best Individual = ", best)
    #print("-- Best Fitness = ", bestFitness)
    
    # Extract the best fit parameters
    A, x0, W, s = best
    return A, x0, W, s, bestFitness


Results = {'A': [], 'x0': [], 'W': [], 's': [], 'Fitness': []}
NumberofRuns = 20
# Run the genetic algorithm 10 times
for i in range(NumberofRuns):
    print(f"Run {i+1}")

    # Run the genetic algorithm
    A, x0, W, s, Fitness = main()
    #print(f"Optimized parameters: A={A}, x0={x0}, W={W}, s={s}")
    Results['A'].append(A)
    Results['x0'].append(x0)
    Results['W'].append(W)
    Results['s'].append(s)
    Results['Fitness'].append(Fitness)


# Calculate the distances for the error bars
yerr_upper = y_upper - y_data  # Distance from the data point to the upper error
yerr_lower = y_data - y_lower  # Distance from the data point to the lower error


# Plot the data points
plt.errorbar(x_data, y_data, yerr=[yerr_lower, yerr_upper], fmt='.', label='Data', capsize=5, elinewidth=1, capthick=1)

# Generate x values for the Weibull fit
# 500 x values between 50% of the minimum x value and 150% of the maximum x value
x_Fit = np.geomspace(min(x_data)/10, max(x_data)*10, 500)

# Sort the resaults by fitness
Results = {key: [value for _, value in sorted(zip(Results['Fitness'], Results[key]))] for key in Results}

# Print the best Result
print(f"Best Result: A={Results['A'][0]:.3e}, x0={Results['x0'][0]:.3e}, W={Results['W'][0]:.3e}, s={Results['s'][0]:.3f}")

for i in range(NumberofRuns):
    y_Fit = weibull_func(x_Fit, Results['A'][i], Results['x0'][i], Results['W'][i], Results['s'][i])
    Fit_label = f'Fit {i+1}: Fitness={Results["Fitness"][i]:.3e}'
    # Plot the Weibull fit
    plt.plot(x_Fit, y_Fit, alpha=1/np.log(i+3))

plt.xlabel('LET [MeV cm2 mg-1]')
plt.ylabel('Cross Section [cm-2] ')
plt.title("SEE cross section vs LET")
#plt.xlim(min(x_data)*0.1, max(x_data)*2)
#plt.ylim(min(y_data)*0.1, max(y_data)*2)
plt.xscale('log')
plt.yscale('log')
plt.grid('both')
plt.legend()

plt.show()