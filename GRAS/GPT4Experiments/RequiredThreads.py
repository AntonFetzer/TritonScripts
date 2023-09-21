

particles_initial = 2e+09
time_initial_hours = 3*24
threads_initial = 1
particles_new = 1.5e+12
time_new_hours = 3*24

# Calculate rates of processing
rate_initial = particles_initial / time_initial_hours * threads_initial  # particles per hour
rate_new = particles_new / time_new_hours  # particles per hour

# Calculate number of threads
threads_required = rate_new / rate_initial

print(f"You need {threads_required} threads.")