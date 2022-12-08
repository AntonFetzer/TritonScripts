import pandas as pd
import argparse

parser = argparse.ArgumentParser()
# One positional and one optional argument
parser.add_argument('input', type=str, default="weather_tapiola.csv", help="The CSV file to load")
parser.add_argument('output', type=str, default="weather.png", help="The name of the result file")
parser.add_argument('-s', '--start', type=str, default="01/06/2021", help="Start date")
parser.add_argument('-e', '--end', type=str, default="01/10/2021", help="End date")

args = parser.parse_args()

#print(args.name + " was born on " + args.date)
print("The input file is ", args.input)
print("The output file is ", args.output)
print("The start date is ", args.start)
print("The end date is ", args.end)

input_file_name = args.input
output_file_name = args.output

# url = "https://raw.githubusercontent.com/AaltoSciComp/python-for-scicomp/master/resources/data/scripts/weather_tapiola.csv"
weather = pd.read_csv(input_file_name, comment='#')

# define the start and end time for the plot 
#start_date = pd.to_datetime('01/06/2021', dayfirst=True)
#end_date = pd.to_datetime('01/10/2021', dayfirst=True)
start_date = pd.to_datetime(args.start, dayfirst=True)
end_date = pd.to_datetime(args.end, dayfirst=True)

# The date format in the file is in a day-first format, which matplotlib does nto understand.
# so we need to convert it.
weather['Local time'] = pd.to_datetime(weather['Local time'], dayfirst=True)
# select the data
weather = weather[weather['Local time'].between(start_date, end_date)]

# Now, we have the data loaded, and adapted to our needs. So lets get plotting

# In[4]:


import matplotlib.pyplot as plt

# start the figure.
fig, ax = plt.subplots()
ax.plot(weather['Local time'], weather['T'])
# label the axes
ax.set_xlabel("Date of observation")
ax.set_ylabel("Temperature in Celsius")
ax.set_title("Temperature Observations")
# adjust the date labels, so that they look nicer
fig.autofmt_xdate()
# save the figure
fig.savefig(output_file_name)

# In[ ]:
