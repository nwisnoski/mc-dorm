# This function sets the parameters to be used in the following scripts:
# `create-landscape.R`,

# Set the dimensions of the landscape (l by l square)
# and the number of patches, m
m = 20
l = 100

# initialize landscape with s species, with prefs for e niche dimensions
s = 10
e = 2

# probability of dispersal
disp = .1

# importance of the environment (from 0 (netural) to 1 (niche))
env.importance = 0

# vital rates
d = 0.1 # probability of death
b = 0.3 # probability of birth


