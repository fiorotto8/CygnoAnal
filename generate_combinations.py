import itertools
import random

param1 = list(range(100, 351, 50))  # start=100, stop=400, step=100
param2 = list(range(1, 5, 1))        # start=1, stop=3, step=1
param3 = list(range(1, 13, 2))       # start=3, stop=9, step=3
param4 = list(range(10, 51, 10))     # start=10, stop=50, step=20

combinations = list(itertools.product(param1, param2, param3, param4))
random.shuffle(combinations)

with open('params.txt', 'w') as f:
    for idx, p in enumerate(combinations, start=1):
        f.write(f"{idx} " + " ".join(map(str, p)) + "\n")