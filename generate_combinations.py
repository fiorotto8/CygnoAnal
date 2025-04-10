import itertools

param1 = list(range(100, 401, 100))  # start=100, stop=400, step=100
param2 = list(range(1, 4, 1))        # start=1, stop=3, step=1
param3 = list(range(3, 10, 3))       # start=3, stop=9, step=3
param4 = list(range(10, 51, 20))     # start=10, stop=50, step=20

with open('params.txt', 'w') as f:
    for p in itertools.product(param1, param2, param3, param4):
        f.write(" ".join(map(str, p)) + "\n")