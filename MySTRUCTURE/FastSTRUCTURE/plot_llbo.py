import matplotlib.pyplot as plt


LLBOS = []
with open('log.txt', 'r') as f:
    all_lines = f.readlines()
    for l in all_lines:
        if l.find('itr #') > 0:
            idx = l.index('LLBO:') + 5
            ll = float(l[idx:])
            LLBOS.append(ll)

x = range(1, len(LLBOS) + 1)
plt.plot(x, LLBOS)
plt.show()

