
def main():
    z = [ [sp.strip() for sp in l.strip().split(',')] for l in open('Zmat.csv') ]
    print(len(z))
    m = []
    for l in z[1:]:
        for r in l[1:]:
            if r not in m:
                m.append(r)
    print(m)


if __name__ == '__main__':
    main()

