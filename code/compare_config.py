

def main():

    f1 = 'Config.sh'
    f2 = 'Config.sh-CDM'

    with open(f1, 'r') as of:
        lines = of.readlines()

    keys1 = set()
    for line in lines:
        if line[0] == '#':
            continue
        if line.strip() == '':
            continue
        keys1.add(line.strip().split()[0])

    with open(f2, 'r') as of:
        lines = of.readlines()
    
    keys2 = set()
    for line in lines:
        if line[0] == '#':
            continue
        if line.strip() == '':
            continue
        keys2.add(line.strip().split()[0])


    diff = keys2 - keys1

    print(f"{f1} has:")
    print(*(keys1 - keys2), sep='\n')
    print()
    print(f"{f2} has:")
    print(*(keys2 - keys1), sep='\n')

    return

if __name__=="__main__":
    main()
