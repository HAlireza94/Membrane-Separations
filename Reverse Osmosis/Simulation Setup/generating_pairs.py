import numpy as np
c = 0
symb, e, sig = [], [], []
with open('potential') as f:
    for line in f:
        p = line.split()
        c += 1
        if c != 1:
            symb.append(p[0])
            e.append(float(p[1]))
            sig.append(float(p[2]))

pairs,e_ij,sig_ij = [], [], []
for i in range(len(sig)):
    for j in range(i+1,len(sig)):
        pairs.append(str(symb[i])+"-"+str(symb[j]))
        e_ij.append(np.sqrt(e[i]*e[j]))
        sig_ij.append(0.5*(sig[i]+sig[j]))

print(f"{'Pairs':<15}{'ε (kcal/mol)':<20}{'σ (Å)':<15}")

symb1 = []
for i in range(len(symb)):
    symb1.append(symb[i]+"-"+symb[i])

for i in range(len(symb1)):
    print(f"{symb1[i]:<15}{e[i]:<15.8f}{sig[i]:<15.8f}")
print("")
for i in range(len(pairs)):
    print(f"{pairs[i]:<15}{e_ij[i]:<15.8f}{sig_ij[i]:<15.8f}")
