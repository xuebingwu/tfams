import pandas as pd
import numpy as np
import os
from params import output_dir

# for plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def hamming(s1,s2): return sum(a!=b for a,b in zip(s1,s2))    

def get_codon_table():
	return dict(zip(codons, amino_acids))
	
def get_inverted_codon_table():
	ct = get_codon_table()
	inv_codon_table = {}
	for k, v in ct.items():
		inv_codon_table[v] = inv_codon_table.get(v, [])
		inv_codon_table[v].append(k)
	return inv_codon_table

def prepare_count_matrix(df):
    matrix = pd.DataFrame(data = 0, index = codons, columns=list('ACDEFGHKLMNPQRSTVWY'),dtype=float)
    df = df[df['DP Decoy']!='+']
    df = df[pd.notnull(df['codon'])]
    df['codon'] = df['codon'].map(lambda x: x.replace('T','U'))
    for label in matrix.index:
        if codon_table[label] == '*':
            matrix.loc[label] = float('NaN')
        for col in matrix.columns:
            if (label in inverted_codon_table[col]) or (codon_table[label] +' to '+col in exact_PTM_spec_list):
                matrix.loc[label, col] = float('NaN')
    subs_agg = pd.DataFrame(np.array(zip(*df.groupby(['protein','position','origin','destination','codon']).groups.keys())).T, columns=['protein','position','origin','destination','codon'])												
    for x, l in subs_agg.groupby(['codon', 'destination']).groups.items():
        codon, destination = x
        if (codon in matrix.index) and pd.notnull(matrix.loc[codon,destination]):
            matrix.loc[codon,destination] = len(l)
    matrix.rename(columns={"L": "I/L"},inplace=True)
    return matrix

def prepare_count_matrix2(df):
    matrix = pd.DataFrame(data = 0, index = codons, columns=list('ACDEFGHKLMNPQRSTVWY'),dtype=float)
    df = df[df['DP Decoy']!='+']
    df = df[pd.notnull(df['codon'])]
    df['codon'] = df['codon'].map(lambda x: x.replace('T','U'))
    for label in matrix.index:
        if codon_table[label] == '*':
            matrix.loc[label] = float('NaN')
        for col in matrix.columns:
            if (label in inverted_codon_table[col]) or (codon_table[label] +' to '+col in exact_PTM_spec_list):
                matrix.loc[label, col] = float('NaN')
    counts = df[['codon','destination','Charge']].groupby(['codon','destination']).count()	
    for codon,destination in counts.index:
        matrix.loc[codon,destination] = counts.Charge[codon][destination]
    matrix.rename(columns={"L": "I/L"},inplace=True)
    return matrix
	
def probe_mismatch(codon1, codon2, pos, spec):
    origin, destination = spec
    for i in range(3):
        if i == pos:
            if codon1[i] != origin or codon2[i] != destination:
                return False
        else:
            if codon1[i] != codon2[i]:
                return False
    return True

def plot_matrix(m,filename):
    fig, ax = plt.subplots()
    im = ax.imshow(m,cmap='magma')

    ax.set_facecolor("gray")

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(m.columns)))
    ax.set_yticks(np.arange(len(m.index)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(m.columns,fontsize=5)
    ax.set_yticklabels(m.index,fontsize=5)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    '''
    for i in range(len(index)):
        for j in range(len(column)):
            text = ax.text(j, i, m[i, j],ha="center", va="center", color="w")
    '''

    #ax.set_title("Title")
    fig.tight_layout()
    #plt.show()
    plt.savefig(filename)  

    '''
    # example
    index = ["AAA", "CCC", "UUU", "GGG", "CGA", "GUG", "UGG"]
    columns = list('ABCDEFG')
    d = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
                            [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
                            [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
                            [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
                            [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
                            [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1],
                            [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])
    m = pd.DataFrame(data = d, index = index, columns=columns,dtype=float)
    plot_matrix(m,'test.pdf')
    '''

bases = 'UCAG'
codons = [a+b+c for a in bases for b in bases for c in bases]

amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
RC = {'A':'U', 'C':'G', 'G':'C', 'U':'A'}

codon_table = get_codon_table()
inverted_codon_table = get_inverted_codon_table()
inverted_codon_table['L'] = inverted_codon_table['L'] + inverted_codon_table['I']
tol = 0.005
MW_dict = {"G": 57.02147, 
            "A" : 71.03712, 
            "S" : 87.03203, 
            "P" : 97.05277, 
            "V" : 99.06842, 
            "T" : 101.04768, 
            "I" : 113.08407, 
            "L" : 113.08407, 
            "N" : 114.04293, 
            "D" : 115.02695, 
            "Q" : 128.05858, 
            "K" : 128.09497, 
            "E" : 129.0426, 
            "M" : 131.04049, 
            "H" : 137.05891,
            "F" : 147.06842, 
            "R" : 156.10112, 
            "C" : 160.030654, #CamCys
            "Y" : 163.0633,
            "W" : 186.07932,
            }

aas_sorted_by_mass = [i[0] for i in sorted(MW_dict.items(),key=lambda x:x[1])]
danger_mods = pd.read_pickle('danger_mods')
exact_PTM_spec = pd.DataFrame(index = aas_sorted_by_mass,
                              columns = aas_sorted_by_mass,
                              dtype = int)

for aa1 in MW_dict.keys():
	for aa2 in MW_dict.keys():
		delta_m = MW_dict[aa2] - MW_dict[aa1]
		exact_PTM_spec.loc[aa1,aa2]=len(danger_mods[(danger_mods['delta_m']<delta_m + 0.0005) & (danger_mods['delta_m']>delta_m - 0.0005) & (danger_mods['site']==aa1)]) > 0

exact_PTM_spec_list = [str(i) + ' to ' + str(j) for i in aas_sorted_by_mass for j in  aas_sorted_by_mass if exact_PTM_spec.loc[i,j]] 

mask = pd.DataFrame(data = False,
                    index = codons,
                    columns = list('ACDEFGHKLMNPQRSTVWY'),
                    dtype = float)	

for label in codons:
	near_cognates = np.array([hamming(i,label)==1 for i in codons])
	reachable_aa = set(np.array(list(amino_acids))[near_cognates])
	mask.loc[label] =[i in reachable_aa for i in 'ACDEFGHKLMNPQRSTVWY']
	
for label in mask.index:
	if codon_table[label] == '*':
		mask.loc[label]=float('NaN')
	for col in mask.columns:
		if (label in inverted_codon_table[col]) or (codon_table[label] +' to '+col in exact_PTM_spec_list):
			mask.loc[label, col] = float('NaN')

subs = pd.read_pickle(os.path.join(output_dir,'subs'))
subs = subs[~subs.decoy]

data = prepare_count_matrix2(subs)

plot_matrix(data,output_dir+'/subs-count-plot.pdf')

data.to_csv(output_dir+'/subs-count.txt', header=True, index=True, sep='\t', float_format='%d')

data.to_pickle(output_dir+'/unique_substitutions_count_matrix')

log_data = np.log2(data+1)
plot_matrix(log_data,output_dir+'/subs-log-count-plot.pdf')

print('- substitution count table: '+output_dir+'/subs-count.txt')
print('- heatmap: '+output_dir+'/subs-log-count-plot.pdf')

#m = log_data.max().max()


