#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import h5py
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
import pandas as pd
import allel; print('scikit-allel', allel.__version__)
from skbio import DistanceMatrix
import allel.chunked
allel.chunked.storage_registry['default'] = allel.chunked.hdf5mem_zlib1_storage
import os


# In[4]:


for chro in list(range(1,30))+['Y','X','MT']:
    command="vcftools --gzvcf /halde/guilherme/Run9/filtered/Chr"+str(chro)+".recode.vcf.gz --recode --keep /halde/guilherme/samples_1kbulls.txt --stdout | bgzip -c > /halde/guilherme/Run9/paper/Chr"+str(chro)+".vcf.gz &"
    os.system(command)


# In[2]:


breeds={'DSN':[]}
bulls_project=pd.read_excel("/home/guilherme/NAS/DNA/1000bulls/Run9-TAUraw-animalDistributionList-20210630.xlsx")
bulls_project.index=bulls_project['SampleID']
callset = allel.read_vcf("/halde/guilherme/Run9/paper/Chr1.vcf.gz",alt_number=1)
for i,sample in enumerate(callset['samples'][:]):
    try:
        sample_breed = bulls_project.loc[sample,'Breed']
        if sample_breed=='FjällCattle':
            sample_breed='Fjäll'
        if sample_breed=='TraditionalAngler':
            sample_breed='GermanRedAngler'
        if sample_breed=='HolsteinFriesian':
            sample_breed='Holstein'
        if sample_breed=='MeuseRhineYssel' or sample_breed=='Rotbunt' :
            sample_breed='RedWhiteDualPurpose'
        if  sample_breed not in list(breeds.keys())+['DeutschesSchwarzbuntesNiederungsrind']:
            breeds[sample_breed]=[]
        if sample_breed=='DeutschesSchwarzbuntesNiederungsrind':
            breeds['DSN'].append(i)
        else:
            breeds[sample_breed].append(i)
    except:
        print(sample)
for key in breeds.keys():
    print(key,len(breeds[key]))


# In[3]:


chroms=list(callset['variants/CHROM'])
positions=list(callset['variants/POS'])
gt = allel.GenotypeChunkedArray(callset['calldata/GT'])

for chro in list(range(2,30)):    
    callset = allel.read_vcf("/halde/guilherme/Run9/paper/Chr"+str(chro)+".vcf.gz",alt_number=1)
    gt = gt.concatenate([allel.GenotypeChunkedArray(callset['calldata/GT'])], axis=0)
    chroms+=list(callset['variants/CHROM'])
    positions+=list(callset['variants/POS'])


# In[37]:


chroms=[]
positions=[]
#target_breeds=['DSN','Angus','GroningenWhiteHeaded','DutchBelted','Holstein','HolsteinRed','ModernAngler','TraditionalAngler','GermanRedAngler','RotesHöhenvieh', 'RedWhiteDualPurpose','Holstein_DNK','Holstein_NLD','Holstein_USA','DutchFriesianRed','Shorthorn','DeepRedCattle','BelgianRedWhiteCampine','DutchImprovedRed','NorwegianRed', 'AyrshireFinnish','HolsteinFriesian','Fjäll','Kholmogory','NorthernFinncattle','ModernDanishRed','SwedishRed','TraditionalDanishRed']
target_breeds=[breed for breed in breeds.keys() if len(breeds[breed])>=10]
diversity={1:{}}
gt = allel.GenotypeChunkedArray(callset['calldata/GT'])
acs = gt.count_alleles_subpops({key:breeds[key] for key in target_breeds})
for breed in target_breeds:
    seg=acs[breed].is_segregating()
    pos=callset['variants/POS'].take(np.where(seg))[0]
    pi, windows, n_bases, counts = allel.windowed_diversity(pos, acs[breed].compress(seg,axis=0), size=10000, start=pos[0], stop=pos[-1])
    diversity[1][breed]=pi.mean()
chroms+=list(callset['variants/CHROM'])
positions+=list(callset['variants/POS'])
print(diversity[1])
    
for chro in list(range(2,30))+['X','MT','Y']:    
    callset = allel.read_vcf("/halde/guilherme/Run9/paper/Chr"+str(chro)+".vcf.gz",alt_number=1)
    gt = gt.concatenate([allel.GenotypeChunkedArray(callset['calldata/GT'])], axis=0)
    acs = allel.GenotypeChunkedArray(callset['calldata/GT']).count_alleles_subpops(breeds)
    diversity[chro]={}
    chroms+=list(callset['variants/CHROM'])
    positions+=list(callset['variants/POS'])
    for breed in target_breeds:
        seg=acs[breed].is_segregating()
        pos=callset['variants/POS'].take(np.where(seg))[0]
        pi, windows, n_bases, counts = allel.windowed_diversity(pos, acs[breed].compress(seg,axis=0), size=10000, start=pos[0], stop=pos[-1])
        diversity[chro][breed]=pi.mean()
        
    print(diversity[chro])
gt


# In[39]:


for breed in target_breeds:
    lista=[]
    for chro in list(range(1,30)):
        lista.append(diversity[chro][breed])
    print(breed,' ',np.mean(lista))


# In[40]:


import json
a_file = open("/halde/guilherme/diversity/nucleotide_diversity.json", "w")
json.dump(diversity, a_file)


# In[47]:


diversity_df.to_excel('/halde/guilherme/diversity/pi.xlsx',index=False)


# In[46]:


diversity_df=pd.DataFrame([item.values() for item in list(diversity.values())],columns=list(list(diversity.values())[0].keys()))
import matplotlib as mpl
mpl.rc('axes',labelsize=25)
mpl.rc('xtick', labelsize=25) 
mpl.rc('ytick', labelsize=25)
mpl.rc('legend',fontsize=25)
f, ax = plt.subplots(figsize=(21, 40))
sns.violinplot(data=diversity_df, palette="Set3", bw=.2, cut=1, linewidth=1,orient="h")
sns.despine()
ax.set_xlabel('Nucleotide Diversity (/bp)')


# In[4]:


df_fsts=pd.DataFrame()
df_fsts['CHROM']=chroms
df_fsts['POS']=positions


# In[5]:


#segregating in DSN
target_breeds=[breed for breed in breeds.keys() if len(breeds[breed])>=10]
acs = gt.count_alleles_subpops({breed:breeds[breed] for breed in target_breeds})
seg=acs['DSN'].is_segregating()
#get positions for DSN and Holstein
df_fsts['DSN_seg']=seg


# In[6]:


len(acs['DSN'].compress(seg,axis=0))


# In[7]:


fst=[]
segs=[]
ns=[]
for breed in target_breeds:
    seg=acs[breed].compress(df_fsts['DSN_seg'],axis=0).is_segregating()
    num, dem = allel.hudson_fst(acs['DSN'].compress(df_fsts['DSN_seg'],axis=0).compress(seg,axis=0), acs[breed].compress(df_fsts['DSN_seg'],axis=0).compress(seg,axis=0))
    fst.append(np.nansum(num)/np.nansum(dem))
    segs.append(len(acs[breed].compress(df_fsts['DSN_seg'],axis=0).compress(seg,axis=0)))
    ns.append(len(breeds[breed]))
    if breed=='Holstein':
        df_fsts['HOL_seg']=acs['Holstein'].is_segregating()
        df_fsts2=df_fsts[df_fsts['DSN_seg']&df_fsts['HOL_seg']]
        df_fsts2['fst'] = num / dem
        df_fsts2['fst_abs']=df_fsts2['fst'].abs()
        df_fsts2=df_fsts2.reset_index(drop=True)
        df_fsts2['ind']=df_fsts2.index
        df_fsts2.to_csv("/halde/guilherme/fsts_DSN_HOL_raw2.tab",sep='\t',index=False,columns=['CHROM','POS','fst','fst_abs'])
fst_df=pd.DataFrame()
fst_df['Fst']=fst
fst_df['Breed']=target_breeds
fst_df['Polymorfic']=segs
fst_df['N']=ns
fst_df=fst_df.sort_values(by=['Fst'])
fst_df


# In[8]:


fst_df.to_excel('/halde/guilherme/fsts_run9_raw.xlsx',index=False)


# In[9]:


df_fsts2.describe()


# In[10]:


df_fsts2['fst_abs'].mean()+2*df_fsts2['fst_abs'].std()


# In[11]:


df_fsts2['fst_abs'].mean()+3*df_fsts2['fst_abs'].std()


# In[55]:


import matplotlib as mpl
mpl.rc('axes',labelsize=25)
mpl.rc('xtick', labelsize=25) 
mpl.rc('ytick', labelsize=25)
mpl.rc('legend',fontsize=25)
fig = plt.figure()
ax = fig.add_subplot(111)
x_labels = []
x_labels_pos = []
df_grouped = df_fsts.groupby(('CHROM'))
for num, (name, group) in enumerate(df_grouped):
    try:        
        if int(name)%2==0:
            indice=0
        else:
            indice=1
    except:
        indice=0
    group.plot(kind='scatter', x='ind', y='fst_abs', ax=ax,color=sns.color_palette("Set2")[indice])
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
#ax.set_xticklabels(x_labels)
ax.set_xticklabels(['1','10','11','12','13','14','15','16','17','18','19','2','20','21', '22','23','24','25','26','27','28','29','3','4','5','6','7','8', '9', '','X MT Y',''])
ax.set_xlim([0, len(df_fsts)])
ax.set_ylim([0, 1])
ax.set_xlabel('Chromosome')
ax.set_ylabel('Fst')
sns.despine()
ax.plot([0, len(df_fsts)], [0.7, 0.7], 'red', lw=2)
fig.set_size_inches(45, 15)


# In[10]:


df_fsts2[df_fsts2['fst_abs']>0.7]


# In[11]:


table=pd.read_csv("/home/guilherme/NAS/DNA/SNPchip/DSN_SNP-Chip/dsn_selection-complete_table.tab",sep='\t')
table=table.astype({'chromosome':str,'POS':int})
table.columns


# In[ ]:


##Annotation from the 1kbulls: "/home/guilherme/NAS/DNA/1000bulls/Run9/Tau/annotations/Chr1-Run9-TAU-beagle-toDistribute.vcf.tab.gz"


# In[20]:


df_fsts_high=df_fsts2[df_fsts2['fst_abs']>0.7]
df_fsts_high=df_fsts_high.astype({'CHROM':str,'POS':int})
df_fsts_high=df_fsts_high.reset_index(drop=True)
df_fsts_high=pd.merge(df_fsts_high,table[['chromosome','POS','annotation']],how='left',left_on=['CHROM','POS'],right_on=['chromosome','POS'])
df_fsts_high


# In[13]:


#get list of genes and impact
genes=[]
impact=[]
for i in df_fsts_high.index:
    print(df_fsts_high.loc[i,'annotation'])
    if 'LOW' in str(df_fsts_high.loc[i,'annotation'])or 'MODERATE' in str(df_fsts_high.loc[i,'annotation']) or 'HIGH' in str(df_fsts_high.loc[i,'annotation']):
        impact.append(True)
    else:
        impact.append(False)
    try:
        for item in df_fsts_high.loc[i,'annotation'].split('|'):
            if 'ENSBTAG' in item and item not in genes:
                genes.append(item)
    except:
        continue
len(genes)


# In[14]:


genes


# In[ ]:


fst=[]
for breed in breeds.keys():
    a, b, c = allel.weir_cockerham_fst(gt, [breeds['DSN'],breeds[breed]])
    fst.append(np.nansum(a) / (np.nansum(a) + np.nansum(b) + np.nansum(c)))
fst_df=pd.DataFrame()
fst_df['Fst']=fst
fst_df['Breed']=list(breeds.keys())[:]
fst_df=fst_df.sort_values(by=['Fst'])
fst_df


# In[6]:


import gc
del gt
gc.collect()


# In[77]:


import os
for chro in list(range(1,30))+['X-R1','X-R2','MT']:
    command="vcftools --gzvcf /halde/guilherme/Run9/filtered/Chr"+str(chro)+".recode.vcf.gz --het --keep /halde/guilherme/samples_1kbulls.txt --out /halde/guilherme/diversity/CHR"+str(chro)+" &"
    os.system(command)


# In[28]:


#Heterosigosity in DSN
het=pd.read_csv("/halde/guilherme/dsn_het_vcftools.het",sep='\t')
het['O(HET)']=(het['N_SITES'] - het['O(HOM)']) / het['N_SITES']
het['E(HET)']=(het['N_SITES'] - het['E(HOM)']) / het['N_SITES']
het


# In[29]:


het.describe()


# In[30]:


from scipy.stats import chisquare
chisquare(f_obs=list(het['O(HET)']),f_exp=list(het['E(HET)']))


# In[ ]:




