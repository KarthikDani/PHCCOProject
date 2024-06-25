#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
import itertools as it


# In[2]:


file_path="C:\\Users\parig\OneDrive\Desktop\grn_sol2.csv"
grn=pd.read_csv(file_path,index_col=None)
grn


# In[3]:


grn.loc[:,'AMPK':'miR34']=(grn.loc[:,'AMPK':'miR34']-grn.loc[:,'AMPK':'miR34'].mean())/grn.loc[:,'AMPK':'miR34'].std()
grn


# In[4]:


plt.figure(figsize=(24,15))
sns.scatterplot(x=grn['AMPK'],y=grn['mTOR'])


# In[5]:


plt.figure(figsize=(24,15))
sns.scatterplot(x=grn['AMPK'],y=grn['ZEB'])


# In[12]:


g= grn.iloc[:,5]
g


# In[14]:


plt.figure(figsize=(24,15))
sns.regplot(x=grn['AMPK'],y=g)


# In[15]:


plt.figure(figsize=(20,10))
sns.regplot(x=grn['mTOR'],y=g)


# In[17]:


grn.describe()


# In[18]:


grn.loc[:,'fasting score']=grn.loc[:,'AMPK']-grn.loc[:,'mTOR']
grn


#  

#  
# 

#  
# 

# In[34]:


pip install seaborn --upgrade --user


# In[35]:


plt.figure(figsize=(20, 30))
#plt.rcParams['figure.dpi']= 200
sns.distplot(grn["fasting score"], bins=120, kde=True, color='purple', hist_kws={'alpha': 0.5, 'fill': True})

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Fasting Score', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')


# In[40]:


plt.figure(figsize=(10, 10))
#sns.set_style("whitegrid")
#plt.grid(False)

plo = sns.scatterplot(data=grn, x="fasting score", y=g, s = 2, color='Black' )

plt.xlabel("fasting score", fontweight='bold', size=14)
plt.ylabel('snail Expression',  fontweight='bold' , size=14 )
plt.xticks(size=12)
plt.yticks(size=12)


plo.plot(grid = False)
plt.axhline(y = -0.38, color = 'r', linestyle = ':')
plt.axhline(y = 0.38, color = 'r', linestyle = ':')
plt.axvline(x = -0.24, color = 'r', linestyle = ':')
plt.axvline(x = 0.64, color = 'r', linestyle = ':')


# In[41]:


plt.figure(figsize=(10, 10))
#sns.set_style("whitegrid")
#plt.grid(False)

plo = sns.scatterplot(data=grn, x="fasting score", y="miR200", s = 2, color='Black' )

plt.xlabel("fasting score", fontweight='bold', size=14)
plt.ylabel('miR200 Expression',  fontweight='bold' , size=14 )
plt.xticks(size=12)
plt.yticks(size=12)


plo.plot(grid = False)
plt.axhline(y = -0.38, color = 'r', linestyle = ':')
plt.axhline(y = 0.38, color = 'r', linestyle = ':')
plt.axvline(x = -0.24, color = 'r', linestyle = ':')
plt.axvline(x = 0.64, color = 'r', linestyle = ':')


# In[42]:


grn.loc[:,'emt score']=g+grn.loc[:,'ZEB']-grn.loc[:,'miR34']-grn.loc[:,'miR200']
grn


# In[43]:


plt.figure(figsize=(20, 30))
#plt.rcParams['figure.dpi']= 200
sns.distplot(grn["emt score"], bins=120, kde=True, color='red', hist_kws={'alpha': 0.5, 'fill': True})

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('emt Score', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')


# In[45]:


plt.figure(figsize=(8,7))
#sns.set_style("whitegrid")
#plt.grid(False)

plo = sns.scatterplot(data=grn, x="fasting score", y="emt score", s = 2, color='Black' )

plt.xlabel("fasting score", fontweight='bold', size=14)
plt.ylabel('emt score',  fontweight='bold' , size=14 )
plt.xticks(size=12)
plt.yticks(size=12)


plo.plot(grid = False)
plt.axhline(y = -0.38, color = 'r', linestyle = ':')
plt.axhline(y = 0.38, color = 'r', linestyle = ':')
plt.axvline(x = -0.24, color = 'r', linestyle = ':')
plt.axvline(x = 0.64, color = 'r', linestyle = ':')


# In[ ]:




