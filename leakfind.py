#!/usr/bin/env python
# coding: utf-8

# # The purpose is to analyze Fortran code such that every allocation has its corresponding deallocation. Otherwise report the array that potentially gives memory leakage.

# In[1]:


import re
import glob


# In[2]:


def find_alloclist(fname, keyword, exclude=None):
    arrlist=[]
    with open(fname) as f:
        for line in f:
            if exclude != None:
                p=re.compile(exclude, re.IGNORECASE)
                m = p.search(line)
                if m:
                    continue

            p=re.compile(keyword, re.IGNORECASE)
            #print(p)
            m = p.search(line)
            if m:
                #print(line, m.end())
                for l in re.finditer("\s?\((.*)\)?", line[m.end():], re.IGNORECASE):
                    #print(l.group(1))
    
                    for n in re.finditer("([\w]+)\s?\(+.", l.group(1)):
                        #print(n.group(1))
                        arrlist.append(str(n.group(1)))
    #print(arrlist)
    return arrlist
                        
find_alloclist('allocate.dat', 'allocate','deallocate');


# In[3]:


def find_dealloclist(fname, keyword):
    arrlist=[]
    with open(fname) as f:
        for line in f:

            p=re.compile(keyword, re.IGNORECASE)
            #print(p)
            m = p.search(line)
            if m:
                #print(line, m.end())
                for l in re.finditer("\(?(\w+),?\)?", line[m.end():], re.IGNORECASE):
                    #print(l.group(1))
    
                    arrlist.append(str(l.group(1)))
    #print(arrlist)
    return arrlist
                        
find_dealloclist('allocate.dat', 'deallocate');


# In[4]:


ffiles=glob.glob("*.f")

of=open("leak.log", 'w+')
for f in ffiles:
    arrleak=[]

    alloclist=find_alloclist(f, 'allocate', 'deallocate')
    dealloclist=find_dealloclist(f, 'deallocate')
    
    #print("%s\n" %(str(f)), file=of)
    of.write("%s\n" %(str(f)))
    for arr in alloclist:
        if arr not in dealloclist:
            of.write("%s\t" %(str(arr)))
    
    of.write("\n\n")

            #print("%s\t" %(str(arr)), file=of)
            #arrleak.append(arr)

of.close()


# In[5]:


text="allOcate ( s0_s(sjl,sjkdlf), b(a))"
for m in re.finditer("allocate\s?\((.*)\)?", text, re.IGNORECASE):
    print(m.group(1))
    
    for n in re.finditer("([\w]+)\s?\([\w,]+\)", m.group(1)):
        print(n.group(1))
    #for arr in m.group(1):
    #    print(arr)

