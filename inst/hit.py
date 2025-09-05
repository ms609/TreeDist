from collections import defaultdict
import re
import numpy as np

def HMI(Ut,Us, depth = 0):
    
    indent = "  " * depth
    print(f"{indent}HMI called with Ut={Ut}, Us={Us}")
    
    """
    This is from https://github.com/jipphysics/hit/blob/master/hit.ipynb
    
    Computes the hierarchical mutual information between two hierarchical partitions.
    
    Returns
    n_ts,HMI(Ut,Us) : where n_ts is the number of common elements between the hierarchical partitions Ut and Us.
    
    NOTE: We label by u,v the children of t,s respectively.
    
    Examples
    >>>"""
    
    Ut = str(Ut)
    Us = str(Us)
    
    Ut.replace(";", "")
    Us.replace(";", "")
    
    Ut = parse_nested(Ut)
    Us = parse_nested(Us)
    
    if isinstance(Ut[0],list):
        if isinstance(Us[0],list):
            # Ut and Us are both internal nodes since they contain other lists.
            n_ts=0.
            H_uv=0.
            H_us=0.
            H_tv=0.
            mean_I_ts=0.0
            n_tv=defaultdict(float)            
            for Uu in Ut:
                n_us=0.
                for v,Uv in enumerate(Us):
                    n_uv,I_uv=HMI(Uu,Uv, depth + 1)
                    print(f"{indent}  n_uv={n_uv}, I_uv={I_uv}")
                    n_ts+=n_uv
                    n_tv[v]+=n_uv
                    n_us+=n_uv                    
                    H_uv+=xlnx(n_uv)
                    mean_I_ts+=n_uv*I_uv
                H_us+=xlnx(n_us)
            for _n_tv in n_tv.values():
                H_tv+=xlnx(_n_tv)
            if n_ts>0.:
                local_I_ts=np.log(n_ts)-(H_us+H_tv-H_uv)/n_ts
                mean_I_ts=mean_I_ts/n_ts
                I_ts=local_I_ts+mean_I_ts
                return n_ts,I_ts
            else:
                return 0.,0.
        else:
            # Ut is internal node and Us is leaf
            return len(set(flattenator(Ut))&set(Us)),0.
    else:
        if isinstance(Us,list):
            # Ut is leaf and Us internal node
            return len(set(flattenator(Us))&set(Ut)),0.          
        else:
            # Both Ut and Us are leaves
            return len(set(Ut)&set(Us)),0.
        

def flattenator(newick):
    """Takes a hierarchical partition represented by nested lists and return a list of all its elements.
    
    Example
    >>> hp = [[3, 4, 5, 6], [[0], [1, 2]], [[7], [8, 9]]]
    >>> sorted(flattenator(hp))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    for e in newick:
        if isinstance(e,list):
            for ee in flattenator(e):
                yield ee
        else:
            yield e
            
def xlnx(x):
    """Returns x*log(x) for x > 0 or returns 0 otherwise."""
    if x <= 0.:
        return 0.
    return x*np.log(x)


def HH(hp):
    """Returns the hierarchical entropy of a hierarchical partition.
    
    Note: this is not the most efficient implementation."""
    return HMI(hp,hp)[1]

def HVI(hp1,hp2):
    """Returns the hierarchical variation of information."""
    return HH(hp1)+HH(hp2)-2.0*HMI(hp1,hp2)[1]

def mean_arit(x,y):
    return .5*(x+y)

def NHMI(hp1,hp2,generalized_mean=mean_arit):
    """Returns the normalized hierarchical mutual information.
    
    By default, it uses the arithmetic mean for normalization. However, another generalized mean can be provided if desired."""
    gm = generalized_mean(HH(hp1),HH(hp2))
    if gm > 0.:
        return HMI(hp1,hp2)[1]/gm
    return 0.

def removeCommas(line):
    newline = line
    removals = 0
    for i in range(len(line)-1):
        if line[i]==")" and line[i+1]==',':
            newline = newline[:i+1-removals] + newline[i+2-removals:]
            removals +=1
    return(str(newline))
    


def parse_nested(text, left=r'[(]', right=r'[)]', sep=r','):
    """Converts a newick string formated tree into a python nested list"""
    text = removeCommas(text)
    text = text.replace(" ", "")
    pat = r'({}|{}|{})'.format(left, right, sep)
    tokens = re.split(pat, text)    
    stack = [[]]
    for x in tokens:
        if not x or re.match(sep, x): continue
        if re.match(left, x):
            stack[-1].append([])
            stack.append(stack[-1][-1])
        elif re.match(right, x):
            stack.pop()
        else:
            stack[-1].append(x)
    return stack.pop()

def d_n(t1,t2,n=1):
    """Computes the distance metric associated to the HVI given by
        d_n(T,S)=1-exp(-n(ln(2)/2)V(T,S))
    """
    ln2d2=0.5*np.log(2.0)
    return 1.0-np.exp(-n*ln2d2*HVI(t1,t2))
