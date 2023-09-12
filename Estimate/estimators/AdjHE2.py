import numpy as np

def AdjHE2(A,y,trA=None,trA2=None, npc= 0):
    y = np.array(y)
    y = y.reshape((y.shape[0],1))
    std_y = (y-np.mean(y))/np.std(y)
    if (trA is None) and (trA2 is None):
        trA = np.sum(np.diag(A))
        trA2 = np.sum(np.multiply(A,A))
    n = A.shape[1]
    yay = np.dot(std_y.T,np.dot(A,std_y))
    yty = np.dot(std_y.T,std_y)
    if (npc==0):
        denominator = trA2 - 2*trA + n
        nominator = n - trA + yay - yty
    else:
        pc = final_PC
        s = np.diag(np.dot(pc.T,np.dot(A,pc)))
        b = s - 1
        c = np.dot(std_y,pc)**2 - 1
        denominator = trA2 - 2*trA + n - np.sum(b**2)
        nominator = n - trA + yay - yty - np.sum(b*c)
    nominator = abs(nominator)
    denominator = abs(denominator)
    h2 = nominator/denominator
    var_ge = 2/denominator
    print(h2)
    return {"h2" : h2[0,0], "var(h2)" : 0}

