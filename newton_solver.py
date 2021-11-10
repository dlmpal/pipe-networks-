def sign(P1,P2):
    if P1 > P2:
        return -1
    else:
        return  1
    
def sign2(i , N1 , N2 , h_vec , BC ):
    if i in BC:
        return 0 
    if i != N1 and i != N2:
        return 0 
    if i == N1:
        return sign(h_vec[N1] , h_vec[N2])
    if i == N2:
        return sign(h_vec[N2] , h_vec[N1])
    else:
        raise ValueError("Error in connectivity.")
        
def sign3(i , j , N1 , N2):
    if i != N1 and i != N2:
        return 0 
    else:
        if i==j:
            return -1 
        else: 
            if i == N1 and j == N2:
                return 1 
            elif i==N2 and j==N1:
                return 1 
            else:
                return 0         
                
    def newton_iter(NNODES,ELE,F,J,H,BC):
    for i in range(NNODES):
        for e in ELE:
            F[i] += np.sqrt(np.abs(H[e["N1"]-1]-H[e["N2"]-1])/e["k"]) * sign2(i , e["N1"]-1 , e["N2"]-1 , H ,BC)
            
    for i in range(NNODES):
        for j in range(NNODES):
            if i in BC:
                if i == j :
                    J[i,j] = 1
                else:
                    J[i,j] = 0 
                continue
            for e in ELE:
                J[i,j] += 1/(2*np.sqrt(np.abs(H[e["N1"]-1]-H[e["N2"]-1])*e["k"])) * sign3(i,j,e["N1"]-1,e["N2"]-1)
   
    return np.linalg.solve(J,F)

def newton(NNODES , ELE , F , H , BC , iters = 5):
    for i in range(iters):
        f = [F[i] for i in range(len(F))]
        J = np.zeros((NNODES,NNODES))
        H = H - newton_iter(NNODES,ELE,f,J,H,BC)
    return H 
    """
    Example element:
    ele = {"N1":N1,"N2":N2,"k":k}
    """
    
