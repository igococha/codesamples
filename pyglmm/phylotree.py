__author__ = 'Igor Siveroni'

# Dendropy library for manipulation of phylogenetic trees
import dendropy
import numpy as np
from scipy.sparse import *

# Phylogeny tree representation useful for matrix generation
class PhyloTree:
     # Static variables
    num_trees = 0
    # Constructor
    def __init__(self,filename,schema='newick'):
        # Read phylogenetic tree using schema name (argument)
        self.dtree = dendropy.Tree(stream=open(filename), schema=schema)
        # preliminaries
        ntips = len(self.dtree.leaf_nodes())
        ninternal = len(self.dtree.internal_nodes())
        ntotal = ntips + ninternal

        # Node -> children information
        self.max_children = max([len(node.child_nodes()) for node in self.dtree.internal_nodes()])
        self.children = np.zeros(ninternal*self.max_children,dtype=np.int16)
        #num_children[nid]==0 iff nid is tip
        self.num_children = np.zeros(ntotal,dtype=np.int8)

        # Initialization
        self.num_tips = 0
        self.num_internal = 0
        self.num_nodes = 0
        self.first_tip = ninternal
        self.height = 0

        self.distances = np.zeros(ntotal)
        self.tipNames = np.empty(ntips,dtype='object')
        self.ultrametric = None
        self.ltip = np.zeros(ntotal,dtype=np.int16)
        self.rtip = np.zeros(ntotal,dtype=np.int16)
        # first_child and rall point to entries in the preorder_nodes array
        self.position_preorder = np.zeros(ntotal,dtype=np.int16)
        self.lastPosition_preorder = np.zeros(ntotal,dtype=np.int16)
        # edge_child[id] = id -- No need for edge_child
        self.edge_parent = np.zeros(ntotal,dtype=np.int16)
        self.edge_length = np.zeros(ntotal)
        self.preorder_nodes = np.zeros(ntotal,dtype=np.int16)
        #self.__init_structures()
        #print self.__computeHeight(self.dtree.seed_node,0)
        self.distances[0] = 0.0
        #print "Root edge:",self.dtree.seed_node.edge_length
        self.hasRootEdge = (self.dtree.seed_node.edge_length is not None)
        (node_id,l,r) = self.__traverse(self.dtree.seed_node,None,0,0)
        if (self.num_nodes != self.num_tips+self.num_internal):
            raise Exception("Computed number of nodes incompatible with original tree")
        # We are done with dendropy tree. Remove it from scope (and eventually from memory)
        del self.dtree

        # Class Counter
        PhyloTree.num_trees += 1

        # Matrix related
        self.__V = None
        self.__logDet = None
        self.__Vinv = None
        self.__S = None
        self.__Sinv = None

    # Preliminary version using iterators
    def __init_structures(self):
         print("Traversing tree")
         # traverse dendropy tree edges
         # indexing internal nodes [0,num_internal-1] preorder
         # indexing tips [num_internal,num_total]
         t = self.dtree
         d = self.distances
         d[0] = 0.0
         internal_counter = 0
         tip_counter = self.num_internal
         for edge in t.preorder_edge_iter():
             child = edge.head_node
             parent = edge.tail_node
             if (child.is_leaf()):
                 child.id = tip_counter
                 tip_counter += 1
             else:
                 child.id = internal_counter
                 internal_counter  += 1
             if (parent is not None): # Not the root
                 d[child.id] = d[parent.id] + edge.length
         #for edge in t.postorder_edge_iter():

    # Example of a simple traversal to compute height
    def __computeHeight(self,node,height):
       if (node.is_leaf()):
            return height
       else:
           return max([self.__computeHeight(child,height+1) for child in node.child_nodes()])

    # Structures initialization using traversal
    def __traverse(self,node,parent,h,distance):
        edge = node.edge
        length = node.edge_length
        #print self.num_nodes,length

        if (parent is None):
            parent_id = -1
            length = (0 if (node.edge_length is None) else node.edge_length)
        else:
            parent_id = parent.id
            length = node.edge_length
            if (length is None):
                raise Exception('Edge length missing. Check tree.')

        new_distance = distance + length
        lsl = []
        lsr = []
        if (node.is_leaf()):
            #msg="TIP:"
            self.tipNames[self.num_tips] = node.get_node_str()
            if (h > self.height):
                self.height=h
            node_id = self.first_tip + self.num_tips
            self.num_tips += 1
        else:
            # pre-traversal
            #msg="NODE:"
            node_id = self.num_internal
            self.num_internal += 1
            # (children_ino) = [self.__traverse(child,node,h+1,new_distance) for child in node.child_nodes()]

        self.edge_length[node_id] = length
        #print "num_nodes",self.num_nodes,msg,node_id

        # annotate tree node
        node.id = node_id
        self.distances[node_id] = new_distance
        # edge has same id as its child node
        self.edge_parent[node_id] = parent_id
        current_idx  = self.num_nodes
        self.preorder_nodes[current_idx] = node_id
        self.position_preorder[node_id] = current_idx
        self.num_nodes += 1

        # traversal
        self.num_children[node_id] = len(node.child_nodes())
        child_idx = self.max_children*node_id
        for child in node.child_nodes():
            (child_id,childl,childr) = self.__traverse(child,node,h+1,new_distance)
            self.children[child_idx] = child_id
            child_idx += 1
            lsl.append(childl)
            lsr.append(childr)

        self.lastPosition_preorder[node_id] = self.num_nodes - 1

        # post-traversal
        if (node.is_leaf()):
            l = r = node_id
        else:
            l = min(lsl)
            r = max(lsr)
        # print msg,"(node,parent,lsl,l,lsr,r)",(node_id,parent_id,lsl,l,lsr,r)
        #node.id = node_id
        self.ltip[node_id] = l
        self.rtip[node_id] = r
        return (node_id,l,r)

    def __testTraverse(self,node):
        self.counter += 1
        if (self.counter > 20):
            return
        n = self.num_children[node]
        if (n==0):
            print "TIP ",node
        else:
            print "INT ", node
            idx = self.max_children*node
            for i in range(idx,idx+n):
                self.__testTraverse(self.children[i])


    # No longer valid - dtree was made null
    #def printTree(self):
    #    self.dtree.print_newick()

    def isUltrametric(self):
        if (self.ultrametric is None):
            self.ultrametric = True
            d = self.distances[self.first_tip]
            for tip in range(self.first_tip+1,self.num_nodes):
                if (abs(d - self.distances[tip]) > 0.0001 ):
                    self.ultrametric = False
                    break
        return self.ultrametric

    # Normalizes V. If V is ultrametric, it divides the matrix by its diagnonal
    # If it is not ultrametric, it divides the matrix by its max diagonal
    # The latter is probabaly nonsense
    def normalize(self):
        V = self.getV()
        if (self.isUltrametric()):
            d = V[0,0]
        else:
            # M.max() whole matrix, M.max(0) max ech column, M.max(1) -> rows
            d = max(np.diag(V))
        V.flags.writeable = True
        V /= d
        V.flags.writeable = False

     # #################### V matrix #####################
    def getV(self):
        if (self.__V is None):
            self.__computeV()
        return self.__V

    def __computeV(self):
        # read structures
        ntips = self.num_tips
        distances = self.distances
        anum_children = self.num_children
        max_children = self.max_children
        ninternal = self.num_internal
        children = self.children
        ltip = self.ltip
        rtip = self.rtip
        # create matrix
        V = np.matrix(np.zeros((ntips,ntips))) # np.asarray(V) if we want to use the ndarray rep.
        # Tips / Diagonals
        node_id = self.first_tip
        while (node_id < self.num_nodes):
            index = node_id - ninternal
            V[index,index] = distances[node_id]
            node_id += 1
        # Internal nodes
        for node_id in range(0,self.num_internal):  # changed from num_tips to num_internal
            d = distances[node_id]
            num_children = anum_children[node_id]
            start = node_id*max_children
            children_id = list(children[range(start,start+num_children)])
            lr = [(ltip[child],rtip[child]) for child in children_id]
            #print "node",node_id,"lr",lr
            for i in range(0,num_children-1):
                (l1,r1) = lr[i]
                #print "**",l1,r1
                for (l2,r2) in lr[i+1:]:
                    #print "***",l2,r2
                    for j in range(l1-ninternal,r1+1-ninternal):
                        for k in range(l2-ninternal,r2+1-ninternal):
                            #print "*****",j,k
                            V[j,k] = d
                            V[k,j] = d
        # Changing V will require changes to the other structures eg distances.
        # Make it read-only but allow changes through interface.
        V.flags.writeable = False
        self.__V = V

    # ################### Log-Determinant of V ######################
    def getLogDet(self):
        #if (self.__logDet is None):
        #    self.__computeLogDet()
        self.__computeLogDet()
        return self.__logDet

    # Determinant computation based on three-point structure
    # Recall: descendant/child(edge_id) = edge_id = node_id
    # Timings
    # tree 1000 np 28-38ms, 3point 27-39ms
    #      2000 np 101-115, 3point 40-49
    #      3000 np 247-255, 3point 61-86
    def __computeLogDet(self):
        # Intermediate results
        vec11 = np.zeros(self.num_nodes) # 1'V^{-1}1
        logd = np.zeros(self.num_nodes)  # logdetV
        for i in range(self.num_nodes-1,0,-1): #perform post-order traversal
            edge_id = node_id = self.preorder_nodes[i]
            el = self.edge_length[edge_id]
            parent_id = self.edge_parent[edge_id]
            #print "visited:", edge_id, self.num_children[node_id],el
            if (node_id >= self.first_tip): # TIP
                #print("tip el=",el)
                if (el > 0):
                    logd[node_id] = np.log(el)
                    vec11[node_id] = 1/el
                else:
                    raise Exception("Found Zero edge lenght at tip")
            else:
                logd[node_id] = logd[node_id]+np.log(1+el*vec11[node_id])
                vec11[node_id] = vec11[node_id]*(1 - 1/(1+1/el/vec11[node_id]))
                # vec11[node_id] = vec11[node_id]*(1 - el/(el+vec11[node_id]))
            # Update/Collect at parent - simple sums in this case
            logd[parent_id] = logd[parent_id] + logd[node_id]
            vec11[parent_id] = vec11[parent_id] + vec11[node_id]
        # Last Node/edge
        edge_id = node_id = self.preorder_nodes[0]
        if (self.hasRootEdge):
            edge_id = node_id = self.preorder_nodes[0]
            el = self.edge_length[edge_id]
            # assert(self.edge_parent[edge_id] = -1
            # Note that el could be zero (if it hasn't been fixed before)
            #print "Root edge",self.edge_parent[node_id],node_id,self.edge_length[node_id]
            logd[node_id] = logd[node_id]+np.log(1+el*vec11[node_id])
        self.__logDet = logd[node_id]


    # #################### V inverse ################################
    def getVinv(self):
        #if (self.__Vinv is None):
        #    self.__computeVinv()
        self.__computeVinv()
        return self.__Vinv

    # Very slow
    # 1000 tips numpy 60ms, 3point 6.34ms
    # 2000 tips numpy 338ms (x5.5), 3point 27.7s (x4.5)
    # 3000 tips numpy 820ms (x13.6), 3point 66s (x10.4)
    def __computeVinv(self):
        ntips = self.num_tips
        ntotal = self.num_nodes
        anum_children = self.num_children
        max_children = self.max_children
        ninternal = self.num_internal
        children = self.children
        ltip = self.ltip
        rtip = self.rtip
        # create matrix
        Vinv = np.matrix(np.zeros((ntips,ntips)))
        vec11 = np.zeros(ntotal) # 1'V^{-1}1
        sumRows = np.zeros(ntips)
        # Main LOOP
        for i in range(self.num_nodes-1,0,-1): #perform post-order traversal
            edge_id = node_id = self.preorder_nodes[i]
            el = self.edge_length[edge_id]
            parent_id = self.edge_parent[edge_id]
            if (node_id >= self.first_tip): # TIP
                tip = node_id - ninternal
                #print "tip",tip,"el",el
                if (el > 0):
                    elinv = 1/el
                    vec11[node_id] = elinv
                    Vinv[tip,tip] = elinv
                    sumRows[tip] = elinv
                else:
                    raise Exception("Found zero edge length at tip")
            else:
                l = ltip[node_id]-ninternal # l in [0,ntips-1]
                r = rtip[node_id]-ninternal # r in [0,ntips-1]
                #print "Internal",node_id,"l",l,"r",r,"el",el
                tPa = el/(1+el*vec11[node_id])
                #print(vec11[node_id])
                #print(tPa)
                for j in range(l,r+1):
                    sumj = sumRows[j]
                    for k in range(l,j+1):
                        v = Vinv[j,k] - sumj*sumRows[k]*tPa
                        Vinv[j,k] = v
                        Vinv[k,j] = v
                # update sumRows=sumCols
                sumRows[l:r+1] = 0.0
                for j in range(l,r+1):
                    sumRows[j] = np.sum(Vinv[j,l:r+1])
                vec11[node_id] = vec11[node_id]*(1 - 1/(1+1/el/vec11[node_id]))
                # IMPORTANT: the formula below is equivalent but quickly blows up
                #vec11[node_id] = vec11[node_id]*(1 - el/(el+vec11[node_id]))
            # ------------- Join --------------
            # Update/Collect at parent - simple sums for vec11
            vec11[parent_id] = vec11[parent_id] + vec11[node_id]
        # Check root
        node_id = edge_id = self.preorder_nodes[0]
        if (self.hasRootEdge):
            l = ltip[node_id]-ninternal # l in [0,ntips-1]
            r = rtip[node_id]-ninternal # r in [0,ntips-1]
            el = self.edge_length[edge_id]
            #print "Internal",node_id,"l",l,"r",r,"el",el
            tPa = el/(1+el*vec11[node_id])
            #print(vec11[node_id])
            #print(tPa)
            for j in range(l,r+1):
                sumj = sumRows[j]
                for k in range(l,j+1):
                    v = Vinv[j,k] - sumj*sumRows[k]*tPa
                    Vinv[j,k] = v
                    Vinv[k,j] = v
        self.__Vinv = Vinv

    # Main result: X'V{-1}y
    # Other: X'V^{-1}X, X'V^{-1}1,y'V{-1}y,y'V^{-1}1,1'V^{-1}1
    # Input: y is a uni-dimensional ndarray , X is a numpy matrix
    # y should be (if bi-dimensional) the transpose of the y-data
    def product(self,y,X):
        ntips = self.num_tips # formerly n
        ntotal = self.num_nodes # formerly N, number edges = number nodes
        ninternal = self.num_internal
        root = self.preorder_nodes[0]
        (nx,p) = np.shape(X)  # p formerly d
        yshape = np.shape(y)
        if (np.size(yshape)==1):
            m = 1
            ny = yshape[0]
        else:
            m = yshape[0]
            ny = yshape[1] # other dimensions ignored
        if (nx != ny):
            raise Exception("Size of X and y don't agree")
        if (nx != ntips):
            raise Exception("Sizes of X and y do not agree with number of tips in three")
        zero = np.zeros(ntotal)
        vec11 = np.zeros(ntotal) # 1'V^{-1}1
        y1 = np.zeros(ntotal) # y'V^{-1}1
        yy = np.zeros(ntotal) # y'V^{-1}y
        logd = np.zeros(ntotal) # logdetV
        X1 = np.matrix(np.zeros((ntotal,p))) # X'V^{-1}1
        XX = np.zeros((ntotal,p,p)) # X'V^{-1}X
        Xy = np.matrix(np.zeros((ntotal,p))) # X'V^{-1}y
        # The LOOP
        for i in range(self.num_nodes-1,0,-1):
            #perform post-order traversal
            edge_id = node_id = self.preorder_nodes[i]
            el = self.edge_length[edge_id]
            parent_id = self.edge_parent[edge_id]
            if (node_id >= self.first_tip): # TIP
                #print("tip el=",el)
                if (el > 0):
                    tip_id =  node_id - ninternal
                    #print "node",node_id,"tip",tip_id
                    logd[node_id] = np.log(el)
                    Xy[node_id,:] = X[tip_id,:] * y[tip_id]/el
                    yy[node_id] = y[tip_id]**2/el
                    y1[node_id] = y[tip_id]/el
                    XX[node_id,:,:] = (X[tip_id,:].T * X[tip_id,:])/el
                    X1[node_id,:] = X[tip_id,:]/el
                    vec11[node_id] = 1/el
                else:
                    raise Exception("Found Zero edge lenght at tip")
            else:
                logd[node_id] = logd[node_id]+np.log(1+el*vec11[node_id])
                k = (1 - 1/(1+1/el/vec11[node_id]))
                k2 = el/(1+el*vec11[node_id])
                Xy[node_id,:] -=  (k2*y1[node_id])*X1[node_id,:]
                yy[node_id] -=  k2 * y1[node_id]**2
                y1[node_id] *= k
                XX[node_id,:,:] -=  k2 * X1[node_id,:].T * X1[node_id,:]
                X1[node_id,:] *= k
                # vec11[node_id] = vec11[node_id]*(1 - el/(el+vec11[node_id]))
                vec11[node_id] = vec11[node_id]*k
            # Update/Collect at parent - simple sums in this case
            logd[parent_id] += logd[node_id]
            Xy[parent_id,:] +=  Xy[node_id,:]
            yy[parent_id] += yy[node_id]
            y1[parent_id] += y1[node_id]
            XX[parent_id,:,:] += XX[node_id,:,:]
            X1[parent_id,:] += X1[node_id,:]
            vec11[parent_id] += vec11[node_id]

        #print "ROOT",root
        return Xy[root,:]

    # ################## Henderson's inverse of full distance matrix ######################

    def getS(self):
        if (self.__S is None):
            self.computeS()
        return self.__S

    # Calculation of a complete V matrix
    def computeS(self):
        max_children = self.max_children
        traversal =  self.preorder_nodes
        root = 0
        # create matrix
        S = np.matrix(np.zeros((self.num_nodes-1,self.num_nodes-1)))
        # Tips / Diagonals -- skipping root=0
        np.fill_diagonal(S,self.distances[1:])
        # Now the rest - internal nodes - no root for now
        for node_id in range(1,self.num_internal):
            d = self.distances[node_id]
            num_children = self.num_children[node_id]
            start = node_id*max_children
            children_id = list(self.children[range(start,start+num_children)])
            lr = [(self.position_preorder[child],self.lastPosition_preorder[child]+1) for child in children_id]
            # Match subtrees' children against each other
            for i in range(0,num_children-1):
                (l1,r1) = lr[i]
                for (l2,r2) in lr[i+1:]:
                    for j in traversal[l1:r1]:
                        for k in traversal[l2:r2]:
                            S[j-1,k-1] = d
                            S[k-1,j-1] = d
            # Match children with parent - all share the same path
            l = self.position_preorder[node_id]
            r = self.lastPosition_preorder[node_id]+1
            for child_id in traversal[l+1:r]:
                S[node_id-1,child_id-1] = d
                S[child_id-1,node_id-1] = d
        self.__S = S

    def getSinv(self):
        if (self.__Sinv is None):
            self.computeSinv()
        return self.__Sinv

    def computeSinv(self):
        Sinv = np.matrix(np.zeros((self.num_nodes-1,self.num_nodes-1)))
        root = self.preorder_nodes[0] # should be zero
        for i in range(1,self.num_nodes,1):
            #perform post-order traversal
            edge_id = self.preorder_nodes[i]
            el = self.edge_length[edge_id]
            node_idx = edge_id-1
            parent_id = self.edge_parent[edge_id]
            parent_idx = parent_id - 1
            # Diagonal
            v = 1.0/el
            Sinv[node_idx,node_idx] = v
            if (parent_id != root):
                Sinv[parent_idx,parent_idx] += v
                Sinv[node_idx,parent_idx] = -v
                Sinv[parent_idx,node_idx] = -v
        self.__Sinv = Sinv

    # Returns S inverse stored as a sparse matrix
    # Avoids the allocation of n^2
    def computeSinvSparse(self, Acorner = "bottom"):
        m = self.num_nodes - 1
        ninternal = self.num_internal
        ntips = self.num_tips
        root = self.preorder_nodes[0] # should be zero
        data_size = 3*(m) - 2*self.num_children[root]
        data = np.zeros(data_size)
        data_diag = np.zeros(m)
        data_ij = np.zeros(2*data_size,dtype=np.int8).reshape(2,data_size)
        data_idx = 0
        for i in range(1,self.num_nodes,1):
            #perform post-order traversal
            edge_id = self.preorder_nodes[i]
            el = self.edge_length[edge_id]
            parent_id = self.edge_parent[edge_id]
            # index calculation
            if (Acorner == "bottom"):
                node_idx = edge_id-1
                parent_idx = parent_id - 1
            else:
                node_idx = edge_id + (ntips-1 if (edge_id < ninternal) else -ninternal)
                parent_idx = parent_id + (ntips-1 if (parent_id < ninternal) else -ninternal)
            # Diagonal
            v = 1.0/el
            # Sinv[node_idx,node_idx] = v
            data_diag[node_idx] = v
            if (parent_id != root):
                # Sinv[parent_idx,parent_idx] += v
                data_diag[parent_idx] += v
                # Sinv[node_idx,parent_idx] = -v
                data[data_idx] = -v
                data_ij[0,data_idx] = node_idx
                data_ij[1,data_idx] = parent_idx
                data_idx += 1
                # Sinv[parent_idx,node_idx] = -v
                data[data_idx] = -v
                data_ij[0,data_idx] = parent_idx
                data_ij[1,data_idx] = node_idx
                data_idx += 1
        print "num nodes: "+str(self.num_nodes)
        print "allocated data: "+str(data_size)
        print "data used: "+str(data_idx)
        # Diagonal data
        data[data_idx:] = data_diag
        data_ij[0,data_idx:] = range(0,m)
        data_ij[1,data_idx:] = range(0,m)
        return csr_matrix((data,data_ij),shape=(m,m))