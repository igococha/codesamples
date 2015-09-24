__author__ = 'Igor Siveroni'

import sys
import phylotree
reload(phylotree)
from phylotree import PhyloTree

import numpy as np

#treeFile = 'data/Spread_N100.phy'
treeFile = 'data/test.newick'
#treeFile = 'data/Birds-01000.newick'

#tree = Phylo.read(treeFile,'newick')
#print type(tree)
#print tree.get_terminals()
#print type(tree.root)


## Dendro example

tree = PhyloTree(treeFile, schema='newick')

#print tree.as_string("newick")
# tree.printTree()

V = tree.getV()

print V

mylog = tree.getLogDet()

nplog = np.linalg.slogdet(V)

print "mylog", mylog
print "nplog", nplog

npinv = np.linalg.inv(V)

# print npinv
# print "my INV"
# myinv = tree.getVinv()
# print myinv

sys.exit()

#print type(tree)
# length: sum of all edge lengths
print "Length: " + str(tree.length())
#If the tree is not rooted, the program sets node.parent_node to None and adds
# a dummy edge with length and tail_node set to None
print "is rooted: " + str(tree.is_rooted)

# Lists the object's attributes
print(tree.__dict__)
# Lists an object's functions
# for ele in dir(tree):
#    if (callable(getattr(tree,ele))):
#        print(ele)

# a---b---c
#     |---d
# edge (b,c) b is tail, c is head   tail----head  internal---tip
# Only heads are TIPS

# EDGES
internal_counter = 0
#for edge in tree.level_order_edge_iter():
for edge in tree.postorder_edge_iter():
    #print "edge:" + str(edge) + "--"
    tail = edge.tail_node
    head = edge.head_node
    if (not head.is_leaf()):
        internal_counter += 1
    else:
        print(head.taxon)





#print node
#print node.taxon
#print node.parent_node
#print node.edge

# It's possible to add fields to objects at runtime
#node.value = 6
#print node.value
#print dtree.seed_node.value

print "last edge"
print "head",head
print "tail", tail
print "length", edge.length

print "seed node:",tree.seed_node


leafs = tree.leaf_nodes()
print "num tips", len(leafs)
print "num internal", len(tree.internal_nodes())
print "computed internal", internal_counter
print "edge set",len(tree.get_edge_set())

print tree.seed_node.parent_node

exit()



print "--- iteration ---"
# Node traversal
msg = "Num children:"
for node in tree.preorder_node_iter():
    msg += str(len(node.child_nodes())) + "  "
#msg = "--"+str(node)
#first_node = node
#if (node.parent_node is None):
#	msg = msg+"--ROOT--"
#if (node.taxon is not None):
#	msg = msg + str(node.taxon)
#if (node.edge is not None):
#	#msg = msg + "NO EDGE!!!!"
#	msg = msg + "--edge:"+str(node.edge)
#if node.is_leaf():
#	msg = msg + "  LEAF"
#msg = msg+" level : " + str(node.level())
#print msg



#edge_counter = 0
#for edge in tree.level_order_edge_iter():
#	print "edge:" + str(edge) + "--"
#	tail = edge.tail_node
#	head = edge.head_node
#	print head
#	print tail
#	print edge.length
#	#if (first_node == head):
#	#	print "OK!!"
#	#msg += str(type(tail))
#	#print edge.length
#	#if (tail.is_leaf()):
#	#	msg += " tail is leaf"
#	#else:
#	#	msg += " tail is NOT leaf"
#	#msg  edge.head_node
#	#print msg
#	break

print "---inspection---"
#node = tree.seed_node
#children = node.child_nodes()
#node = children[0]
#children = node.child_nodes()
#for child in children:
#	print child
#	print len(child.child_nodes())
#	if child.is_leaf():
#		print "Child"
#		print child.taxon
#	else:
#		print "Internal"
#	print child.edge.length

#tree.calc_node_ages(check_prec=False)

# Compute path lengths
sum = 0.0
for node in tree.preorder_node_iter():
    parent = node.parent_node
    if (parent is not None):
        node.path_length = node.edge.length + parent.path_length
    else:
        node.path_length = 0.0

#	if node.edge.length is not None:
#		sum += float(node.edge.length)
#
#print sum
#print tree.length()





