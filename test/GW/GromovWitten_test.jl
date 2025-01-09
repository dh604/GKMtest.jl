# example for P^3
P3 = GKMproj_space(3)
R = equivariant_cohomology_ring(P3)
con = build_GKM_connection(P3)

# build decorated tree, which is 1-2-3, mapping 1 -> 1, 2-> 4, 3->3, and marked points (1, 2) at 1, and (3) at 2.
tree = Graph{Undirected}(3)
add_edge!(tree, 1, 2)
add_edge!(tree, 2, 3)

vDict = Dict(1=>1, 2=>4, 3=>3)

edgeMult = Dict(Edge(1, 2) => 1, Edge(2, 3) => 1)

marks = [1, 1, 2]

dt = decoratedTree(P3, tree, vDict, edgeMult, marks)

C = R.coeffRing
(t1, t2, t3, t4) = gens(C)
H = R.cohomRing
c1 = pointClass(1, R) #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(1, R) #gens(H)[1] # pointClass(1, R) # H([C(0), C(0), C(0), C(1)]) # gens(H)[2]
c3 = pointClass(4, R) #gens(H)[4] # pointClass(4, R)
classes = [c1, c2, c3]

GWTreeContribution(dt, R, con, classes)