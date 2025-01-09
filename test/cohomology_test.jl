P3 = GKMproj_space(3)
R = equivariant_cohomology_ring(P3)
Q = R.coeffRing
(t1,t2,t3,t4) = gens(Q)

@test isGKMclass(R.cohomRing([t1, t2, t3, t4]), R)
@test !isGKMclass(R.cohomRing([Q(1), Q(0), Q(0), Q(0)]), R)

println(R)

p1 = pointClass(1, R)
p2 = pointClass(2, R)

@test p1 == eulerClass(1,R)*R.cohomRing([Q(1), Q(0), Q(0), Q(0)])
@test multiply(p1, p2, R) == zero(R)

println("Integral of p1 is $(integrateGKMClass(p1, R)).")
println("Integral of p1^2 is $(integrateGKMClass(multiply(p1,p1,R), R)).")
println("Integral of (1,0,0,0) is $(integrateClass(R.cohomRing([Q(1), Q(0), Q(0), Q(0)]), R)).")