# test build_GKM_connection in all three cases when the connection is uniquely defined:
#  1. valency >= 3 and 3-independent
#  2. valency == 2 and 2-independent
#  3. valency == 1

P1 = GKMproj_space(1)
c1 = build_GKM_connection(P1)

P2 = GKMproj_space(2)
c2 = build_GKM_connection(P2)

P3 = GKMproj_space(3)
c3 = build_GKM_connection(P3)