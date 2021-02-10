### FIGURES

#plotSimuScale(z1_1, w1_1, z2_1, w2_1, z3_1, w3_1)

##################

u1 <- colScale(min = min(z1_2,z1_3), max = max(z1_2,z1_3),
               epsilon = abs(max(z1_2,z1_3))/2, nb = 16)
v1 <- grayScale(max = max(w1_2,w1_3),  nb = 20)
u2 <- colScale(min = min(z2_2,z2_3), max = max(z2_2,z2_3), epsilon = 0.1, nb = 16)
v2 <- grayScale(max = max(w2_2,w2_3),  nb = 20)
u3 <- colScale(min = min(z3_2,z3_3), max = max(z3_2,z3_3), epsilon = abs(max(z3_2,z3_3)/4), nb = 16)
v3 <- grayScale(max = max(w3_2,w3_3),  nb = 20)

u1 <- colScale(min = min(z1_1,z1_2,z1_3), max = max(z1_1,z1_2,z1_3),
               epsilon = abs(max(z1_1,z1_2,z1_3))/2, nb = 16)
v1 <- grayScale(max = max(w1_1,w1_2,w1_3),  nb = 20)
u2 <- colScale(min = min(z2_1,z2_2,z2_3), max = max(z2_1,z2_2,z2_3), epsilon = 0.1, nb = 16)
v2 <- grayScale(max = max(w2_1,w2_2,w2_3),  nb = 20)
u3 <- colScale(min = min(z3_1,z3_2,z3_3), max = max(z3_1,z3_2,z3_3), epsilon = abs(max(z3_1,z3_2,z3_3)/4), nb = 16)
v3 <- grayScale(max = max(w3_1,w3_2,w3_3),  nb = 20)
plotSimu(z1_1, w1_1, z2_1, w2_1, z3_1, w3_1, u1, v1, u2, v2, u3, v3)
plotSimu(z1_2, w1_2, z2_2, w2_2, z3_2, w3_2, u1, v1, u2, v2, u3, v3)
plotSimu(z1_3, w1_3, z2_3, w2_3, z3_3, w3_3, u1, v1, u2, v2, u3, v3)


####################

u1 <- colScale(min = min(z1,z1_3), max = max(z1,z1_3),
               epsilon = abs(max(z1,z1_3))/2, nb = 16)
v1 <- grayScale(max = max(w1,w1_3),  nb = 20)
u2 <- colScale(min = min(z2,z2_3), max = max(z2,z2_3), epsilon = 0.1, nb = 16)
v2 <- grayScale(max = max(w2,w2_3),  nb = 20)
u3 <- colScale(min = min(z3,z3_3), max = max(z3,z3_3), epsilon = abs(max(z3,z3_3)/4), nb = 16)
v3 <- grayScale(max = max(w3,w3_3),  nb = 20)


plotSimu(z1_3, w1_3, z2_3, w2_3, z3_3, w3_3, u1, v1, u2, v2, u3, v3)
plotSimu(z1, w1, z2, w2, z3, w3, u1, v1, u2, v2, u3, v3)



