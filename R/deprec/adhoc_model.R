data("residential", package="RMoE")
df_ex <- residential
X_ex <- as.matrix(df_ex[-9])
Xs_ex <- scale_column_subset(X_ex, exclude=1)
Z_ex <- df_ex$V9
K <- 3
Lambda <- 15
gamma <- 5
model_ex <- RMoE::GaussRMoE(Xs_ex, as.vector(scale(Z_ex)), K, Lambda, gamma)
