# Transform kernel matrix to centered matrix
kernel_centerer <- function(X){
  N_k <- nrow(X)
  one_N <- matrix(rep(1, N_k^2), nrow = N_k)/N_k
  return(X - one_N%*%X -X%*%one_N + one_N%*%X%*%one_N)
}

k_PCA <- function(X, gamma_k=1, n_components=0, kernel_='rbf',
                  teta=0, p=2){
  if(!require(purrr))
    install.packages("purrr")
  library(purrr)
  X <- as.matrix(X)
  K <- switch (kernel_,
               # rbf kernel
    'rbf' = (- gamma_k * (X %>% stats::dist(
    ) %>% as.matrix())^2) %>% exp(),
    # sigmoid kernel
    'sigmoid' = (teta * X%*%t(X) + gamma_k) %>% tanh(),
    # polynomial kernel
    'poly' = (X%*%t(X) + gamma_k)^p,
    # cosine kernel
    'cosine' = (X%*%t(X)) %>% cov()
    ) %>% kernel_centerer() %>% eigen(symmetric = TRUE)
  # récupération des indices des valeurs propres non nulles
  non_zeros_ind <- which(K$values != 0)
  # récupération des vecteurs propres dont les valeurs propre
  # sont non nulles
  vects <- K$vectors[, non_zeros_ind]
  n_ <- ncol(X)
  # Labelisation des composantes principales
  colnames(vects) <- paste('CP', 1:nrow(X), sep = '')
  if(n_components == 0)
    return(vects[, 1:n_])
  return(vects[, 1:min(n_, n_components)])
}

make_moons <- read.csv("~/Desktop/M2-SAAD/stage_proj/make_moons.csv", sep = ";")
X <- make_moons[, 1:2]
y <- make_moons[, 3]
plot(X[y==0, ],
     pch = puces[1], col = colors_[1], lwd = 1, ylab = "", xlab = "",
     xlim = c(min(X[, 1]) -.1,
              max(X[, 1])+.1),
     ylim = c(min(X[, 2]) -.1,
              max(X[, 2])+.2)
     )
points(X[y==1, ], pch = puces[2], col = colors_[2], lwd=1)
grid()

moons_eig_vecs <- k_PCA(X, gamma_k = 15, kernel_ = 'rbf',
                        teta = 1, p=5)[, 1:2]

plot(moons_eig_vecs[y==0, ],
     pch = puces[1], col = colors_[1], lwd = 1, ylab = "", xlab = "",
     xlim = c(min(moons_eig_vecs[, 1]) ,
              max(moons_eig_vecs[, 1])),
     ylim = c(min(moons_eig_vecs[, 2]) ,
              max(moons_eig_vecs[, 2])))
points(moons_eig_vecs[y==1, ], pch = puces[2], col = colors_[2], lwd=1)
grid()
