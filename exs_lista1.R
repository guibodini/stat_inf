# Implemente a função f(x)=x−−√ e desenhe seu gráfico no intervalo (0,3).
x <- seq(0, 3, l = 100)
fx <- function(x) {
  sqrt(x)
}
y <- fx(x = x)
plot(y ~ x, type = "l")

# Implemente a função f(x)=log(x) e desenhe seu gráfico no intervalo (−5,5).

x <- seq(-5, 5, l = 100)
fx <- function(x){
  log(x)
}
y <- fx(x = x)
plot(y ~ x, type = 'l')


#Calcule as derivadas parciais numéricas de primeira e segunda ordem da seguinte função de duas variáveis:
#u=a2x+3y.
#Use o método de Richardson com h=0.00001. 
#Avalie cada uma das derivadas no ponto x=1, y=2 e a=0.5. 
#A resposta são seis valores numéricos na seguinte ordem: ux;uy;uxx;uyy;uxy e uyx. 
#Use arredondamento com três casas decimais para a sua resposta.
# Funções de derivada central
deriv_funcs <- list(
  ux  = function(f, x, y, a, h) (f(x + h, y, a) - f(x - h, y, a)) / (2 * h),
  uy  = function(f, x, y, a, h) (f(x, y + h, a) - f(x, y - h, a)) / (2 * h),
  uxx = function(f, x, y, a, h) (f(x + h, y, a) - 2*f(x, y, a) + f(x - h, y, a)) / (h^2),
  uyy = function(f, x, y, a, h) (f(x, y + h, a) - 2*f(x, y, a) + f(x, y - h, a)) / (h^2),
  uxy = function(f, x, y, a, h) (
    f(x + h, y + h, a) - f(x + h, y - h, a) - f(x - h, y + h, a) + f(x - h, y - h, a)
  ) / (4 * h^2)
)

# Parâmetros
x <- 1; y <- 2; a <- 0.5; h <- 1e-5

# Extrapolação de Richardson
richardson <- function(D, D2) (4 * D2 - D) / 3

# Resultado em lista
resultados <- list()
for (nome in names(deriv_funcs)) {
  D  <- deriv_funcs[[nome]](fxya, x, y, a, h)
  D2 <- deriv_funcs[[nome]](fxya, x, y, a, h/2)
  resultados[[nome]] <- round(richardson(D, D2), 3)
}

# Mostrar resultados
print(resultados)


# Calcule∫422dx.

fa <- function(x) rep(2, length((x)))
resultado <- integrate(fa, lower = 2, upper = 4)
print(resultado$value)



# Instale o pacote se necessário: install.packages("statmod")
install.packages("fastGHQuad")
# Instale se necessário: install.packages("fastGHQuad")
library(fastGHQuad)
lambda <- 0.091

# Obtenha pontos e pesos de Gauss-Hermite de ordem 21
herm <- gaussHermiteData(21)
x <- herm$x
w <- herm$w

# Calcule a integral usando a mudança de variável
resultado <- lambda * exp(-lambda^2) * sum(w * exp(2 * lambda * x))

# Arredonde para três casas decimais
resultado <- round(resultado, 3)
cat("Resultado da integral:", resultado, "\n")


# Define a função
d <- function(y, mu){y * log(y / mu) - y + mu + 10}

# Valores de y
y <- seq(0.01, 20, length = 200) # começa em 0.01 para evitar log(0)

# Valores de mu para o gráfico
mu_vals <- c(1 ,3, 5, 8, 17)

# Plota
plot(y, d(y, mu_vals[1]), type = "l", col = 1, lwd = 2, ylim = c(min(d(y, mu_vals)), max(d(y, mu_vals))),
     ylab = "d(y, mu)", xlab = "y", main = "d(y, mu) para vários valores de mu")
for(i in 2:length(mu_vals)){
  lines(y, d(y, mu_vals[i]), col = i, lwd = 2)
}
legend("topright", legend = paste("mu =", mu_vals), col = 1:length(mu_vals), lty = 1, lwd = 2)

# Exemplo 1: função simples
y <- function(x)(-x+4)/2
y(0) #intercepto
a <- deriv((-x+4)/2, "x")
a
cat("O intercepto vertical é:", intercepto, "\n")

matriz <- matrix(c(1.0, 0.8, 0.8, 0.7, 0.8, 1.0, 0.7, 0.8,
                   0.8, 0.7, 1.0, 0.8, 0.7, 0.8, 0.8, 1.0), nrow=4, byrow=T)
print(matriz)
sum(diag(matriz))

# Instale se necessário:
# install.packages("fastGHQuad")

library(fastGHQuad)

# Função para integrar
gauss_hermite <- function(lambda, n=20) {
  integrando <- function(x) lambda * exp(-abs(x - lambda))
  pontos <- gaussHermiteData(n)
  f_x <- integrando(pontos$x)
  integral <- sum(pontos$w * f_x / exp(-pontos$x^2))
  return(integral)
}

# Valor de lambda
lambda <- 0.172

# Calculando a integral
resultado <- gauss_hermite(lambda)
resultado <- round(resultado, 3)
cat("Resultado da integral:", resultado, "\n")

A <- matrix(
  c(8, 7, 11, 15, 8,
    12, 5, 9, 11, 10,
    13, 11, 6, 8, 8,
    14, 9, 9, 10, 13,
    15, 9, 9, 10, 14),
  nrow=5, ncol=5, byrow=FALSE
)
A

B <- matrix(
  c(14, 13, 9, 9, 6,
    11, 11, 7, 12, 6,
    12, 10, 9, 14, 6,
    11, 3, 12, 9, 15,
    10, 9, 11, 9, 16),
  nrow=5, ncol=5, byrow=FALSE
)
B

C = A + B
C

total <- sum(C)
cat(total)


fb <- function(x) x^5
resultado <- integrate(fb, lower = 1, upper = 5)
print(resultado$value)


##
set.seed(123)
kappa <- 0.5
n <- 985

monte.carlo <- function(integrando, n) {
  pontos <- rnorm(n)
  norma <- dnorm(pontos)
  integral <- mean(integrando(pontos) / norma)
  return(integral)
}


# Definindo a função a ser integrada (com kappa fixo)
integrando <- function(x) {
  1/(kappa + 1/kappa) * exp(-x / kappa)
}


resultado <- monte.carlo(integrando, n)
round(resultado, 3)

