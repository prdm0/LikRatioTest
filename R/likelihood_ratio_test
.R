#' @title Teste da razão de verossimilhança generalizada
#' @author Pedro Rafael D. Marinho
#' @description Calcula a estatística da razão de verossimilhança generalizada.
#' @details A função recebe como argumento uma função densidade de probabilidade ou uma função de probabilidade
#' que é passada como argumento à \code{f}. O objeto passado como argumento de \code{f} deverá ser implementado conforme os
#' exemplos. Note que se a função for implementada segundo os exemplos, não haverá a necessidade de implementar a função sob a
#' hipótese nula, visto que são as mesmas funções. Para especificar os parâmetros que serão fixados, i.e, para especificar
#' a distribuição sob a hipótese nula, utiliza-se o argumento \code{par0} que receberá uma lista formada por dois vetores.
#' O primeiro vetor da lista de verá ser um vetor de strings com os nomes das variáveis que deseja-se fixar e o segundo vetor
#' deverá conter os valores que serão atribuidos à cada uma das variáveis passada ao primeiro vetor.
#' @export
#' @param f Função densidade de probabilidade considerada no teste da razão de verossimilhança. Essa função deverá ser
#' implementada conforme os exemplos.
#' @param data Um conjunto de dados que será considerado para a realização do teste.
#' @param par0 Uma lista contendo como primeiro elemento um vetor com o nomes das variáveis que serão fixadas como hipótese nula
#' (variáveis aos quais desejamos testar) e um vetor com os valores fixados para cada uma das respectivas variáveis.
#' @param ... Lista de argumentos adicionais que serão passados à função \code{optim()} otimizada para otimização. Por exemplo,
#' será possível escolher o método de otimização a ser utilizado.
#' @return Valor da estatística do teste da razão de verossimilhança generalizada.
#' @examples
#'  pdf_w <- function(par, x, var = NULL){
#'      alpha <- par[1]
#'      beta <- par[2]
#'      if (is.list(var)) eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
#'      dweibull(x, shape = alpha, scale = beta)
#'   }
#'
#'   rw <- function(n = 1L, alpha, beta){
#'    rweibull(n = n, shape = alpha, scale = beta)
#'   }
#'
#'   data <- rw(n = 100L, alpha = 1, beta = 1)
#'
#'   likelihood_ratio_test(f = pdf_w, data = data, kicks = c(1, 1), par0 = list("beta", 1))

likelihood_ratio_test <- function(f, data, kicks, par0 = NULL, ...){

  if (is.null(par0)) stop("Informar uma lista informando o parâmetro e o valor sob a hipótese nula.")

  # Log-Likelihood under the null hypothesis. -----------------------------------
  log_lik_h0 <- function(par, x){
    -sum(log(f(par, x, var = par0)))
  }

  # Log-Likelihood under the null hypothesis. -----------------------------------
  log_lik <- function(par, x){
    -sum(log(f(par, x)))
  }

  myoptim <- function(...) tryCatch(expr = optim(...), error = function(e) NA)

  par_h0 <- myoptim(par = kicks, fn = log_lik_h0, x = data, ...)

  if(!is.list(par_h0) || par_h0$convergence != 0L) return(NA)

  par_h <- myoptim(par = kicks, fn = log_lik, x = data, ...)

  if(!is.list(par_h) || par_h$convergence != 0L) return(NA)

  lambda <- 2 * (log_lik_h0(par = par_h0$par, x = data) - log_lik(par_h$par, x = data))

  lambda[lambda < 0] <- 0

  # Estatística de razão de verossimilhança:
  lambda
}

#' @title Simulações de Monte-Carlo para o teste da razão de verossimilhança generalizado.
#' @author Pedro Rafael D. Marinho
#' @description Realiza uma única iteração de um procedimento de Monte-Carlo para o teste da razão de verossimilhança generalizado. Dado um nível
#' de significância, será retornado 1 (um) se a estatística de teste está acima do quantil da distribuição qui-quadrado e 0 (zero),
#' caso contrário.
#' @details bla bla
#' @param i Para uma única realização de um procedimento de Monte-Carlo, esse parâmetro não terá importância. Esse perâmetro será
#' utilizado em situações em que utiliza-se de um funcional para reproduzir a função \code{mc}.
#' @param n Tamanho da amostra a ser considerada.
#' @param sig Nível de significância adotado.
#' @param f Função densidade de probabilidade considerada no teste. Essa função deverá ser implementada conforme o exemplo abaixo.
#' @param q Função responsável pela geração de observações de uma variável aleatório com função densidade passada para \code{f}.
#' @param kicks Vetor com os chutes iniciais utilizados para a otimização.
#' @param par0 Lista com dois elementos, sendo o primeiro um vetor com os nomes das variáveis que receberão valores fixos sob a
#' hipótese nula e o segundo elemendo é um outro vetor com os valores impostos às variáveis.
#' @param ... Lista de argumetos que serão passados para a função passada à \code{q}.
#' @export
#' @return Retornará 0 (zero) se a estatística calculado não estiver acima do quantil da distribuição qui-quadrado e 1 (um),
#' caso contrário.
#' @examples
#'pdf_ew <- function(par, x, var = NULL){
#'  alpha <- par[1]
#'  sigma <- par[2]
#'  theta <- par[3]
#'
#'  if (is.list(var)) eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
#'
#'  alpha * theta / sigma * (1 - exp(-(x / sigma) ^ alpha)) ^ (theta - 1) * exp(-(x / sigma) ^ alpha) * (x / sigma) ^ (alpha - 1)
#'}
#'
#'rew <- function(n, alpha, sigma, theta){
#'  u <- runif(n, 0, 1)
#'  sigma * (-log(1 - u ^ (1 / theta))) ^ (1 / alpha)
#'}
#'
#'set.seed(1L, kind = "L'Ecuyer-CMRG")
#'
#'tictoc::tic()
#'rejeicao <- unlist(pbmcapply::pbmclapply(X = 1L:1e3L, FUN = mc,
#'                                         mc.cores = parallel::detectCores(), f = pdf_ew, q = rew,
#'                                         sig = 0.05, n = 100L, kicks = c(1, 1, 1),
#'                                         par0 = list("beta", 1.5),
#'                                         alpha = 1.7, beta = 1.5))
#'
#'# Proporção de rejeição ---------------------------------------------------
#'sum(rejeicao)/length(rejeicao)
#'tictoc::toc()
# Simulação de Monte-Carlo ------------------------------------------------
mc <- function(i, n, sig = 0.05, f, q, kicks, par0, ...){
  amostra <- q(n, ...)
  result <- likelihood_ratio_test(f = f, data = amostra, kicks = kicks,
                                  par0 = par0)

  # Selecionando uma amostra que não gere erro nos chutes iniciais ----------
  repeat{
    amostra <- q(n, ...)
    result <- likelihood_ratio_test(f = f, data = amostra, kicks = kicks,
                                    par0 = par0)
    if (!is.na(result)) break
  }

  q_teorico <- qchisq(p = 1 - sig, df = length(par0[[1]]))
  ifelse(result > q_teorico, 1, 0) # Contando a rejeição.
}
