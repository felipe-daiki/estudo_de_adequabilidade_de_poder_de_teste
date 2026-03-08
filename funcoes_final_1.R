suppressPackageStartupMessages(library(Rlab))
options(warn = -1)
options(repr.plot.width = 16, repr.plot.height = 8)
razao <- seq(from = 0.1, to = 10, by = 0.1)
library(ggplot2)
library(survival)
library(tidyr)
library(survminer)
library(ggsurvfit)

### =============== CRIAÇÃO DO EXEMPLO ===================== ###
exemplo <- function(distribuicao, n, pi,
                    parametro_falha_experimental, parametro_falha_controle,
                    parametro_censura_experimental, parametro_censura_controle,
                    shape_falha_experimental_weibull = NULL, shape_censura_experimental_weibull = NULL,
                    shape_falha_controle_weibull = NULL, shape_censura_controle_weibull = NULL,
                    sdlog_falha_experimental = NULL, sdlog_censura_experimental = NULL,
                    sdlog_falha_controle = NULL, sdlog_censura_controle = NULL){
  
  # seja X a bernoulli para indicar o grupo experimental e grupo controle
  X <- rbern(n,pi)
  n1_inicial <- sum(X == 1)
  n2_inicial <- sum(X == 0)
  
  # ================= GRUPO EXPERIMENTAL ==================== #
  
  if (distribuicao == "weibull") {
    T1 <- rweibull(n1_inicial, shape = shape_falha_experimental_weibull, scale = parametro_falha_experimental)
    C1 <- rweibull(n1_inicial, shape = shape_censura_experimental_weibull, scale = parametro_censura_experimental)
  } 
  
  if (distribuicao == "log_normal") {
    T1 <- rlnorm(n1_inicial, meanlog = parametro_falha_experimental, sdlog = sdlog_falha_experimental)
    C1 <- rlnorm(n1_inicial, meanlog = parametro_censura_experimental, sdlog = sdlog_censura_experimental)
  }
  
  grupo_experimental <- rep('A',n1_inicial)
  
  # lembre que zj = min(Tj, Cj)
  z <- rep(x = 0, n1_inicial)
  for (j in 1:n1_inicial){
    z[j] <- min(T1[j], C1[j])
  }
  
  # e que deltaj = I(Tj <= Cj)
  delta1 <- rep(x = 0, n1_inicial)
  for (j in 1:n1_inicial){
    if (T1[j] <= C1[j]){
      delta1[j] <- 1
    }
  }
  # data frame do grupo experimental
  experimental <- data.frame(z = z, delta = delta1, grupo = grupo_experimental)
  
  # =================== GRUPO CONTROLE ====================== #
  if (distribuicao == "weibull") {
    T2 <- rweibull(n2_inicial, shape = shape_falha_controle_weibull, scale = parametro_falha_controle)
    C2 <- rweibull(n2_inicial, shape = shape_censura_controle_weibull, scale = parametro_censura_controle)
  } 
  
  if (distribuicao == "log_normal") {
    T2 <- rlnorm(n2_inicial, meanlog = parametro_falha_controle, sdlog = sdlog_falha_controle)
    C2 <- rlnorm(n2_inicial, meanlog = parametro_censura_controle, sdlog = sdlog_censura_controle)
  }
  grupo_controle <- rep('B',n2_inicial)
  
  # lembre que zj = min(Tj, Cj)
  z2 <- rep(x = 0, n2_inicial)
  for (j in 1:n2_inicial){
    z2[j] <- min(T2[j], C2[j])
  }
  
  # e que deltaj = I(Tj <= Cj)
  delta2 <- rep(x = 0, n2_inicial)
  for (j in 1:n2_inicial){
    if (T2[j] <= C2[j]){
      delta2[j] <- 1
    }
  }
  # data frame do grupo experimental
  controle <- data.frame(z = z2, delta = delta2, grupo = grupo_controle)
  
  # ================= RELATÓRIO DE CENSURA (NOVO) ================= #
  
  # delta = 1 (evento), delta = 0 (censura)
  # Portanto, a proporção de censura é a média de (delta == 0)
  prop_real_exp <- mean(experimental$delta == 0)
  prop_real_ctrl <- mean(controle$delta == 0)
  
  cat("\n=== Relatório da Simulação ===\n")
  cat(sprintf("Grupo Experimental (A) - Censura Observada: %.1f%%\n", 
              prop_real_exp * 100))
  cat(sprintf("Grupo Controle     (B) - Censura Observada: %.1f%%\n", 
              prop_real_ctrl * 100))
  cat("==============================\n")
  
  # ================= JUNÇÃO E ORDENAÇÃO ====================== #
  dados <- rbind(experimental, controle)
  dados <- dados[order(dados$z), ] # Ordena pelo tempo
  rownames(dados) <- NULL
  rownames(dados) <- NULL
  return(dados)
}

### =============== FUNÇÕES ================================= ###

## numero de individuos do grupo experimental no indice de evento j
n_e <- function(j){
  tamanho_experimental <- sum(dados$grupo == 'A')
  
  for (i in 1:(j)){
    if (i != 1){
      if (dados$grupo[i-1] == 'A'){
        tamanho_experimental <- tamanho_experimental - 1
      } 
    }}
  return(tamanho_experimental)
}

## numero de individuos do grupo controle no indice de evento j
n_c <- function(j){
  tamanho_controle <- sum(dados$grupo == 'B')
  for (i in 1:(j)){
    if (i != 1){
      if (dados$grupo[i-1] == 'B'){
        tamanho_controle <- tamanho_controle - 1
      } 
    }}
  return(tamanho_controle)
}

## numero de eventos do grupo experimental no indice de evento j
d_e <- function(j){
  
  if (dados$grupo[j] == 'A' && dados$delta[j] == 1){
    return(1)
  } else{
    return(0)
  }
}

## numero de eventos do grupo controle no indice de evento j
d_c <- function(j){
  if (dados$grupo[j] == 'B' && dados$delta[j] == 1){
    return(1)
  } else{
    return(0)
  }
}

## funcao densidade experimental no tempo t
f_e <- function(t){
  
  if (distribuicao == "weibull") {
    
    termo1 <- shape_falha_experimental_weibull / parametro_falha_experimental
    termo2 <- (t / parametro_falha_experimental)^(shape_falha_experimental_weibull - 1)
    termo3 <- exp(- (t / parametro_falha_experimental)^shape_falha_experimental_weibull)
    
    valor <- termo1 * termo2 * termo3
    return(valor)
  } 
  
  if (distribuicao == "log_normal") {
    
    termo1 <- 1 / (t * sdlog_falha_experimental * sqrt(2 * pi))
    
    numerador_exp <- (log(t) - parametro_falha_experimental)^2
    denominador_exp <- 2 * (sdlog_falha_experimental^2)
    termo2 <- exp(- (numerador_exp / denominador_exp))
    
    return(termo1 * termo2)
  }
}

## funcao densidade controle no tempo t
f_c <- function(t){
  
  if (distribuicao == "weibull") {
    
    termo1 <- shape_falha_controle_weibull / parametro_falha_controle
    termo2 <- (t / parametro_falha_controle)^(shape_falha_controle_weibull - 1)
    termo3 <- exp(- (t / parametro_falha_controle)^shape_falha_controle_weibull)
    valor <- termo1 * termo2 * termo3
    return(valor)
  } 
  
  if (distribuicao == "log_normal") {
    termo1 <- 1 / (t * sdlog_falha_controle * sqrt(2 * pi))
    
    numerador_exp <- (log(t) - parametro_falha_controle)^2
    denominador_exp <- 2 * (sdlog_falha_controle^2)
    termo2 <- exp(- (numerador_exp / denominador_exp))
    return(termo1 * termo2)
  }
}

## funcao sobrevivencia experimental no tempo t
S_e <- function(t){
  dados_grupoA <- dados[dados$grupo == 'A',]
  s1 <- survfit(Surv(dados_grupoA$z, dados_grupoA$delta)~1)
  if (is.null(summary(s1, times = t)$surv)){
    return(0)
  } else{
    
    return(summary(s1, times = t)$surv)
  }
}

## funcao sobrevivencia controle no tempo t
S_c <- function(t){
  dados_grupoB <- dados[dados$grupo == 'B',]
  s2 <- survfit(Surv(dados_grupoB$z, dados_grupoB$delta)~1)
  valor <- summary(s2, times = t)$surv
  if (is.null(valor)){
    return(0)
  } else{ 
    return(valor)
  }}

## funcao acumulada dos tempos de censura experimental
H_e <- function(t){
  
  if (distribuicao == "weibull") {
    valor <- pweibull(t, shape = shape_censura_experimental_weibull, scale = parametro_censura_experimental)
  }
  if (distribuicao == "log_normal") {
    valor <- plnorm(t, meanlog = parametro_censura_experimental, sdlog = sdlog_censura_experimental)
  }  
  
  return(valor)
}

## funcao acumulada dos tempos de censura controle
H_c <- function(t){
  if (distribuicao == "weibull") { 
    valor <- pweibull(t, shape = shape_censura_controle_weibull, scale = parametro_censura_controle)
  }
  if (distribuicao == "log_normal") {
    valor <- plnorm(t, meanlog = parametro_censura_controle, sdlog = sdlog_censura_controle)
  }  
  
  return(valor)
}

## =========================================================================================================================== ###

### FUNCOES TEORICAS

## Funcao sobrevivencia teorica experimental
S_e_teorica <- function(t) {
  if (distribuicao == "weibull") {
    return(pweibull(t, shape = shape_falha_experimental_weibull, 
                    scale = parametro_falha_experimental, lower.tail = FALSE))
  }
  if (distribuicao == "log_normal") {
    return(plnorm(t, meanlog = parametro_falha_experimental, 
                  sdlog = sdlog_falha_experimental, lower.tail = FALSE))
  }
}

## Funcao sobrevivencia teorica controle
S_c_teorica <- function(t) {
  if (distribuicao == "weibull") {
    return(pweibull(t, shape = shape_falha_controle_weibull, 
                    scale = parametro_falha_controle, lower.tail = FALSE))
  }
  if (distribuicao == "log_normal") {
    return(plnorm(t, meanlog = parametro_falha_controle, 
                  sdlog = sdlog_falha_controle, lower.tail = FALSE))
  }
}

## funcao razao de riscos
log_hazard_ratio <- function(t) {
  
  S_exp  <- S_e_teorica(t)
  S_ctrl <- S_c_teorica(t)
  
  # proteger contra divisão por zero
  S_exp[S_exp <= 1e-12]   <- NA
  S_ctrl[S_ctrl <= 1e-12] <- NA
  
  h_exp  <- f_e(t) / S_exp
  h_ctrl <- f_c(t) / S_ctrl
  
  log_hr <- log(h_exp / h_ctrl)
  
  # substituir não finitos por 0
  log_hr[!is.finite(log_hr)] <- 0
  
  return(log_hr)
}

## funcao que indica o processo de observacao dos eventos
V_teorica <- function(t){
  conta <- (1-pi)*f_c(t)*(1 - H_c(t)) + pi*f_e(t)*(1 - H_e(t))
  return(conta)
}

## funcao de parametro de nao-centralidade
pit_teorica <- function(t){
  conta <- (pi*S_e_teorica(t)*(1-H_e(t)))/((1-pi)*S_c_teorica(t)*(1-H_c(t)) + pi*S_e_teorica(t)*(1-H_e(t)))
  conta <- ifelse(is.na(conta), 0, conta)
  return(conta)
}

## funcao da média
phi_teorica <- function(rho, eta, razao = 0){
  #curva_sobrevivencia <- survfit(Surv(z, delta) ~ 1, data = dados)
  if (razao == 0){
    int_num <- function(t){
      termo_peso <- ((S_c_teorica(t)^rho)*((1-S_c_teorica(t))^eta))
      valor_num <- termo_peso*log_hazard_ratio(t)*pit_teorica(t)*(1-pit_teorica(t))*V_teorica(t)
      return(valor_num)
    }
    numerador <- sqrt(n)*integrate(int_num, lower = 0, upper = 2000, subdivisions = 10000, rel.tol = 1e-4)$value 
    
    int_den <- function(t){
      termo_peso <- (((S_c_teorica(t))^rho)*((1-S_c_teorica(t))^eta))
      valor_den <- (termo_peso^2)*pit_teorica(t)*(1-pit_teorica(t))*V_teorica(t)
      return(valor_den)
    }
    denominador <- integrate(int_den, lower = 0, upper = 2000, subdivisions = 10000, rel.tol = 1e-4)$value
    
  } else{
    
      int_num <- function(t){
        termo_peso <- ((S_c_teorica(t)^rho)*((1-S_c_teorica(t))^eta))
        valor_num <- termo_peso*log(razao)*pit_teorica(t)*(1-pit_teorica(t))*V_teorica(t)
        return(valor_num)
      }
      numerador <- sqrt(n)*integrate(int_num, lower = 0, upper = 2000, subdivisions = 10000, rel.tol = 1e-4)$value 
      
      int_den <- function(t){
        termo_peso <- (((S_c_teorica(t))^rho)*((1-S_c_teorica(t))^eta))
        valor_den <- (termo_peso^2)*pit_teorica(t)*(1-pit_teorica(t))*V_teorica(t)
        return(valor_den)
      }
      denominador <- integrate(int_den, lower = 0, upper = 2000, subdivisions = 10000, rel.tol = 1e-4)$value
      }

  return(numerador/sqrt(denominador))
}

## funcao de poder
poder_teorica <- function(rho,eta, lambda = 0, alpha = 0.05) {
  media <- phi_teorica(rho,eta, lambda)
  z_critico <- qnorm(1 - alpha/2)
  poder <- (1 - pnorm(z_critico, mean = media, sd = 1)) + 
    pnorm(-z_critico, mean = media, sd = 1)
  return(poder)
}

## funcao do grafico de poder do exemplo
grafico_poder_teorica <- function(rho1 = 0, eta1 = 0, rho2 = 1, eta2 = 0, razao = 0, alpha = 0.05){
  vetor_poder_logrank <- c()
  vetor_poder_wilcoxon <- c()
  
  # poder para log-rank
  for (i in razao){
    poder <- poder(rho1,eta1,i, alpha)
    vetor_poder_logrank <- c(vetor_poder_logrank, poder)
  }
  for (j in razao){
    poder <- poder(rho2,eta2,j, alpha)
    vetor_poder_wilcoxon <- c(vetor_poder_wilcoxon, poder)
  }
  library(ggplot2)
  
  # Criar um dataframe combinando os dados de poder
  dados_poder <- data.frame(
    razao = razao,
    logrank = vetor_poder_logrank,
    wilcoxon = vetor_poder_wilcoxon
  )
  
  # Converter para formato longo para ggplot
  dados_poder_long <- reshape2::melt(dados_poder, id.vars = "razao", 
                                     variable.name = "Teste", 
                                     value.name = "Poder")
  
  # Criar o gráfico com ggplot2
  ggplot(data = dados_poder_long, 
         aes(x = razao, y = Poder, color = Teste)) +
    
    geom_line(linewidth = 1.3) +
    geom_point(size = 2, alpha = 0.8) +
    
    geom_hline(yintercept = alpha,
               linetype = "dashed",
               color = "gray40",
               linewidth = 0.8) +
    
    scale_color_manual(values = c("logrank" = "#D55E00",
                                  "wilcoxon" = "#0072B2"),
                       labels = c("Log-Rank", "Wilcoxon")) +
    
    scale_y_continuous(limits = c(0,1),
                       breaks = seq(0,1,0.1)) +
    
    labs(title = "Curvas de Poder dos Testes",
         subtitle = paste("Nível de significância =", alpha),
         x = "Razão de Risco (Hazard Ratio)",
         y = "Poder Estatístico",
         color = "Teste") +
    
    theme_minimal(base_size = 14) +
    
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold",
                                size = 18,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 12,
                                   hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

## simulacao monte carlo
simulacao_final_teorica <- function(n_sim,
                            rho1 = 0, eta1 = 0,
                            rho2 = 1, eta2 = 0, alpha = 0.05) {
  
  nome_teste1 <- paste0("G^(", rho1, ",", eta1, ")")
  nome_teste2 <- paste0("G^(", rho2, ",", eta2, ")")
  
  # matrizes para armazenar resultados
  matriz1 <- matrix(0, nrow = n_sim, ncol = length(razao))
  matriz2 <- matrix(0, nrow = n_sim, ncol = length(razao))
  
  for (i in 1:n_sim) {
    
    cat("Simulação", i, "de", n_sim, "\n")
    
    # gerar nova base
    dados <- exemplo(distribuicao, n, pi,
                      parametro_falha_experimental,
                      parametro_falha_controle,
                      parametro_censura_experimental,
                      parametro_censura_controle,
                      shape_falha_experimental_weibull,
                      shape_censura_experimental_weibull,
                      shape_falha_controle_weibull,
                      shape_censura_controle_weibull,
                      sdlog_falha_experimental,
                      sdlog_censura_experimental,
                      sdlog_falha_controle,
                      sdlog_censura_controle)
    
    for (j in 1:length(razao)) {
      lambda <- razao[j]
      
      matriz1[i, j] <- poder_teorica(rho1, eta1, lambda, alpha)
      matriz2[i, j] <- poder_teorica(rho2, eta2, lambda, alpha)
    }
  }
  
  # médias
  poder_medio1 <- colMeans(matriz1)
  poder_medio2 <- colMeans(matriz2)
  
  # dataframe final
  dados_poder_final <- data.frame(
    razao = razao,
    Teste1 = poder_medio1,
    Teste2 = poder_medio2
  )
  
  # formato longo
  dados_long <- tidyr::pivot_longer(
    dados_poder_final,
    cols = -razao,
    names_to = "Teste",
    values_to = "Poder"
  )
  
  # substituir nomes na legenda
  levels_legenda <- c("Teste1", "Teste2")
  labels_legenda <- c(nome_teste1, nome_teste2)
  
  # gráfico
  ggplot(dados_long,
         aes(x = razao, y = Poder, color = Teste)) +
    
    geom_line(linewidth = 1.3) +
    
    geom_hline(yintercept = alpha,
               linetype = "dashed",
               color = "gray40") +
    
    scale_color_manual(
      values = c("#D55E00", "#0072B2"),
      breaks = levels_legenda,
      labels = labels_legenda
    ) +
    
    scale_y_continuous(limits = c(0,1)) +
    
    labs(
      title = paste("Poder Médio -", n_sim, "Simulações"),
      subtitle = paste("Teste 1:", nome_teste1,
                       "| Teste 2:", nome_teste2),
      x = "Razão de Risco (HR)",
      y = "Poder",
      color = "Teste"
    ) +
    
    theme_minimal(base_size = 14) +
    
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
}

## grafico de poder x tamanho amostral
simulacao_amostral_teorica <- function(n_sim = 100, 
                                       ns = seq(20, 100, by = 20),
                                       rho1 = 0, eta1 = 0,
                                       rho2 = 1, eta2 = 0, 
                                       alpha = 0.05,
                                       razao = 0) { 
  nome_teste1 <- paste0("G^(", rho1, ",", eta1, ")")
  nome_teste2 <- paste0("G^(", rho2, ",", eta2, ")")
  
  matriz1 <- matrix(0, nrow = n_sim, ncol = length(ns))
  matriz2 <- matrix(0, nrow = n_sim, ncol = length(ns))
  
  for (k in 1:length(ns)) {
    n_atual <- ns[k]
    
    cat("Processando n =", n_atual, "\n")
    
    for (i in 1:n_sim) {
      n <<- n_atual 
      dados <<- exemplo(distribuicao, n_atual, pi,
                        parametro_falha_experimental, parametro_falha_controle,
                        parametro_censura_experimental, parametro_censura_controle,
                        shape_falha_experimental_weibull, shape_censura_experimental_weibull,
                        shape_falha_controle_weibull, shape_censura_controle_weibull,
                        sdlog_falha_experimental, sdlog_censura_experimental,
                        sdlog_falha_controle, sdlog_censura_controle)
      # Poder para Teste 1
      matriz1[i, k] <- poder_teorica(rho1, eta1, alpha)
      
      # Poder para Teste 2
      matriz2[i, k] <- poder_teorica(rho2, eta2, alpha)
      cat("mu for n=", n, ":", phi_teorica(rho1, eta1), "\n")
    }
  }
  # Poder médio por n
  poder_medio1 <- colMeans(matriz1)
  poder_medio2 <- colMeans(matriz2)
  
  # Dataframe final
  dados_poder <- data.frame(
    n = ns,
    Teste1 = poder_medio1,
    Teste2 = poder_medio2
  )
  
  # Formato longo para ggplot
  dados_long <- tidyr::pivot_longer(
    dados_poder,
    cols = -n,
    names_to = "Teste",
    values_to = "Poder"
  )
  
  # Substituir nomes na legenda
  levels_legenda <- c("Teste1", "Teste2")
  labels_legenda <- c(nome_teste1, nome_teste2)
  
  # Gráfico
  g <- ggplot(dados_long, aes(x = n, y = Poder, color = Teste)) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 2, alpha = 0.8) +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "blue", linewidth = 0.8) +  # Linha de poder alvo (80%)
    scale_color_manual(
      values = c("#D55E00", "#0072B2"),
      breaks = levels_legenda,
      labels = labels_legenda
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(breaks = ns) +
    labs(
      title = paste("Poder Médio vs. Tamanho Amostral -", n_sim, "Simulações"),
      subtitle = paste("Teste 1:", nome_teste1, "| Teste 2:", nome_teste2),
      x = "Tamanho da Amostra (n)",
      y = "Poder Estatístico",
      color = "Teste"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(g)
  return(dados_poder)  # Retorna dataframe para análise adicional
}

## =========================================================================================================================== ###

### FUNCOES EMPIRICAS

## funcao estatistica de teste
T <- function(rho,eta){
  numerador <- 0
  denominador <- 0
  curva_sobrevivencia <- survfit(Surv(z, delta) ~ 1, data = dados)
  
  indices_eventos <- which(dados$delta == 1)
  
  if (length(indices_eventos) == 0) return(0)
  
  for (x in indices_eventos){
    # captura o tempo t que está na posição x
    tempo_no_indice_x <- dados$z[x]
    
    # usa esse tempo no summary 
    sum_obj <- summary(curva_sobrevivencia, times = tempo_no_indice_x)
    valor_sobrevivencia <- sum_obj$surv
    n <- n_e(x) + n_c(x)
    d <- d_e(x) + d_c(x)
    numerador <- numerador + ((valor_sobrevivencia^rho)*(1-valor_sobrevivencia)^eta)*(d_e(x) - ((n_e(x)*d)/n))
    if ((n - 1) != 0) {
      denominador <- denominador + ((((valor_sobrevivencia^rho)*(1-valor_sobrevivencia)^eta))^2)*((n_e(x) * n_c(x) * d * (n - d)) / ((n^2) * (n - 1)))
      
    }
  }  
  return(numerador/sqrt(denominador)) 
}

phi_empirica <- function(rho, eta, razao = 0){ # razao > 0 para PH constante
  
  # ajuste da KM para os pesos
  curva_km <- survfit(Surv(z, delta) ~ 1, data = dados)
  
  # identificar índices de eventos
  indices_eventos <- which(dados$delta == 1)
  
  if (length(indices_eventos) == 0) return(0)
  
  soma_numerador <- 0
  soma_denominador <- 0
  
  for (i in indices_eventos){
    

    ti <- dados$z[i]
    s_t <- summary(curva_km, times = ti, extend = TRUE)$surv
    peso <- (s_t^rho) * ((1 - s_t)^eta)
    
    n_exp <- n_e(i)
    n_ctrl <- n_c(i)
    n_total_risco <- n_exp + n_ctrl
    
    if (n_total_risco > 1) {
      p_t <- n_exp / n_total_risco
      
      var_termo <- p_t * (1 - p_t)
      
      if (razao == 0) {
        beta_t <- log_hazard_ratio(ti) 
      } else {
        beta_t <- log(razao)
      }
      
      soma_numerador <- soma_numerador + (peso * beta_t * var_termo)
      
      soma_denominador <- soma_denominador + (peso^2 * var_termo)
    }
  }
  
  if (soma_denominador <= 0) return(0)
  

  return(soma_numerador / sqrt(soma_denominador))
}

## funcao de poder empirico 
poder_empirica <- function(rho, eta, razao = 0, alpha = 0.05) {

  estatistica <- T(rho, eta)
  z_critico <- qnorm(1 - alpha/2)
  
  return(ifelse(abs(estatistica) > z_critico, 1, 0))
}

## Simulacao Poder Empirico vs Razao de Risco
simulacao_final_empirica <- function(n_sim,
                                     rho1 = 0, eta1 = 0,
                                     rho2 = 1, eta2 = 0, alpha = 0.05){
  
  nome_teste1 <- paste0("G^(", rho1, ",", eta1, ")")
  nome_teste2 <- paste0("G^(", rho2, ",", eta2, ")")
  
  # matrizes para armazenar resultados de rejeicao (0 ou 1)
  matriz1 <- matrix(0, nrow = n_sim, ncol = length(razao))
  matriz2 <- matrix(0, nrow = n_sim, ncol = length(razao))
  
  for (j in 1:length(razao)) {
    lambda_atual <- razao[j]
    
    for (i in 1:n_sim) {

      dados <<- exemplo(distribuicao, n, pi,
                        parametro_falha_experimental = parametro_falha_controle / lambda_atual,             ## SERVE APENAS PARA O CONTEXTO DE RISCOS PROPORCIONAIS EXPONENCIAIS
                        parametro_falha_controle,
                        parametro_censura_experimental, parametro_censura_controle,
                        shape_falha_experimental_weibull, shape_censura_experimental_weibull,
                        shape_falha_controle_weibull, shape_censura_controle_weibull,
                        sdlog_falha_experimental, sdlog_censura_experimental,
                        sdlog_falha_controle, sdlog_censura_controle)
      
      matriz1[i, j] <- poder_empirica(rho1, eta1, lambda_atual, alpha)
      matriz2[i, j] <- poder_empirica(rho2, eta2, lambda_atual, alpha)
    }
  }
  poder_medio1 <- colMeans(matriz1)
  poder_medio2 <- colMeans(matriz2)
  
  dados_poder_final <- data.frame(razao = razao, Teste1 = poder_medio1, Teste2 = poder_medio2)
  dados_long <- tidyr::pivot_longer(dados_poder_final, cols = -razao, names_to = "Teste", values_to = "Poder")
  
  # gerar o grafico
  ggplot(dados_long, aes(x = razao, y = Poder, color = Teste)) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 2) +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") + # Linha de 80% de poder
    scale_color_manual(values = c("#D55E00", "#0072B2"), labels = c(nome_teste1, nome_teste2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(title = paste("Poder Empírico (Monte Carlo) -", n_sim, "Simulações"),
         subtitle = paste("n =", n, "| Testes:", nome_teste1, "vs", nome_teste2),
         x = "Razão de Risco (Hazard Ratio)", y = "Poder (Taxa de Rejeição)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))
}

simulacao_amostral_empirica <- function(n_sim = 200, 
                                        ns = seq(20, 100, by = 20),
                                        rho1 = 0, eta1 = 0,
                                        rho2 = 1, eta2 = 0, 
                                        alpha = 0.05) {
  
  resultados <- data.frame()
  z_critico <- qnorm(1 - alpha/2)
  
  for (n_atual in ns) {
    cat("Processando n =", n_atual, "\n")
    rejeicoes1 <- 0
    rejeicoes2 <- 0
    n <<- n_atual
    for (i in 1:n_sim) {

      dados <<- exemplo(distribuicao, n, pi, 
                        parametro_falha_experimental, parametro_falha_controle, 
                        parametro_censura_experimental, parametro_censura_controle,
                        shape_falha_experimental_weibull, shape_censura_experimental_weibull, 
                        shape_falha_controle_weibull, shape_censura_controle_weibull, 
                        sdlog_falha_experimental, sdlog_censura_experimental, 
                        sdlog_falha_controle, sdlog_censura_controle)
      
      # teste empirico
      if (abs(T(rho1, eta1)) > z_critico) rejeicoes1 <- rejeicoes1 + 1
      if (abs(T(rho2, eta2)) > z_critico) rejeicoes2 <- rejeicoes2 + 1
    }
    
    resultados <- rbind(resultados, data.frame(
      n = n_atual,
      Poder = rejeicoes1 / n_sim,
      Teste = paste0("G(", rho1, ",", eta1, ")")
    ))
    resultados <- rbind(resultados, data.frame(
      n = n_atual,
      Poder = rejeicoes2 / n_sim,
      Teste = paste0("G(", rho2, ",", eta2, ")")
    ))
  }
  
  # Gráfico de Poder Empírico
  ggplot(resultados, aes(x = n, y = Poder, color = Teste)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(title = "Poder Empírico: Log-Rank vs Wilcoxon",
         subtitle = paste("Simulação de Monte Carlo (", n_sim, " repetições)"),
         x = "Tamanho Amostral (n)", y = "Poder (Taxa de Rejeição)") +
    theme_minimal()
}

grafico_phi_comparativo <- function(razao, rho = 0, eta = 0){
  
  phi_emp <- numeric(length(razao))
  phi_teo <- numeric(length(razao))
  
  # calcular as medias
  for (i in 1:length(razao)){
    phi_emp[i] <- phi_empirica(rho = rho, eta = eta, razao = razao[i])
    phi_teo[i] <- phi_teorica(rho = rho, eta = eta, razao = razao[i])
  }
  
  library(ggplot2)
  
  # criacao dos dadso
  dados_phi <- data.frame(
    HR = rep(razao, 2),
    Phi = c(phi_emp, phi_teo),
    Tipo = rep(c("Empírica", "Teórica"), each = length(razao)) # repete os pontos
  )
  
  g <- ggplot(dados_phi,
              aes(x = HR, y = Phi, color = Tipo, linetype = Tipo)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "gray50") +
    scale_color_manual(
      values = c("Empírica" = "#D55E00", 
                 "Teórica"  = "#0072B2"),
      labels = c("φ Empírica", "φ Teórica")
    ) +
    scale_linetype_manual(
      values = c("Empírica" = "solid",
                 "Teórica"  = "dashed"),
      labels = c("φ Empírica", "φ Teórica")
    ) +
    labs(
      title = paste0("Comparação entre φ Empírica e Teórica (ρ = ",
                     rho, ", η = ", eta, ")"),
      x = "Razão de Risco (HR)",
      y = "Valor de φ"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  print(g)
}

