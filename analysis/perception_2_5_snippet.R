# ============================================================================
# PERCEPTION 2.5 - ANÁLISE DE DADOS CCLE, PRISM E GDSC (CORRIGIDO)
# ============================================================================

# Configuração inicial
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║                 PERCEPTION 2.5 - ANÁLISE                ║\n")
cat("║          Dados CCLE, PRISM e GDSC                      ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Caminhos dos arquivos
ccle_expr_path <- "CCLE_expression.csv"       # Expressão bulk das linhagens celulares (log2 TPM)
prism_resp_path <- "PRISM_drug_response.csv"   # Resposta a fármacos PRISM (AUC ou IC50)
gdsc_resp_path <- "GDSC_drug_response.csv"     # Resposta a fármacos GDSC

# ============================================================================
# FUNÇÕES AUXILIARES CORRIGIDAS
# ============================================================================

# Função para carregar arquivos CSV com possíveis row.names duplicados
safe_read_csv <- function(file_path, has_rownames = TRUE) {
  cat(sprintf("  Lendo: %s\n", basename(file_path)))
  
  # Primeiro, ler sem row.names para verificar
  temp_data <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Verificar se há coluna com nomes das linhas
  if (has_rownames) {
    # Verificar se há duplicatas na primeira coluna
    first_col <- temp_data[, 1]
    duplicates <- duplicated(first_col)
    
    if (any(duplicates)) {
      warning(sprintf("  ⚠️  Encontradas %d linhas duplicadas em %s. Removendo duplicatas...", 
                      sum(duplicates), basename(file_path)))
      
      # Manter apenas a primeira ocorrência
      temp_data <- temp_data[!duplicated(first_col), ]
      first_col <- first_col[!duplicated(first_col)]
    }
    
    # Definir row.names
    rownames(temp_data) <- first_col
    temp_data <- temp_data[, -1, drop = FALSE]  # Remover primeira coluna
  }
  
  cat(sprintf("    • %d linhas, %d colunas\n", nrow(temp_data), ncol(temp_data)))
  return(temp_data)
}

# Função para carregar dados de expressão com tratamento especial
load_expression_data <- function(file_path) {
  cat("  1. Carregando expressão gênica CCLE...\n")
  
  # Ler dados de expressão
  expr_data <- safe_read_csv(file_path, has_rownames = TRUE)
  
  # Verificar se os dados estão em log2 TPM
  # Se não estiverem, aplicar log2 transformação (adicionando 1 para evitar log(0))
  if (max(expr_data, na.rm = TRUE) > 100) {
    cat("     • Aplicando transformação log2(TPM + 1)...\n")
    expr_data <- log2(expr_data + 1)
  }
  
  # Converter para matriz numérica
  expr_matrix <- as.matrix(expr_data)
  
  # Remover genes com todos os valores iguais a zero ou NA
  zero_genes <- colSums(expr_matrix == 0, na.rm = TRUE) == nrow(expr_matrix)
  if (any(zero_genes)) {
    cat(sprintf("     • Removendo %d genes com todos os valores zero\n", sum(zero_genes)))
    expr_matrix <- expr_matrix[, !zero_genes, drop = FALSE]
  }
  
  cat(sprintf("     • %d linhagens celulares, %d genes\n", 
              nrow(expr_matrix), ncol(expr_matrix)))
  
  return(expr_matrix)
}

# Função para carregar dados de resposta a fármacos
load_response_data <- function(file_path, dataset_name = "") {
  cat(sprintf("  2. Carregando dados %s...\n", dataset_name))
  
  # Ler dados de resposta
  resp_data <- safe_read_csv(file_path, has_rownames = TRUE)
  
  # Verificar tipo de dados (AUC ou IC50)
  # AUC: valores entre 0-1, onde baixo = sensível, alto = resistente
  # IC50: valores > 0, onde baixo = sensível, alto = resistente
  sample_values <- as.numeric(as.matrix(resp_data))
  sample_values <- sample_values[!is.na(sample_values) & is.finite(sample_values)]
  
  if (length(sample_values) > 0) {
    if (max(sample_values) <= 1 && min(sample_values) >= 0) {
      cat("     • Dados detectados como AUC (0-1)\n")
    } else if (max(sample_values) > 1) {
      cat("     • Dados detectados como IC50 (valores positivos)\n")
    }
  }
  
  cat(sprintf("     • %d linhagens, %d fármacos\n", 
              nrow(resp_data), ncol(resp_data)))
  
  return(resp_data)
}

# ============================================================================
# FUNÇÃO PARA CARREGAR E INTEGRAR DADOS (CORRIGIDA)
# ============================================================================

load_and_integrate_datasets <- function(expression_file, prism_file, gdsc_file, 
                                        max_genes = 5000, min_samples = 20) {
  
  cat("Carregando e integrando dados...\n")
  
  # 1. Carregar dados de expressão CCLE
  ccle_expr <- load_expression_data(expression_file)
  
  # 2. Carregar dados PRISM
  prism_data <- load_response_data(prism_file, "PRISM")
  
  # 3. Carregar dados GDSC (se o arquivo existir)
  if (file.exists(gdsc_file)) {
    gdsc_data <- load_response_data(gdsc_file, "GDSC")
  } else {
    cat("  3. Arquivo GDSC não encontrado, usando apenas PRISM...\n")
    gdsc_data <- NULL
  }
  
  # 4. Encontrar linhagens celulares comuns
  if (!is.null(gdsc_data)) {
    common_cell_lines <- Reduce(intersect, list(
      rownames(ccle_expr),
      rownames(prism_data),
      rownames(gdsc_data)
    ))
  } else {
    common_cell_lines <- intersect(
      rownames(ccle_expr),
      rownames(prism_data)
    )
  }
  
  cat(sprintf("     • Linhagens comuns: %d\n", length(common_cell_lines)))
  
  if (length(common_cell_lines) < min_samples) {
    warning(sprintf("Número baixo de linhagens comuns: %d (mínimo recomendado: %d)", 
                    length(common_cell_lines), min_samples))
  }
  
  # 5. Filtrar dados para linhagens comuns
  ccle_expr <- ccle_expr[common_cell_lines, , drop = FALSE]
  prism_data <- prism_data[common_cell_lines, , drop = FALSE]
  
  if (!is.null(gdsc_data)) {
    gdsc_data <- gdsc_data[common_cell_lines, , drop = FALSE]
  }
  
  # 6. Selecionar genes mais variáveis
  cat("  4. Selecionando genes mais variáveis...\n")
  
  # Calcular variância por gene
  gene_variance <- apply(ccle_expr, 2, var, na.rm = TRUE)
  
  # Ordenar por variância (decrescente)
  sorted_genes <- names(sort(gene_variance, decreasing = TRUE))
  
  # Selecionar top genes
  n_genes_to_keep <- min(max_genes, length(sorted_genes))
  selected_genes <- sorted_genes[1:n_genes_to_keep]
  
  ccle_expr <- ccle_expr[, selected_genes, drop = FALSE]
  
  cat(sprintf("     • %d genes selecionados (variância mais alta)\n", n_genes_to_keep))
  
  # 7. Processar dados de resposta
  cat("  5. Processando dados de resposta...\n")
  
  # Função para processar matriz de resposta
  process_response_matrix <- function(resp_matrix) {
    # Identificar tipo de dados
    sample_vals <- as.numeric(as.matrix(resp_matrix))
    sample_vals <- sample_vals[!is.na(sample_vals) & is.finite(sample_vals)]
    
    if (length(sample_vals) == 0) return(resp_matrix)
    
    # Se for AUC (0-1), converter para sensibilidade
    if (max(sample_vals) <= 1 && min(sample_vals) >= 0) {
      cat("     • Convertendo AUC para sensibilidade (1 - AUC)\n")
      resp_matrix <- 1 - resp_matrix
    } 
    # Se for IC50, aplicar log transformação
    else if (max(sample_vals) > 10) {
      cat("     • Aplicando transformação log10(IC50 + 1)\n")
      resp_matrix <- log10(resp_matrix + 1)
    }
    
    return(resp_matrix)
  }
  
  # Processar PRISM
  prism_processed <- process_response_matrix(prism_data)
  
  # Processar GDSC se existir
  if (!is.null(gdsc_data)) {
    gdsc_processed <- process_response_matrix(gdsc_data)
  } else {
    gdsc_processed <- NULL
  }
  
  # 8. Combinar dados de resposta (se ambos existirem)
  if (!is.null(gdsc_processed)) {
    cat("  6. Combinando dados de resposta PRISM e GDSC...\n")
    
    # Encontrar fármacos únicos
    prism_drugs <- colnames(prism_processed)
    gdsc_drugs <- colnames(gdsc_processed)
    
    common_drugs <- intersect(prism_drugs, gdsc_drugs)
    unique_prism_drugs <- setdiff(prism_drugs, gdsc_drugs)
    unique_gdsc_drugs <- setdiff(gdsc_drugs, prism_drugs)
    
    cat(sprintf("     • Fármacos comuns: %d\n", length(common_drugs)))
    cat(sprintf("     • Fármacos únicos PRISM: %d\n", length(unique_prism_drugs)))
    cat(sprintf("     • Fármacos únicos GDSC: %d\n", length(unique_gdsc_drugs)))
    
    # Para fármacos comuns, usar média de PRISM e GDSC
    combined_response <- matrix(NA, 
                                nrow = length(common_cell_lines), 
                                ncol = length(prism_drugs) + length(unique_gdsc_drugs))
    
    rownames(combined_response) <- common_cell_lines
    colnames(combined_response) <- c(prism_drugs, unique_gdsc_drugs)
    
    # Preencher com dados PRISM
    combined_response[, prism_drugs] <- as.matrix(prism_processed[, prism_drugs, drop = FALSE])
    
    # Para fármacos comuns, calcular média
    if (length(common_drugs) > 0) {
      for (drug in common_drugs) {
        combined_response[, drug] <- rowMeans(
          cbind(prism_processed[, drug], gdsc_processed[, drug]),
          na.rm = TRUE
        )
      }
    }
    
    # Adicionar fármacos únicos do GDSC
    if (length(unique_gdsc_drugs) > 0) {
      combined_response[, unique_gdsc_drugs] <- as.matrix(gdsc_processed[, unique_gdsc_drugs, drop = FALSE])
    }
    
    response_data <- combined_response
    
  } else {
    # Usar apenas PRISM
    cat("  6. Usando apenas dados PRISM (GDSC não disponível)...\n")
    response_data <- as.matrix(prism_processed)
  }
  
  # 9. Filtrar fármacos com muitos valores missing
  cat("  7. Filtrando fármacos...\n")
  
  missing_percent <- colSums(is.na(response_data)) / nrow(response_data) * 100
  valid_drugs <- names(missing_percent[missing_percent < 50])
  
  if (length(valid_drugs) > 0) {
    response_data <- response_data[, valid_drugs, drop = FALSE]
    cat(sprintf("     • Fármacos válidos (missing < 50%%): %d/%d\n", 
                length(valid_drugs), length(missing_percent)))
  } else {
    stop("Nenhum fármaco válido após filtragem (todos têm mais de 50% valores missing)")
  }
  
  # 10. Normalizar expressão gênica
  cat("  8. Normalizando dados de expressão...\n")
  
  # Aplicar normalização Z-score por gene
  ccle_expr_norm <- scale(ccle_expr, center = TRUE, scale = TRUE)
  
  # Verificar se há NAs resultantes da normalização
  na_count <- sum(is.na(ccle_expr_norm))
  if (na_count > 0) {
    warning(sprintf("  ⚠️  %d valores NA após normalização. Substituindo por 0.", na_count))
    ccle_expr_norm[is.na(ccle_expr_norm)] <- 0
  }
  
  # 11. Criar objeto de dados integrado
  integrated_data <- list(
    expression = ccle_expr_norm,
    response = response_data,
    cell_lines = common_cell_lines,
    genes = colnames(ccle_expr_norm),
    drugs = colnames(response_data),
    metadata = list(
      n_cell_lines = length(common_cell_lines),
      n_genes = ncol(ccle_expr_norm),
      n_drugs = ncol(response_data),
      genes_selected_method = "high_variance",
      response_type = ifelse(max(response_data, na.rm = TRUE) <= 1, "AUC-based", "IC50-based")
    )
  )
  
  cat("  ✅ Dados integrados com sucesso!\n")
  cat(sprintf("  • Linhagens finais: %d\n", integrated_data$metadata$n_cell_lines))
  cat(sprintf("  • Genes finais: %d\n", integrated_data$metadata$n_genes))
  cat(sprintf("  • Fármacos finais: %d\n", integrated_data$metadata$n_drugs))
  cat("\n")
  
  return(integrated_data)
}

# ============================================================================
# FUNÇÃO PRINCIPAL SIMPLIFICADA
# ============================================================================

run_perception_analysis <- function(expression_file, response_file, 
                                    output_dir = "perception_results",
                                    max_drugs = 50,
                                    n_cores = 1) {
  
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║                    PERCEPTION 2.5                       ║\n")
  cat("║   Sistema de Predição de Resistência a Drogas           ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")
  
  # Configurar ambiente
  cat("Configurando ambiente PERCEPTION 2.5...\n")
  
  required_packages <- c(
    "glmnet",    # Modelos regularizados
    "caret",     # Machine learning
    "dplyr",     # Manipulação de dados
    "ggplot2",   # Visualizações
    "pROC",      # Curvas ROC
    "ROCR",      # Avaliação de modelos
    "reshape2",  # Transformação de dados
    "corrplot",  # Matrizes de correlação
    "viridis",   # Paleta de cores
    "patchwork"  # Combinação de gráficos
  )
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Instalando pacote: %s\n", pkg))
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
      library(pkg, character.only = TRUE)
      cat(sprintf("  ✅ %s instalado\n", pkg))
    } else {
      cat(sprintf("  ✅ %s carregado\n", pkg))
    }
  }
  
  cat("\n✅ Ambiente configurado com sucesso!\n")
  
  # Carregar dados
  cat("\nCarregando dados...\n")
  
  # Carregar expressão gênica
  cat("  • Expressão gênica...\n")
  expr_data <- load_expression_data(expression_file)
  
  # Carregar resposta a fármacos
  cat("  • Resposta a fármacos...\n")
  resp_data <- load_response_data(response_file, "PRISM")
  
  # Encontrar amostras comuns
  common_samples <- intersect(rownames(expr_data), rownames(resp_data))
  
  if (length(common_samples) < 20) {
    stop(sprintf("Número insuficiente de amostras comuns: %d (mínimo: 20)", length(common_samples)))
  }
  
  cat(sprintf("  • Amostras comuns: %d\n", length(common_samples)))
  
  # Filtrar dados
  expr_data <- expr_data[common_samples, , drop = FALSE]
  resp_data <- resp_data[common_samples, , drop = FALSE]
  
  # Selecionar genes mais variáveis (5000 genes)
  cat("  • Selecionando genes mais variáveis...\n")
  gene_variance <- apply(expr_data, 2, var, na.rm = TRUE)
  top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:min(5000, length(gene_variance))]
  expr_data <- expr_data[, top_genes, drop = FALSE]
  
  # Normalizar expressão
  cat("  • Normalizando dados de expressão...\n")
  expr_norm <- scale(expr_data, center = TRUE, scale = TRUE)
  expr_norm[is.na(expr_norm)] <- 0
  
  # Processar dados de resposta
  cat("  • Processando dados de resposta...\n")
  
  # Verificar se é AUC ou IC50
  sample_vals <- as.numeric(as.matrix(resp_data))
  sample_vals <- sample_vals[!is.na(sample_vals) & is.finite(sample_vals)]
  
  if (length(sample_vals) > 0) {
    if (max(sample_vals) <= 1 && min(sample_vals) >= 0) {
      cat("    • Convertendo AUC para sensibilidade (1 - AUC)\n")
      resp_data <- 1 - resp_data
    } else if (max(sample_vals) > 1) {
      cat("    • Aplicando transformação log10\n")
      resp_data <- log10(resp_data + 1)
    }
  }
  
  # Limitar número de fármacos se necessário
  if (!is.null(max_drugs) && max_drugs < ncol(resp_data)) {
    cat(sprintf("  • Limitando a %d fármacos\n", max_drugs))
    resp_data <- resp_data[, 1:max_drugs, drop = FALSE]
  }
  
  # Filtrar fármacos com muitos valores missing
  missing_percent <- colSums(is.na(resp_data)) / nrow(resp_data) * 100
  valid_drugs <- names(missing_percent[missing_percent < 50])
  
  if (length(valid_drugs) == 0) {
    stop("Nenhum fármaco válido após filtragem (todos têm mais de 50% valores missing)")
  }
  
  resp_data <- resp_data[, valid_drugs, drop = FALSE]
  cat(sprintf("  • Fármacos válidos: %d\n", length(valid_drugs)))
  
  # Preparar dados para PERCEPTION
  data_for_perception <- list(
    expression = expr_norm,
    response = resp_data,
    metadata = NULL,
    samples = common_samples,
    genes = colnames(expr_norm),
    drugs = colnames(resp_data)
  )
  
  # Criar diretório de saída
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Executar pipeline PERCEPTION
  cat("\n" + paste(rep("=", 70), collapse = "") + "\n")
  cat("EXECUTANDO PIPELINE PERCEPTION 2.5\n")
  cat(paste(rep("=", 70), collapse = "") + "\n\n")
  
  # Carregar funções PERCEPTION (do script original)
  # Nota: As funções train_drug_model, run_perception_pipeline, etc.
  # devem estar disponíveis no ambiente
  
  # Executar pipeline
  results <- run_perception_pipeline(
    data = data_for_perception,
    output_dir = output_dir,
    parallel = (n_cores > 1),
    n_cores = n_cores
  )
  
  cat("\n" + paste(rep("=", 70), collapse = "") + "\n")
  cat("✅ ANÁLISE CONCLUÍDA!\n")
  cat(paste(rep("=", 70), collapse = "") + "\n")
  cat(sprintf("Resultados salvos em: %s\n", output_dir))
  
  return(results)
}

# ============================================================================
# EXECUÇÃO PASSO A PASSO
# ============================================================================

# Primeiro, carregar as funções do PERCEPTION 2.5 que você já tem
# (certifique-se de que todas as funções estão no ambiente)

# 1. Testar carregamento dos dados
cat("Testando carregamento dos dados...\n")

# Verificar se arquivos existem
if (!file.exists(ccle_expr_path)) {
  stop(sprintf("Arquivo não encontrado: %s", ccle_expr_path))
}

if (!file.exists(prism_resp_path)) {
  stop(sprintf("Arquivo não encontrado: %s", prism_resp_path))
}

# 2. Carregar e integrar dados (versão simplificada)
cat("\n1. Carregando dados de expressão...\n")
ccle_expr <- load_expression_data(ccle_expr_path)

cat("\n2. Carregando dados de resposta PRISM...\n")
prism_resp <- load_response_data(prism_resp_path, "PRISM")

# 3. Preparar dados para análise
cat("\n3. Preparando dados para análise...\n")

# Encontrar amostras comuns
common_cell_lines <- intersect(rownames(ccle_expr), rownames(prism_resp))
cat(sprintf("   • Linhagens comuns: %d\n", length(common_cell_lines)))

# Filtrar dados
ccle_expr_filtered <- ccle_expr[common_cell_lines, ]
prism_resp_filtered <- prism_resp[common_cell_lines, ]

# Selecionar 5000 genes mais variáveis
cat("   • Selecionando genes mais variáveis...\n")
gene_variance <- apply(ccle_expr_filtered, 2, var, na.rm = TRUE)
top_5000_genes <- names(sort(gene_variance, decreasing = TRUE))[1:5000]
ccle_expr_filtered <- ccle_expr_filtered

# Normalizar
cat("   • Normalizando expressão gênica...\n")
ccle_expr_norm <- scale(ccle_expr_filtered, center = TRUE, scale = TRUE)
ccle_expr_norm[is.na(ccle_expr_norm)] <- 0

# Processar resposta (assumindo que é AUC)
cat("   • Processando dados de resposta...\n")
prism_resp_processed <- 1 - prism_resp_filtered  # Converter AUC para sensibilidade

# Filtrar fármacos com menos de 50% missing
missing_percent <- colSums(is.na(prism_resp_processed)) / nrow(prism_resp_processed) * 100
valid_drugs <- names(missing_percent[missing_percent < 50])

if (length(valid_drugs) == 0) {
  stop("Nenhum fármaco válido após filtragem")
}

prism_resp_final <- prism_resp_processed[, valid_drugs, drop = FALSE]

cat(sprintf("   • Dados finais: %d linhagens, %d genes, %d fármacos\n",
            nrow(ccle_expr_norm), ncol(ccle_expr_norm), ncol(prism_resp_final)))

# 4. Executar PERCEPTION com um subconjunto para teste
cat("\n4. Executando PERCEPTION (modo teste com 5 fármacos)...\n")

# Criar objeto de dados
test_data <- list(
  expression = ccle_expr_norm,
  response = prism_resp_final[, 1:min(5, ncol(prism_resp_final)), drop = FALSE],
  metadata = NULL,
  samples = rownames(ccle_expr_norm),
  genes = colnames(ccle_expr_norm),
  drugs = colnames(prism_resp_final[, 1:min(5, ncol(prism_resp_final)), drop = FALSE])
)

# Executar pipeline
test_results <- run_perception_pipeline(
  data = test_data,
  output_dir = "test_results",
  parallel = FALSE,
  max_drugs = 5
)

# 5. Se funcionar, executar análise completa
cat("\n5. Preparando para análise completa...\n")

full_data <- list(
  expression = ccle_expr_norm,
  response = prism_resp_final,
  metadata = NULL,
  samples = rownames(ccle_expr_norm),
  genes = colnames(ccle_expr_norm),
  drugs = colnames(prism_resp_final)
)

# Perguntar se deve executar análise completa
run_full <- readline("Executar análise completa com todos os fármacos? (s/n): ")

if (tolower(run_full) == "s") {
  cat("\n6. Executando análise completa...\n")
  
  full_results <- run_perception_pipeline(
    data = full_data,
    output_dir = "full_analysis_results",
    parallel = TRUE,
    n_cores = parallel::detectCores() - 1,
    max_drugs = min(50, ncol(prism_resp_final))
  )
  
  cat("\n" + paste(rep("=", 70), collapse = "") + "\n")
  cat("✅ ANÁLISE COMPLETA CONCLUÍDA!\n")
  cat(paste(rep("=", 70), collapse = "") + "\n")
  
  # Resumo
  summary_df <- create_summary_dataframe(full_results$models)
  cat(sprintf("\nResumo da análise:\n"))
  cat(sprintf("  • Modelos treinados: %d\n", nrow(summary_df)))
  cat(sprintf("  • Correlação média (teste): %.3f\n", mean(summary_df$cor_test, na.rm = TRUE)))
  cat(sprintf("  • Correlação mediana (teste): %.3f\n", median(summary_df$cor_test, na.rm = TRUE)))
  cat(sprintf("  • Top 3 fármacos melhor preditos:\n"))
  
  top_3 <- head(summary_df[order(-summary_df$cor_test), ], 3)
  for (i in 1:nrow(top_3)) {
    cat(sprintf("      %d. %s: r=%.3f\n", i, top_3$drug[i], top_3$cor_test[i]))
  }
  
} else {
  cat("\nAnálise interrompida pelo usuário.\n")
}

# ============================================================================
# FUNÇÃO DE RESOLUÇÃO DE PROBLEMAS
# ============================================================================

troubleshoot_data_loading <- function() {
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║                 DIAGNÓSTICO DE PROBLEMAS                ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n\n")
  
  # Verificar arquivos
  cat("1. Verificando arquivos...\n")
  files <- c(ccle_expr_path, prism_resp_path, gdsc_resp_path)
  for (f in files) {
    if (file.exists(f)) {
      cat(sprintf("   ✅ %s: encontrado (%.1f MB)\n", f, file.size(f)/1024/1024))
    } else {
      cat(sprintf("   ❌ %s: NÃO ENCONTRADO\n", f))
    }
  }
  
  # Verificar conteúdo dos arquivos
  cat("\n2. Examinando arquivo PRISM...\n")
  
  # Ler primeiras linhas
  prism_head <- readLines(prism_resp_path, n = 5)
  cat("   Primeiras linhas:\n")
  for (i in 1:length(prism_head)) {
    cat(sprintf("   [%d] %s\n", i, prism_head[i]))
  }
  
  # Contar colunas
  first_line <- prism_head[1]
  n_commas <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  cat(sprintf("   • Número de colunas (pela primeira linha): %d\n", n_commas + 1))
  
  # Verificar duplicatas na primeira coluna
  cat("\n3. Verificando duplicatas...\n")
  
  # Ler apenas a primeira coluna
  temp <- read.csv(prism_resp_path, check.names = FALSE, stringsAsFactors = FALSE)
  first_col <- temp[, 1]
  duplicates <- duplicated(first_col)
  
  if (any(duplicates)) {
    cat(sprintf("   ⚠️  Encontradas %d linhas duplicadas\n", sum(duplicates)))
    cat("   Exemplos de duplicatas:\n")
    dup_values <- unique(first_col[duplicates])
    for (i in 1:min(5, length(dup_values))) {
      cat(sprintf("      • %s\n", dup_values[i]))
    }
  } else {
    cat("   ✅ Nenhuma duplicata encontrada na primeira coluna\n")
  }
  
  # Sugestões
  cat("\n4. SUGESTÕES:\n")
  cat("   • Se há duplicatas, remover linhas duplicadas manualmente\n")
  cat("   • Verificar se o arquivo tem cabeçalho correto\n")
  cat("   • Verificar codificação do arquivo (UTF-8 recomendado)\n")
  cat("   • Usar a função safe_read_csv() para carregar dados\n")
}

# Para executar diagnóstico:
# troubleshoot_data_loading()

# ============================================================================
# EXECUTAR ANÁLISE DIRETAMENTE (DESCOMENTE SE NECESSÁRIO)
# ============================================================================

# Opção 1: Executar diagnóstico primeiro
# troubleshoot_data_loading()

# Opção 2: Executar análise diretamente
results <- run_perception_analysis(
#   expression_file = ccle_expr_path,
#   response_file = prism_resp_path,
#   output_dir = "prism_analysis",
#   max_drugs = 30,
#   n_cores = 4
# )

cat("\n" + paste(rep("=", 70), collapse = "") + "\n")
cat("SCRIPT PRONTO PARA USO\n")
cat(paste(rep("=", 70), collapse = "") + "\n")
cat("\nPara executar a análise, use:\n")
cat("  results <- run_perception_analysis(ccle_expr_path, prism_resp_path)\n")
cat("\nPara diagnóstico de problemas:\n")
cat("  troubleshoot_data_loading()\n")
