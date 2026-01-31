# ============================================================================
# PERCEPTION 2.7 (PRO) - Pipeline robusto para prever resposta/resistência a drogas
# - Sem leakage: scaling apenas no treino e aplicado no teste/novos dados
# - Split estratificado por quantis + múltiplas tentativas para evitar teste sem variância
# - Filtro de genes por top variância (ex.: 3000)
# - Filtro de drogas por amostras/variância/valores únicos
# - Um modelo por droga (glmnet) com busca de alpha + lambda via CV
# - Redução de dimensionalidade por top_k_genes (|coef|) após ajuste
# - Paralelo opcional (doParallel + foreach)
# ============================================================================

# -----------------------------
# 0) Utilitários: logging e erros
# -----------------------------
.log <- function(level = "INFO", ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %-5s | %s\n", ts, level, paste0(..., collapse = "")))
}
stopf <- function(...) stop(sprintf(...), call. = FALSE)
warnf <- function(...) warning(sprintf(...), call. = FALSE)

.safe_cor <- function(a, b, method = "pearson") {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = method))
}

# -----------------------------
# 1) Setup de dependências
# -----------------------------
setup_perception <- function(install_missing = TRUE) {
  .log("INFO", "Configurando ambiente PERCEPTION 2.7...")
  required_packages <- c(
    "glmnet", "caret", "dplyr", "data.table"
  )
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!install_missing) stopf("Pacote necessário não instalado: %s", pkg)
      .log("INFO", "Instalando pacote: ", pkg)
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  .log("INFO", "✅ Ambiente pronto.")
  invisible(TRUE)
}

# -----------------------------
# 2) Leitura (csv/tsv/txt/rds) + rownames seguras
# -----------------------------
.read_table_auto <- function(path) {
  if (!file.exists(path)) stopf("Arquivo não encontrado: %s", path)

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }

  read_with_dt <- requireNamespace("data.table", quietly = TRUE)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    if (read_with_dt) {
      return(as.data.frame(data.table::fread(path, data.table = FALSE, check.names = FALSE)))
    }
    return(read.csv(path, check.names = FALSE))
  }
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
    if (read_with_dt) {
      return(as.data.frame(data.table::fread(path, data.table = FALSE, check.names = FALSE)))
    }
    return(read.delim(path, check.names = FALSE))
  }

  stopf("Formato não suportado: %s (use .csv/.tsv/.txt/.rds)", path)
}

.coerce_numeric_df <- function(df) {
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  df
}

.set_rownames_safe <- function(df, row_col = 1) {
  if (is.null(df) || ncol(df) < row_col) return(df)
  row_ids <- df[[row_col]]
  if (any(duplicated(row_ids))) {
    warnf("Encontrados %d IDs duplicados. Mantendo primeira ocorrência.", sum(duplicated(row_ids)))
    keep <- !duplicated(row_ids)
    df <- df[keep, , drop = FALSE]
    row_ids <- row_ids[keep]
  }
  rownames(df) <- row_ids
  df <- df[, -row_col, drop = FALSE]
  df
}

# -----------------------------
# 3) Carregar e preparar dados
#    - Filtro de genes: top N por variância
#    - Detecção opcional AUC/IC50
# -----------------------------
load_and_prepare_data <- function(expression_file, response_file,
                                  metadata_file = NULL,
                                  top_var_genes = 3000,
                                  drop_all_na_drugs = TRUE,
                                  response_transform = c("auto", "auc", "ic50", "none")) {
  response_transform <- match.arg(response_transform)

  .log("INFO", "Carregando dados...")

  expr <- if (is.character(expression_file)) .read_table_auto(expression_file) else expression_file
  resp <- if (is.character(response_file)) .read_table_auto(response_file) else response_file

  expr <- .set_rownames_safe(as.data.frame(expr, check.names = FALSE))
  resp <- .set_rownames_safe(as.data.frame(resp, check.names = FALSE))

  expr <- .coerce_numeric_df(expr)
  resp <- .coerce_numeric_df(resp)

  .log("INFO", sprintf("  • Bruto: %d amostras x %d genes", nrow(expr), ncol(expr)))

  common <- intersect(rownames(expr), rownames(resp))
  if (length(common) < 20) stopf("Poucas amostras em comum: %d", length(common))

  expr <- expr[common, , drop = FALSE]
  resp <- resp[common, , drop = FALSE]

  # Filtro genes por variância (top N)
  gv <- apply(expr, 2, var, na.rm = TRUE)
  gv[!is.finite(gv)] <- 0
  if (!is.null(top_var_genes) && top_var_genes < ncol(expr)) {
    ord <- order(gv, decreasing = TRUE)
    keep <- ord[seq_len(top_var_genes)]
    expr <- expr[, keep, drop = FALSE]
    .log(
      "INFO",
      sprintf(
        "  • Filtro de Variância: Reduzido de %d para %d genes (Top %d)",
        length(gv), ncol(expr), top_var_genes
      )
    )
  } else {
    .log("INFO", "  • Filtro de Variância: não aplicado (top_var_genes NULL ou >= total)")
  }

  # Remover drogas 100% NA
  if (drop_all_na_drugs) {
    all_na <- sapply(resp, function(x) all(is.na(x)))
    if (any(all_na)) {
      .log("WARN", sprintf("  • Removendo %d drogas 100%% NA", sum(all_na)))
      resp <- resp[, !all_na, drop = FALSE]
    }
  }

  # Transformação de resposta
  if (response_transform != "none") {
    y_vals <- as.numeric(as.matrix(resp))
    y_vals <- y_vals[is.finite(y_vals)]

    detected <- "none"
    if (response_transform == "auto") {
      if (length(y_vals) > 0 && max(y_vals) <= 1 && min(y_vals) >= 0) detected <- "auc"
      if (length(y_vals) > 0 && max(y_vals) > 1) detected <- "ic50"
    } else {
      detected <- response_transform
    }

    if (detected == "auc") {
      .log("INFO", "  • Resposta detectada como AUC (0-1). Convertendo para sensibilidade (1 - AUC)")
      resp <- 1 - resp
    }
    if (detected == "ic50") {
      .log("INFO", "  • Resposta detectada como IC50. Aplicando log10(IC50 + 1)")
      resp <- log10(resp + 1)
    }
  }

  meta <- NULL
  if (!is.null(metadata_file)) {
    meta <- if (is.character(metadata_file)) .read_table_auto(metadata_file) else metadata_file
    meta <- .set_rownames_safe(as.data.frame(meta, check.names = FALSE))
    meta <- meta[common, , drop = FALSE]
  }

  list(
    expression = expr,
    response = resp,
    metadata = meta,
    samples = common,
    genes = colnames(expr),
    drugs = colnames(resp)
  )
}

# -----------------------------
# 4) Split estratificado + múltiplas tentativas
# -----------------------------
.train_test_split_stratified <- function(y, train_frac = 0.8, seed = 123, n_bins = 5,
                                         jitter_sd = 0) {
  set.seed(seed)
  ok <- is.finite(y)
  idx <- which(ok)
  y_ok <- y[idx]
  if (jitter_sd > 0) {
    y_ok <- y_ok + stats::rnorm(length(y_ok), mean = 0, sd = jitter_sd)
  }
  if (length(y_ok) < 20) return(NULL)

  probs <- seq(0, 1, length.out = n_bins + 1)
  qs <- unique(quantile(y_ok, probs = probs, na.rm = TRUE))
  if (length(qs) < 3) {
    train_idx <- sample(idx, size = floor(train_frac * length(idx)))
    test_idx <- setdiff(idx, train_idx)
    if (length(test_idx) == 0 || length(train_idx) == 0) return(NULL)
    return(list(train = train_idx, test = test_idx))
  }

  bins <- cut(y_ok, breaks = qs, include.lowest = TRUE, labels = FALSE)
  tr_local <- caret::createDataPartition(bins, p = train_frac, list = FALSE)
  train_idx <- idx[tr_local]
  test_idx <- setdiff(idx, train_idx)
  if (length(test_idx) == 0 || length(train_idx) == 0) return(NULL)
  list(train = train_idx, test = test_idx)
}

.find_good_split <- function(y, train_frac, seed, split_tries, min_test_size, min_sd_y,
                             min_train_size, min_train_sd, n_bins = 5, jitter_sd = 0) {
  for (i in seq_len(split_tries)) {
    sp <- .train_test_split_stratified(
      y, train_frac = train_frac, seed = seed + i, n_bins = n_bins, jitter_sd = jitter_sd
    )
    if (is.null(sp)) next
    ytr <- y[sp$train]
    yte <- y[sp$test]
    if (sum(is.finite(ytr)) < min_train_size) next
    if (sum(is.finite(yte)) < min_test_size) next
    if (sd(ytr[is.finite(ytr)]) < min_train_sd) next
    if (sd(yte[is.finite(yte)]) < min_sd_y) next
    return(sp)
  }
  NULL
}

.scale_train_apply <- function(X_train, X_test) {
  mu <- colMeans(X_train, na.rm = TRUE)
  sdv <- apply(X_train, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  Xtr <- sweep(sweep(X_train, 2, mu, "-"), 2, sdv, "/")
  Xte <- sweep(sweep(X_test, 2, mu, "-"), 2, sdv, "/")

  list(X_train = Xtr, X_test = Xte, center = mu, scale = sdv)
}

# -----------------------------
# 5) Métricas
# -----------------------------
calculate_metrics <- function(y_train, y_pred_train, y_test, y_pred_test) {
  pear_tr <- .safe_cor(y_train, y_pred_train, "pearson")
  pear_te <- .safe_cor(y_test, y_pred_test, "pearson")
  sp_tr <- .safe_cor(y_train, y_pred_train, "spearman")
  sp_te <- .safe_cor(y_test, y_pred_test, "spearman")

  rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
  mae <- function(a, b) mean(abs(a - b), na.rm = TRUE)

  r2 <- function(y, yp) {
    ok <- is.finite(y) & is.finite(yp)
    if (sum(ok) < 3) return(NA_real_)
    ss_total <- sum((y[ok] - mean(y[ok]))^2)
    ss_res <- sum((y[ok] - yp[ok])^2)
    if (ss_total == 0) return(NA_real_)
    1 - ss_res / ss_total
  }

  list(
    pearson_train = pear_tr,
    pearson_test = pear_te,
    spearman_train = sp_tr,
    spearman_test = sp_te,
    rmse_train = rmse(y_train, y_pred_train),
    rmse_test = rmse(y_test, y_pred_test),
    mae_train = mae(y_train, y_pred_train),
    mae_test = mae(y_test, y_pred_test),
    r2_train = r2(y_train, y_pred_train),
    r2_test = r2(y_test, y_pred_test),
    overfitting = pear_tr - pear_te,
    n_train = sum(is.finite(y_train)),
    n_test = sum(is.finite(y_test))
  )
}

# -----------------------------
# 6) Treino por droga (glmnet)
# -----------------------------
train_drug_model <- function(X, y, drug_name,
                             seed = 123,
                             train_frac = 0.8,
                             n_folds = 5,
                             alpha_values = seq(0, 1, 0.1),
                             lambda_mode = c("1se", "min"),
                             min_samples = 30,
                             resistance_quantile = 0.75,
                             min_test_size = 10,
                             min_sd_y = 1e-6,
                             min_train_sd = 1e-6,
                             split_tries = 100,
                             top_k_genes = 200,
                             n_bins = 5,
                             y_round = NULL,
                             jitter_sd = 0,
                             min_pred_sd = 1e-6) {
  lambda_mode <- match.arg(lambda_mode)

  if (!is.null(y_round)) y <- round(y, y_round)
  ok <- is.finite(y)
  n_ok <- sum(ok)
  if (n_ok < min_samples) {
    .log("WARN", sprintf("%s: amostras insuficientes (%d). Pulando.", drug_name, n_ok))
    return(NULL)
  }

  split <- .find_good_split(
    y,
    train_frac = train_frac,
    seed = seed,
    split_tries = split_tries,
    min_test_size = min_test_size,
    min_sd_y = min_sd_y,
    min_train_size = max(5, min_samples - min_test_size),
    min_train_sd = min_train_sd,
    n_bins = n_bins,
    jitter_sd = jitter_sd
  )
  if (is.null(split)) {
    .log("WARN", sprintf("%s: não consegui split válido (teste sem variância/tamanho). Pulando.", drug_name))
    return(NULL)
  }

  X_train <- as.matrix(X[split$train, , drop = FALSE])
  X_test <- as.matrix(X[split$test, , drop = FALSE])
  y_train <- as.numeric(y[split$train])
  y_test <- as.numeric(y[split$test])

  sc <- .scale_train_apply(X_train, X_test)
  X_train <- sc$X_train
  X_test <- sc$X_test

  best <- list(alpha = NA_real_, lambda = NA_real_, score = -Inf)

  for (a in alpha_values) {
    cv <- tryCatch(
      glmnet::cv.glmnet(
        x = X_train, y = y_train,
        family = "gaussian",
        alpha = a,
        nfolds = min(n_folds, sum(is.finite(y_train))),
        type.measure = "mse",
        standardize = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(cv)) next

    lam <- if (lambda_mode == "1se") cv$lambda.1se else cv$lambda.min
    mse <- cv$cvm[which.min(abs(cv$lambda - lam))]
    score <- -mse

    if (is.finite(score) && score > best$score) {
      best <- list(alpha = a, lambda = lam, score = score)
    }
  }

  if (!is.finite(best$alpha)) {
    .log("WARN", sprintf("%s: falha na otimização alpha/lambda.", drug_name))
    return(NULL)
  }

  fit_full <- glmnet::glmnet(
    x = X_train, y = y_train,
    family = "gaussian",
    alpha = best$alpha,
    lambda = best$lambda,
    standardize = FALSE
  )

  # Seleção top_k_genes por |coef|
  coefs <- as.matrix(coef(fit_full))
  coefs <- coefs[setdiff(rownames(coefs), "(Intercept)"), 1, drop = TRUE]
  abscoef <- abs(coefs)
  abscoef[!is.finite(abscoef)] <- 0
  nonzero <- names(abscoef)[abscoef > 0]

  if (length(nonzero) > 0) {
    ord <- order(abscoef[nonzero], decreasing = TRUE)
    keep_genes <- nonzero[ord][seq_len(min(top_k_genes, length(ord)))]
  } else {
    keep_genes <- character(0)
  }

  if (length(keep_genes) < 2) {
    keep_genes <- colnames(X_train)[seq_len(min(top_k_genes, ncol(X_train)))]
  }

  Xtr2 <- X_train[, keep_genes, drop = FALSE]
  Xte2 <- X_test[, keep_genes, drop = FALSE]

  fit <- glmnet::glmnet(
    x = Xtr2, y = y_train,
    family = "gaussian",
    alpha = best$alpha,
    lambda = best$lambda,
    standardize = FALSE
  )

  pred_train <- as.numeric(predict(fit, Xtr2))
  pred_test <- as.numeric(predict(fit, Xte2))

  metrics <- calculate_metrics(y_train, pred_train, y_test, pred_test)
  if (is.na(metrics$pearson_test) && sd(pred_test, na.rm = TRUE) < min_pred_sd) {
    alt_lambda <- if (lambda_mode == "1se") "min" else "1se"
    .log("INFO", sprintf("%s: teste sem variância. Tentando lambda=%s.", drug_name, alt_lambda))
    cv_alt <- tryCatch(
      glmnet::cv.glmnet(
        x = X_train, y = y_train,
        family = "gaussian",
        alpha = best$alpha,
        nfolds = min(n_folds, sum(is.finite(y_train))),
        type.measure = "mse",
        standardize = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(cv_alt)) {
      alt_lam <- if (alt_lambda == "1se") cv_alt$lambda.1se else cv_alt$lambda.min
      fit_alt <- glmnet::glmnet(
        x = Xtr2, y = y_train,
        family = "gaussian",
        alpha = best$alpha,
        lambda = alt_lam,
        standardize = FALSE
      )
      pred_train <- as.numeric(predict(fit_alt, Xtr2))
      pred_test <- as.numeric(predict(fit_alt, Xte2))
      metrics <- calculate_metrics(y_train, pred_train, y_test, pred_test)
      fit <- fit_alt
      best$lambda <- alt_lam
      best$lambda_mode <- alt_lambda
    }
  }

  thr <- as.numeric(stats::quantile(pred_train, probs = resistance_quantile, na.rm = TRUE))

  .log(
    "INFO",
    sprintf(
      "%s | alpha=%.1f r_test=%s sp_test=%s rmse=%.3f top%d=%d",
      drug_name,
      best$alpha,
      ifelse(is.na(metrics$pearson_test), "NA", sprintf("%.3f", metrics$pearson_test)),
      ifelse(is.na(metrics$spearman_test), "NA", sprintf("%.3f", metrics$spearman_test)),
      metrics$rmse_test,
      top_k_genes, length(keep_genes)
    )
  )

  list(
    drug = drug_name,
    model = fit,
    metrics = metrics,
    predictions = list(train = pred_train, test = pred_test),
    actual = list(train = y_train, test = y_test),
    genes = list(important = keep_genes, all = colnames(X)),
    split = split,
    scale_params = list(center = sc$center, scale = sc$scale),
    parameters = list(alpha = best$alpha, lambda = best$lambda, lambda_mode = lambda_mode),
    resistance = list(
      quantile = resistance_quantile,
      threshold = thr,
      interpretation = "score_alto => mais_resistente (assumindo AUC)"
    )
  )
}

# -----------------------------
# 7) Sumário
# -----------------------------
create_summary_dataframe <- function(results) {
  if (length(results) == 0) return(data.frame())
  do.call(rbind, lapply(names(results), function(drug) {
    r <- results[[drug]]
    m <- r$metrics
    p <- r$parameters
    data.frame(
      drug = drug,
      pearson_train = m$pearson_train,
      pearson_test = m$pearson_test,
      spearman_train = m$spearman_train,
      spearman_test = m$spearman_test,
      rmse_test = m$rmse_test,
      mae_test = m$mae_test,
      r2_test = m$r2_test,
      overfitting = m$overfitting,
      n_genes = length(r$genes$important),
      n_train = m$n_train,
      n_test = m$n_test,
      alpha = p$alpha,
      lambda = p$lambda,
      lambda_mode = p$lambda_mode,
      resistance_q = r$resistance$quantile,
      resistance_thr = r$resistance$threshold,
      stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------
# 8) Pipeline multi-drogas (com filtro de drogas + paralelo)
# -----------------------------
run_perception_pipeline <- function(data,
                                    max_drugs = NULL,
                                    parallel = FALSE,
                                    n_cores = NULL,
                                    output_dir = "perception_results",
                                    seed = 123,
                                    train_frac = 0.8,
                                    n_folds = 5,
                                    alpha_values = seq(0, 1, 0.1),
                                    lambda_mode = "1se",
                                    min_samples = 30,
                                    resistance_quantile = 0.75,
                                    min_test_size = 10,
                                    min_sd_y = 1e-6,
                                    min_train_sd = 1e-6,
                                    split_tries = 100,
                                    top_k_genes = 200,
                                    n_bins = 5,
                                    jitter_sd = 0,
                                    min_pred_sd = 1e-6,
                                    # filtro de drogas
                                    min_sd_y_global = 0.02,
                                    min_unique_y = 10,
                                    y_round = 3) {
  .log("INFO", paste0("\n", paste(rep("=", 70), collapse = "")))
  .log("INFO", "EXECUTANDO PIPELINE PERCEPTION 2.7 (PRO) - FILTRO DE DROGAS")
  .log("INFO", paste(rep("=", 70), collapse = ""))

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  X <- data$expression
  Y <- data$response

  if (!is.null(max_drugs)) {
    max_drugs <- min(max_drugs, ncol(Y))
    Y <- Y[, seq_len(max_drugs), drop = FALSE]
    .log("INFO", sprintf("Limitando análise a %d drogas", ncol(Y)))
  }

  # ---- filtro drogas ----
  drugs0 <- colnames(Y)
  keep <- logical(length(drugs0))

  for (i in seq_along(drugs0)) {
    d <- drugs0[i]
    y <- as.numeric(Y[[d]])
    y <- y[is.finite(y)]
    if (!is.null(y_round)) y <- round(y, y_round)

    if (length(y) < min_samples) { keep[i] <- FALSE; next }
    if (sd(y) < min_sd_y_global) { keep[i] <- FALSE; next }
    if (length(unique(y)) < min_unique_y) { keep[i] <- FALSE; next }

    keep[i] <- TRUE
  }

  Y <- Y[, keep, drop = FALSE]
  drugs <- colnames(Y)

  .log("INFO", sprintf("Drogas antes do filtro: %d", length(drugs0)))
  .log("INFO", sprintf("Drogas após  do filtro: %d", length(drugs)))

  if (length(drugs) == 0) {
    stopf("Nenhuma droga passou no filtro. Ajuste min_sd_y_global/min_unique_y/min_samples.")
  }

  # Paralelo
  if (parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE) || !requireNamespace("foreach", quietly = TRUE)) {
      .log("WARN", "doParallel/foreach não instalados. Rodando sequencial.")
      parallel <- FALSE
    }
  }

  start_time <- Sys.time()
  results <- list()

  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    .log("INFO", sprintf("Paralelo ON (%d cores)", n_cores))

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    res_list <- foreach::foreach(
      d = drugs,
      .packages = c("glmnet", "caret"),
      .errorhandling = "pass"
    ) %dopar% {
      y <- as.numeric(Y[[d]])

      train_drug_model(
        X = X, y = y, drug_name = d,
        seed = seed,
        train_frac = train_frac,
        n_folds = n_folds,
        alpha_values = alpha_values,
        lambda_mode = lambda_mode,
        min_samples = min_samples,
        resistance_quantile = resistance_quantile,
        min_test_size = min_test_size,
        min_sd_y = min_sd_y,
        min_train_sd = min_train_sd,
        split_tries = split_tries,
        top_k_genes = top_k_genes,
        n_bins = n_bins,
        y_round = y_round,
        jitter_sd = jitter_sd,
        min_pred_sd = min_pred_sd
      )
    }
    names(res_list) <- drugs
    results <- res_list[!sapply(res_list, is.null)]
  } else {
    for (d in drugs) {
      y <- as.numeric(Y[[d]])

      fit <- train_drug_model(
        X = X, y = y, drug_name = d,
        seed = seed,
        train_frac = train_frac,
        n_folds = n_folds,
        alpha_values = alpha_values,
        lambda_mode = lambda_mode,
        min_samples = min_samples,
        resistance_quantile = resistance_quantile,
        min_test_size = min_test_size,
        min_sd_y = min_sd_y,
        min_train_sd = min_train_sd,
        split_tries = split_tries,
        top_k_genes = top_k_genes,
        n_bins = n_bins,
        y_round = y_round,
        jitter_sd = jitter_sd,
        min_pred_sd = min_pred_sd
      )
      if (!is.null(fit)) results[[d]] <- fit
    }
  }

  saveRDS(results, file.path(output_dir, "all_models.rds"))
  summary_df <- create_summary_dataframe(results)
  write.csv(summary_df, file.path(output_dir, "models_summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  duration <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  .log("INFO", sprintf("✅ Concluído em %.1f min | modelos=%d", duration, length(results)))
  if (nrow(summary_df) > 0) {
    .log("INFO", sprintf("Média Pearson teste (ign. NA): %.3f", mean(summary_df$pearson_test, na.rm = TRUE)))
  }

  list(models = results, summary = summary_df, output_dir = output_dir)
}

# -----------------------------
# 9) Predição em novas amostras
# -----------------------------
predict_new_samples <- function(models, new_expression_data,
                                output_file = "perception_predictions.csv") {
  if (is.null(models) || length(models) == 0) stopf("Nenhum modelo disponível para predição.")

  new_expression_data <- as.data.frame(new_expression_data, check.names = FALSE)
  new_expression_data <- .coerce_numeric_df(new_expression_data)

  required_genes <- names(models[[1]]$scale_params$center)

  missing <- setdiff(required_genes, colnames(new_expression_data))
  if (length(missing) > 0) {
    warnf("%d genes ausentes nos novos dados. Preenchendo com 0.", length(missing))
    for (g in missing) new_expression_data[[g]] <- 0
  }
  new_expression_data <- new_expression_data[, required_genes, drop = FALSE]

  preds <- data.frame(Sample = rownames(new_expression_data))

  for (drug in names(models)) {
    mdl <- models[[drug]]
    Xn <- as.matrix(new_expression_data)

    mu <- mdl$scale_params$center
    sdv <- mdl$scale_params$scale
    Xn <- sweep(sweep(Xn, 2, mu, "-"), 2, sdv, "/")

    keep_genes <- mdl$genes$important
    Xn2 <- Xn[, keep_genes, drop = FALSE]

    score <- as.numeric(predict(mdl$model, Xn2))
    preds[[drug]] <- score

    thr <- mdl$resistance$threshold
    preds[[paste0(drug, "_resistant")]] <- score > thr
  }

  if (!is.null(output_file)) {
    write.csv(preds, output_file, row.names = FALSE)
    .log("INFO", "Predições salvas em: ", output_file)
  }
  preds
}

# -----------------------------
# 10) Wrapper principal
# -----------------------------
perception <- function(expression_file, response_file, metadata_file = NULL,
                       output_dir = "perception_results",
                       max_drugs = NULL,
                       parallel = FALSE,
                       n_cores = NULL,
                       seed = 123,
                       train_frac = 0.8,
                       n_folds = 5,
                       alpha_values = seq(0, 1, 0.1),
                       lambda_mode = "1se",
                       min_samples = 30,
                       resistance_quantile = 0.75,
                       min_test_size = 10,
                       min_sd_y = 1e-6,
                       min_train_sd = 1e-6,
                       split_tries = 100,
                       top_var_genes = 3000,
                       top_k_genes = 200,
                       n_bins = 5,
                       jitter_sd = 0,
                       min_pred_sd = 1e-6,
                       min_sd_y_global = 0.02,
                       min_unique_y = 10,
                       y_round = 3,
                       install_missing = TRUE,
                       response_transform = "auto") {
  setup_perception(install_missing = install_missing)

  data <- load_and_prepare_data(
    expression_file = expression_file,
    response_file = response_file,
    metadata_file = metadata_file,
    top_var_genes = top_var_genes,
    response_transform = response_transform
  )

  run_perception_pipeline(
    data = data,
    max_drugs = max_drugs,
    parallel = parallel,
    n_cores = n_cores,
    output_dir = output_dir,
    seed = seed,
    train_frac = train_frac,
    n_folds = n_folds,
    alpha_values = alpha_values,
    lambda_mode = lambda_mode,
    min_samples = min_samples,
    resistance_quantile = resistance_quantile,
    min_test_size = min_test_size,
    min_sd_y = min_sd_y,
    min_train_sd = min_train_sd,
    split_tries = split_tries,
    top_k_genes = top_k_genes,
    n_bins = n_bins,
    jitter_sd = jitter_sd,
    min_pred_sd = min_pred_sd,
    min_sd_y_global = min_sd_y_global,
    min_unique_y = min_unique_y,
    y_round = y_round
  )
}

# -----------------------------
# 11) Exemplo de uso (comente se não for executar)
# -----------------------------
# results <- perception(
#   expression_file = "dados/expressao.csv",
#   response_file = "dados/resposta.csv",
#   output_dir = "meus_resultados_lung_v2",
#   parallel = FALSE,
#   min_samples = 50,
#   min_test_size = 15,
#   top_k_genes = 200,
#   top_var_genes = 3000,
#   split_tries = 200
# )
