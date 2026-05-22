#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  options(encoding = "UTF-8")
})

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, args_vec) {
  idx <- which(args_vec == flag)
  if (length(idx) == 0 || idx[length(idx)] == length(args_vec)) {
    return(NULL)
  }
  args_vec[idx[length(idx)] + 1]
}

parse_bool_arg <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }
  val <- tolower(trimws(as.character(x[[1]])))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

config_path <- get_arg_value("--config", args)
if (is.null(config_path) || !file.exists(config_path)) {
  stop("Missing --config <analysis_config.json> argument.")
}
run_multinomial_flag <- parse_bool_arg(get_arg_value("--run_multinomial", args), default = FALSE)

required_packages <- c(
  "jsonlite", "readxl", "dplyr", "glmmTMB", "DHARMa",
  "emmeans", "performance", "ggplot2", "broom.mixed"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ", paste(missing_packages, collapse = ", "),
      ". Run install_r_packages.R first."
    )
  )
}

library(jsonlite)
library(readxl)
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(performance)
library(ggplot2)
library(broom.mixed)

message_line <- function(...) {
  cat(paste0(..., "\n"))
}

as_factor_if_needed <- function(df, vars) {
  for (v in vars) {
    if (!is.factor(df[[v]])) {
      df[[v]] <- as.factor(df[[v]])
    }
    df[[v]] <- droplevels(df[[v]])
  }
  df
}

pick_first_existing <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

list_get <- function(lst, key, default = NULL) {
  val <- lst[[key]]
  if (is.null(val) || length(val) == 0) {
    return(default)
  }
  val
}

build_random_candidates <- function(within_vars, animal_id) {
  if (length(within_vars) == 0) {
    return(c(paste0("(1|`", animal_id, "`)") ))
  }

  slope_part <- paste0("`", within_vars, "`", collapse = " + ")
  c(
    paste0("(1 + ", slope_part, " | `", animal_id, "`)") ,
    paste0("(1 + ", slope_part, " || `", animal_id, "`)") ,
    paste0("(1|`", animal_id, "`)")
  )
}

quote_terms <- function(vars) {
  if (length(vars) == 0) return(character(0))
  paste0("`", vars, "`")
}

build_random_candidates_v2 <- function(animal_id, random_effect_vars, random_slope_vars) {
  base_intercept <- paste0("(1|`", animal_id, "`)")

  extra_re <- character(0)
  if (length(random_effect_vars) > 0) {
    extra_re <- paste0("(1|`", random_effect_vars, "`)")
  }

  if (length(random_slope_vars) == 0) {
    return(unique(c(
      paste(c(base_intercept, extra_re), collapse = " + "),
      base_intercept
    )))
  }

  slope_terms <- paste0("`", random_slope_vars, "`", collapse = " + ")
  cands <- c(
    paste(c(paste0("(1 + ", slope_terms, " | `", animal_id, "`)"), extra_re), collapse = " + "),
    paste(c(paste0("(1 + ", slope_terms, " || `", animal_id, "`)"), extra_re), collapse = " + "),
    paste(c(base_intercept, extra_re), collapse = " + "),
    base_intercept
  )
  unique(cands)
}

report_image_markdown <- function(image_path) {
  if (is.null(image_path) || !nzchar(image_path) || !file.exists(image_path)) {
    return(character(0))
  }

  rel_path <- file.path("figures", basename(image_path))
  c("", paste0("![](", rel_path, ")"))
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.4f", p)
}

coerce_numeric_response <- function(x) {
  if (is.numeric(x)) {
    return(list(values = as.numeric(x), converted_bad = 0L))
  }

  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  # Handle French-style decimal commas when present.
  x_chr <- gsub(",", ".", x_chr, fixed = TRUE)

  y_num <- suppressWarnings(as.numeric(x_chr))
  converted_bad <- sum(!is.na(x_chr) & is.na(y_num))

  list(values = y_num, converted_bad = converted_bad)
}

sig_code <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

safe_json <- function(x) {
  toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
}

fit_plan_b <- function(df, response, fixed_formula_text, family_obj, animal_id) {
  formula_obj <- as.formula(paste0("`", response, "` ~ ", fixed_formula_text))
  glm_fit <- glm(formula_obj, data = df, family = family_obj)

  robust_table <- NULL
  robust_note <- "Cluster-robust SE unavailable (missing sandwich/lmtest)."

  if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
    vc <- sandwich::vcovCL(glm_fit, cluster = df[[animal_id]])
    ct <- lmtest::coeftest(glm_fit, vcov. = vc)
    robust_table <- data.frame(
      term = rownames(ct),
      estimate = ct[, 1],
      std.error = ct[, 2],
      statistic = ct[, 3],
      p.value = ct[, 4],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    robust_note <- "Cluster-robust SE applied (sandwich::vcovCL)."
  }

  list(model = glm_fit, robust_table = robust_table, robust_note = robust_note)
}

analyze_single_dv <- function(
  df,
  dv_name,
  dv_type,
  animal_id,
  iv_names,
  iv_types,
  output_dir,
  iv_roles = NULL,
  iv_scopes = NULL,
  observation_cycle = 30,
  non_visible_column = NULL
) {
  warnings <- c()
  status <- "ok"

  if (is.null(iv_roles) || length(iv_roles) != length(iv_names)) {
    iv_roles <- rep("fixed_effect", length(iv_names))
  }
  if (is.null(iv_scopes) || length(iv_scopes) != length(iv_names)) {
    iv_scopes <- ifelse(tolower(iv_types) == "between", "between_subject", "within_subject")
  }

  fixed_vars <- unique(iv_names[iv_roles %in% c("fixed_effect", "random_slope")])
  random_effect_vars <- unique(iv_names[iv_roles %in% c("random_effect", "random_intercept")])
  random_slope_vars <- unique(iv_names[iv_roles == "random_slope"])
  if (length(random_slope_vars) == 0) {
    random_slope_vars <- unique(iv_names[tolower(iv_scopes) == "within_subject"])
  }

  if (length(fixed_vars) == 0) {
    fixed_vars <- iv_names
    warnings <- c(warnings, "No fixed_effect/random_slope IV was selected; all IVs were included as fixed effects.")
  }

  needed_cols <- unique(c(animal_id, fixed_vars, random_effect_vars, random_slope_vars, dv_name, non_visible_column))
  model_df <- df %>% dplyr::select(all_of(needed_cols))
  n_total <- nrow(model_df)

  coerced <- coerce_numeric_response(model_df[[dv_name]])
  model_df[[dv_name]] <- coerced$values

  if (coerced$converted_bad > 0) {
    warnings <- c(
      warnings,
      paste0(
        "Dependent variable '", dv_name, "' had ", coerced$converted_bad,
        " non-numeric value(s) that were converted to NA and removed."
      )
    )
  }

  if (dv_type == "proportion_without_non_visible") {
    if (is.null(non_visible_column) || !nzchar(non_visible_column) || !(non_visible_column %in% names(model_df))) {
      return(list(
        dv = dv_name,
        type = dv_type,
        status = "failed",
        message = "Missing non_visible_column for proportion_without_non_visible.",
        warnings = warnings,
        n_total = n_total,
        n_na = sum(is.na(model_df[[dv_name]]))
      ))
    }

    obs_cycle_num <- suppressWarnings(as.numeric(observation_cycle))
    if (!is.finite(obs_cycle_num) || obs_cycle_num <= 0) {
      obs_cycle_num <- 30
      warnings <- c(warnings, "Invalid observation_cycle in config; defaulting to 30.")
    }

    nv_coerced <- coerce_numeric_response(model_df[[non_visible_column]])
    model_df[[non_visible_column]] <- nv_coerced$values
    if (nv_coerced$converted_bad > 0) {
      warnings <- c(
        warnings,
        paste0(
          "Non-visible column '", non_visible_column, "' had ", nv_coerced$converted_bad,
          " non-numeric value(s) converted to NA."
        )
      )
    }

    n_na <- sum(is.na(model_df[[dv_name]]))
    model_df <- model_df %>% dplyr::filter(!is.na(.data[[dv_name]]))

    factor_vars <- unique(c(animal_id, fixed_vars, random_effect_vars))
    model_df <- as_factor_if_needed(model_df, factor_vars)

    fixed_formula_text <- if (length(fixed_vars) > 0) {
      paste(quote_terms(fixed_vars), collapse = " + ")
    } else {
      "1"
    }
    random_candidates <- build_random_candidates_v2(animal_id, random_effect_vars, random_slope_vars)

    run_binomial_variant <- function(label, subtract_non_visible) {
      local_warnings <- c()
      local_df <- model_df

      local_df$.__success <- local_df[[dv_name]]
      if (subtract_non_visible) {
        local_df$.__trials <- obs_cycle_num - local_df[[non_visible_column]]
      } else {
        local_df$.__trials <- obs_cycle_num
      }
      local_df$.__failure <- local_df$.__trials - local_df$.__success

      invalid <- is.na(local_df$.__success) |
        is.na(local_df$.__trials) |
        local_df$.__success < 0 |
        local_df$.__trials <= 0 |
        local_df$.__success > local_df$.__trials |
        local_df$.__failure < 0

      if (any(invalid)) {
        local_warnings <- c(
          local_warnings,
          paste0(label, ": dropped ", sum(invalid), " row(s) with invalid binomial counts/trials.")
        )
        local_df <- local_df[!invalid, , drop = FALSE]
      }

      if (nrow(local_df) < 10) {
        return(list(
          name = label,
          status = "failed",
          message = "Not enough usable rows (<10) for binomial modeling.",
          n_used = nrow(local_df),
          warnings = local_warnings
        ))
      }

      model <- NULL
      selected_random <- NULL
      fit_errors <- c()
      used_family <- "binomial"
      used_zi <- FALSE

      for (rand in random_candidates) {
        model_formula <- as.formula(paste0("cbind(.__success, .__failure) ~ ", fixed_formula_text, " + ", rand))
        fit_try <- tryCatch({
          glmmTMB(
            formula = model_formula,
            data = local_df,
            family = binomial(link = "logit"),
            ziformula = ~0,
            control = glmmTMBControl(
              optimizer = optim,
              optArgs = list(method = "BFGS"),
              profile = FALSE
            )
          )
        }, error = function(e) {
          fit_errors <<- c(fit_errors, paste0("[", rand, "] ", conditionMessage(e)))
          NULL
        })

        if (!is.null(fit_try)) {
          model <- fit_try
          selected_random <- rand
          break
        }
      }

      plan_b_used <- FALSE
      plan_b_note <- NULL
      robust_table <- NULL

      if (is.null(model)) {
        plan_b_used <- TRUE
        local_warnings <- c(local_warnings, paste0(label, ": mixed model failed; switched to GLM + cluster-robust SE (Plan B)."))
        pb <- tryCatch({
          glm_formula <- as.formula(paste0("cbind(.__success, .__failure) ~ ", fixed_formula_text))
          glm_fit <- glm(glm_formula, data = local_df, family = quasibinomial(link = "logit"))

          robust_table <- NULL
          robust_note <- "Cluster-robust SE unavailable (missing sandwich/lmtest)."
          if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
            vc <- sandwich::vcovCL(glm_fit, cluster = local_df[[animal_id]])
            ct <- lmtest::coeftest(glm_fit, vcov. = vc)
            robust_table <- data.frame(
              term = rownames(ct),
              estimate = ct[, 1],
              std.error = ct[, 2],
              statistic = ct[, 3],
              p.value = ct[, 4],
              stringsAsFactors = FALSE,
              row.names = NULL
            )
            robust_note <- "Cluster-robust SE applied (sandwich::vcovCL)."
          }

          list(model = glm_fit, robust_table = robust_table, robust_note = robust_note)
        }, error = function(e) NULL)

        if (is.null(pb) || is.null(pb$model)) {
          return(list(
            name = label,
            status = "failed",
            message = "Plan B GLM failed for this binomial variant.",
            warnings = local_warnings,
            fit_errors = fit_errors,
            n_used = nrow(local_df)
          ))
        }

        model <- pb$model
        robust_table <- pb$robust_table
        plan_b_note <- pb$robust_note
        selected_random <- "Plan B (GLM)"
        used_family <- "quasi-binomial"
      }

      dh <- tryCatch({
        sim <- DHARMa::simulateResiduals(fittedModel = model, n = 500)
        list(
          sim = sim,
          uniformity = DHARMa::testUniformity(sim),
          dispersion = DHARMa::testDispersion(sim),
          zi = DHARMa::testZeroInflation(sim)
        )
      }, error = function(e) {
        local_warnings <<- c(local_warnings, paste0(label, ": DHARMa diagnostics unavailable: ", conditionMessage(e)))
        NULL
      })

      if (!plan_b_used && !is.null(dh)) {
        disp_p <- tryCatch(dh$dispersion$p.value, error = function(e) NA_real_)
        if (is.finite(disp_p) && disp_p < 0.05) {
          bb_formula <- formula(model)
          bb_refit <- tryCatch({
            glmmTMB(
              formula = bb_formula,
              data = local_df,
              family = betabinomial(link = "logit"),
              ziformula = ~0,
              control = glmmTMBControl(
                optimizer = optim,
                optArgs = list(method = "BFGS"),
                profile = FALSE
              )
            )
          }, error = function(e) NULL)

          if (!is.null(bb_refit)) {
            model <- bb_refit
            used_family <- "beta-binomial"
            local_warnings <- c(local_warnings, paste0(label, ": overdispersion detected; switched from binomial to beta-binomial."))
          } else {
            local_warnings <- c(local_warnings, paste0(label, ": overdispersion detected but beta-binomial refit failed; kept binomial."))
          }
        }

        zi_p <- tryCatch(dh$zi$p.value, error = function(e) NA_real_)
        if (used_family == "binomial" && is.finite(zi_p) && zi_p < 0.05) {
          zi_formula <- formula(model)
          zi_refit <- tryCatch({
            glmmTMB(
              formula = zi_formula,
              data = local_df,
              family = binomial(link = "logit"),
              ziformula = ~1,
              control = glmmTMBControl(
                optimizer = optim,
                optArgs = list(method = "BFGS"),
                profile = FALSE
              )
            )
          }, error = function(e) NULL)
          if (!is.null(zi_refit)) {
            model <- zi_refit
            used_zi <- TRUE
            selected_random <- paste0(selected_random, " + ziformula(~1)")
          }
        }
      }

      fixed_tbl <- NULL
      if (!plan_b_used) {
        fixed_tbl <- tryCatch({
          broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
        }, error = function(e) NULL)
      } else if (!is.null(robust_table)) {
        fixed_tbl <- robust_table
      } else {
        fixed_tbl <- tryCatch({
          broom::tidy(model, conf.int = TRUE)
        }, error = function(e) NULL)
      }

      r2_values <- tryCatch({
        r2 <- performance::r2_nakagawa(model)
        list(marginal = unname(r2$R2_marginal), conditional = unname(r2$R2_conditional))
      }, error = function(e) {
        list(marginal = NA_real_, conditional = NA_real_)
      })

      emmeans_tbl <- NULL
      contrasts_tbl <- NULL
      emmeans_plot_path <- NULL
      emmeans_within_plot_path <- NULL

      between_fixed_ivs <- unique(iv_names[iv_roles %in% c("fixed_effect", "random_slope") & tolower(iv_scopes) == "between_subject"])
      within_fixed_ivs  <- unique(iv_names[iv_roles %in% c("fixed_effect", "random_slope") & tolower(iv_scopes) == "within_subject"])
      between_iv <- if (length(between_fixed_ivs) > 0) between_fixed_ivs[1] else if (length(fixed_vars) > 0) fixed_vars[1] else NA_character_
      within_iv  <- if (length(within_fixed_ivs)  > 0) within_fixed_ivs[1]  else NA_character_
      first_iv   <- between_iv

      if (!is.na(first_iv) && nzchar(first_iv)) {
        emm <- tryCatch(emmeans::emmeans(model, specs = first_iv, type = "response"), error = function(e) NULL)
        if (!is.null(emm)) {
          emmeans_tbl <- tryCatch(as.data.frame(emm), error = function(e) NULL)
          contrasts_tbl <- tryCatch(as.data.frame(contrast(emm, method = "pairwise", adjust = "none")), error = function(e) NULL)
          if (!is.null(contrasts_tbl) && "estimate" %in% names(contrasts_tbl)) {
            contrasts_tbl$effect_size <- exp(contrasts_tbl$estimate)
          }

          if (!is.null(emmeans_tbl) && first_iv %in% names(emmeans_tbl)) {
            emmeans_plot_path <- file.path(output_dir, "figures", paste0("emmeans_", dv_name, "_", label, ".png"))
            estimate_col <- pick_first_existing(c("response", "prob", "emmean"), names(emmeans_tbl))
            lower_col <- pick_first_existing(c("lower.CL", "asymp.LCL", "LCL"), names(emmeans_tbl))
            upper_col <- pick_first_existing(c("upper.CL", "asymp.UCL", "UCL"), names(emmeans_tbl))

            if (!is.null(estimate_col)) {
              p <- ggplot(emmeans_tbl, aes(x = .data[[first_iv]], y = .data[[estimate_col]])) +
                geom_point(size = 2) +
                labs(
                  title = paste0("Estimated Marginal Means: ", dv_name, " (", label, ") \u2014 between-subject"),
                  x = first_iv,
                  y = dv_name
                ) +
                theme_minimal(base_size = 12) +
                coord_cartesian(ylim = c(0, NA))
              if (!is.null(lower_col) && !is.null(upper_col)) {
                p <- p + geom_errorbar(aes(ymin = .data[[lower_col]], ymax = .data[[upper_col]]), width = 0.2)
              }
              tryCatch({
                ggsave(emmeans_plot_path, plot = p, width = 7, height = 4.5, dpi = 150)
              }, error = function(e) {
                emmeans_plot_path <<- NULL
              })
            }
          }
        }
      }

      if (!is.na(within_iv) && nzchar(within_iv)) {
        emm_w <- tryCatch(emmeans::emmeans(model, specs = within_iv, type = "response"), error = function(e) NULL)
        if (!is.null(emm_w)) {
          emmeans_w_tbl <- tryCatch(as.data.frame(emm_w), error = function(e) NULL)
          if (!is.null(emmeans_w_tbl) && within_iv %in% names(emmeans_w_tbl)) {
            emmeans_within_plot_path <- file.path(output_dir, "figures", paste0("emmeans_", dv_name, "_", label, "_within.png"))
            estimate_col_w <- pick_first_existing(c("response", "prob", "emmean"), names(emmeans_w_tbl))
            lower_col_w    <- pick_first_existing(c("lower.CL", "asymp.LCL", "LCL"), names(emmeans_w_tbl))
            upper_col_w    <- pick_first_existing(c("upper.CL", "asymp.UCL", "UCL"), names(emmeans_w_tbl))

            if (!is.null(estimate_col_w)) {
              pw <- ggplot(emmeans_w_tbl, aes(x = .data[[within_iv]], y = .data[[estimate_col_w]])) +
                geom_point(size = 2) +
                labs(
                  title = paste0("Estimated Marginal Means: ", dv_name, " (", label, ") \u2014 within-subject"),
                  x = within_iv,
                  y = dv_name
                ) +
                theme_minimal(base_size = 12) +
                coord_cartesian(ylim = c(0, NA))
              if (!is.null(lower_col_w) && !is.null(upper_col_w)) {
                pw <- pw + geom_errorbar(aes(ymin = .data[[lower_col_w]], ymax = .data[[upper_col_w]]), width = 0.2)
              }
              tryCatch({
                ggsave(emmeans_within_plot_path, plot = pw, width = 7, height = 4.5, dpi = 150)
              }, error = function(e) {
                emmeans_within_plot_path <<- NULL
              })
            }
          }
        }
      }

      qq_plot_path <- file.path(output_dir, "figures", paste0("dharma_qq_", dv_name, "_", label, ".png"))
      if (!is.null(dh)) {
        png(qq_plot_path, width = 1000, height = 700)
        plot(dh$sim)
        dev.off()
      } else {
        qq_plot_path <- NULL
      }

      list(
        name = label,
        status = "ok",
        n_used = nrow(local_df),
        selected_random = selected_random,
        family_used = used_family,
        used_zi = used_zi,
        plan_b_used = plan_b_used,
        plan_b_note = plan_b_note,
        fit_errors = fit_errors,
        warnings = local_warnings,
        fixed_effects = fixed_tbl,
        r2 = r2_values,
        emmeans = emmeans_tbl,
        contrasts = contrasts_tbl,
        effect_metric = "Odds Ratio",
        qq_plot_path = qq_plot_path,
        emmeans_plot_path = emmeans_plot_path,
        emmeans_within_plot_path = emmeans_within_plot_path,
        dharma = if (is.null(dh)) NULL else list(
          uniformity_p = tryCatch(unname(dh$uniformity$p.value), error = function(e) NA_real_),
          dispersion_p = tryCatch(unname(dh$dispersion$p.value), error = function(e) NA_real_),
          zi_p = tryCatch(unname(dh$zi$p.value), error = function(e) NA_real_)
        )
      )
    }

    without_nv <- run_binomial_variant("without_non_visible", subtract_non_visible = TRUE)

    final_status <- if (identical(without_nv$status, "failed")) "failed" else "ok"

    return(list(
      dv = dv_name,
      type = dv_type,
      status = final_status,
      n_total = n_total,
      n_na = n_na,
      n_used = list_get(without_nv, "n_used", 0),
      non_visible_column = non_visible_column,
      observation_cycle = obs_cycle_num,
      selected_random = "See sub-analysis",
      method_text = paste0(
        "Count-based proportions were modeled with binomial GLMMs (logit link). ",
        "The denominator was adjusted to exclude non-visible observations: denominator = ",
        obs_cycle_num, " - Non-visible. ",
        "When DHARMa indicated overdispersion, a beta-binomial refit was attempted. ",
        "Zero-inflation was checked after fitting; ZOIB was not used for this count-trial setup. ",
        "Post-hoc contrasts are FDR-corrected across all DVs (cross-DV Benjamini-Hochberg)."
      ),
      warnings = unique(c(warnings, list_get(without_nv, "warnings", character(0)))),
      analyses = list(
        without_non_visible = without_nv
      )
    ))
  }

  n_na <- sum(is.na(model_df[[dv_name]]))
  model_df <- na.omit(model_df)

  if (nrow(model_df) < 10) {
    return(list(
      dv = dv_name,
      type = dv_type,
      status = "failed",
      message = "Not enough non-NA rows (<10)",
      warnings = c(),
      n_total = n_total,
      n_na = n_na
    ))
  }

  model_df <- as_factor_if_needed(model_df, unique(c(animal_id, fixed_vars, random_effect_vars, random_slope_vars)))

  # Type-specific preprocessing
  family_obj <- NULL
  link_name <- ""
  effect_metric <- ""
  zero_one_detected <- FALSE
  used_zi <- FALSE
  boundary_method_note <- "No boundary adjustment needed (no exact 0/1 detected)."
  rt_log_lmm <- NULL
  response_col <- dv_name

  if (dv_type == "counts") {
    if (any(model_df[[dv_name]] < 0, na.rm = TRUE)) {
      bad_n <- sum(model_df[[dv_name]] < 0, na.rm = TRUE)
      warnings <- c(warnings, paste0("Counts DV has ", bad_n, " negative value(s); those rows were dropped."))
      model_df <- model_df %>% dplyr::filter(.data[[dv_name]] >= 0)
    }

    family_obj <- nbinom2(link = "log")
    link_name <- "log"
    effect_metric <- "IRR"
  } else if (dv_type == "rt") {
    if (any(model_df[[dv_name]] <= 0, na.rm = TRUE)) {
      bad_n <- sum(model_df[[dv_name]] <= 0, na.rm = TRUE)
      warnings <- c(warnings, paste0("RT has ", bad_n, " zero/negative values; those rows were dropped for Gamma model."))
      model_df <- model_df %>% dplyr::filter(.data[[dv_name]] > 0)
    }
    family_obj <- Gamma(link = "log")
    link_name <- "log"
    effect_metric <- "Response Ratio"
  } else if (dv_type == "proportion") {
    y <- model_df[[dv_name]]

    if (max(y, na.rm = TRUE) > 1 && max(y, na.rm = TRUE) <= 100) {
      model_df[[dv_name]] <- y / 100
      warnings <- c(warnings, "Proportion values appeared to be 0-100; automatically divided by 100.")
    }

    # Keep only valid proportion values for beta/binomial-based models.
    bad_prop <- sum(model_df[[dv_name]] < 0 | model_df[[dv_name]] > 1, na.rm = TRUE)
    if (bad_prop > 0) {
      warnings <- c(
        warnings,
        paste0("Proportion DV had ", bad_prop, " value(s) outside [0,1]; those rows were dropped.")
      )
      model_df <- model_df %>% dplyr::filter(.data[[dv_name]] >= 0, .data[[dv_name]] <= 1)
    }

    y2 <- model_df[[dv_name]]
    zero_one_detected <- any(y2 %in% c(0, 1), na.rm = TRUE)
    if (zero_one_detected) {
      warnings <- c(warnings, "Exact 0/1 values detected; applying Smithson & Verkuilen-style interior transform for beta modeling.")
    }

    n <- nrow(model_df)
    model_df[[dv_name]] <- (model_df[[dv_name]] * (n - 1) + 0.5) / n
    if (zero_one_detected) {
      boundary_method_note <- "Exact 0/1 values were transformed to interior values using the Smithson & Verkuilen-style adjustment y* = (y*(n-1)+0.5)/n."
    }

    family_obj <- beta_family(link = "logit")
    link_name <- "logit"
    effect_metric <- "Odds Ratio"
  } else {
    return(list(
      dv = dv_name,
      type = dv_type,
      status = "failed",
      message = paste0("Unsupported dv_type: ", dv_type),
      warnings = warnings,
      n_total = n_total,
      n_na = n_na
    ))
  }

  if (nrow(model_df) < 10) {
    return(list(
      dv = dv_name,
      type = dv_type,
      status = "failed",
      message = "Not enough rows after preprocessing (<10)",
      warnings = warnings,
      n_total = n_total,
      n_na = n_na
    ))
  }

  fixed_formula_text <- if (length(fixed_vars) > 0) paste0("`", fixed_vars, "`", collapse = " + ") else "1"
  random_candidates <- build_random_candidates_v2(animal_id, random_effect_vars, random_slope_vars)

  model <- NULL
  selected_random <- NULL
  fit_errors <- c()

  for (rand in random_candidates) {
    model_formula <- as.formula(paste0("`", response_col, "` ~ ", fixed_formula_text, " + ", rand))

    fit_try <- tryCatch({
      glmmTMB(
        formula = model_formula,
        data = model_df,
        family = family_obj,
        ziformula = if (used_zi) ~1 else ~0,
        control = glmmTMBControl(
          optimizer = optim,
          optArgs = list(method = "BFGS"),
          profile = FALSE
        )
      )
    }, error = function(e) {
      fit_errors <<- c(fit_errors, paste0("[", rand, "] ", conditionMessage(e)))
      NULL
    })

    if (!is.null(fit_try)) {
      model <- fit_try
      selected_random <- rand
      break
    }
  }

  plan_b_used <- FALSE
  plan_b_note <- NULL
  robust_table <- NULL

  if (is.null(model)) {
    warnings <- c(warnings, "Mixed model did not converge with random-effects funnel; switched to GLM + cluster-robust SE (Plan B).")
    plan_b_used <- TRUE

    pb <- tryCatch({
      if (dv_type == "counts") {
        fit_plan_b(model_df, response_col, fixed_formula_text, MASS::negative.binomial(theta = 1), animal_id)
      } else if (dv_type == "rt") {
        fit_plan_b(model_df, response_col, fixed_formula_text, Gamma(link = "log"), animal_id)
      } else {
        fit_plan_b(model_df, response_col, fixed_formula_text, quasibinomial(link = "logit"), animal_id)
      }
    }, error = function(e) {
      NULL
    })

    if (is.null(pb) || is.null(pb$model)) {
      return(list(
        dv = dv_name,
        type = dv_type,
        status = "failed",
        message = "Plan B GLM fitting failed after preprocessing. Check DV type and value range.",
        warnings = warnings,
        n_total = n_total,
        n_na = n_na
      ))
    }

    model <- pb$model
    robust_table <- pb$robust_table
    plan_b_note <- pb$robust_note
    selected_random <- "Plan B (GLM)"
  }

  # DHARMa diagnostics
  dh <- tryCatch({
    sim <- DHARMa::simulateResiduals(fittedModel = model, n = 500)
    list(
      sim = sim,
      uniformity = DHARMa::testUniformity(sim),
      dispersion = DHARMa::testDispersion(sim),
      zi = DHARMa::testZeroInflation(sim)
    )
  }, error = function(e) {
    warnings <<- c(warnings, paste0("DHARMa diagnostics unavailable: ", conditionMessage(e)))
    NULL
  })

  if (!is.null(dh)) {
    if (!is.null(dh$uniformity$p.value) && is.finite(dh$uniformity$p.value) && dh$uniformity$p.value < 0.05) {
      warnings <- c(warnings, paste0("DHARMa KS/uniformity test significant (p=", format_p(dh$uniformity$p.value), ")."))
    }
  }

  # If counts model was not ZI initially, check zero inflation and refit once.
  if (dv_type == "counts" && !plan_b_used && !used_zi && !is.null(dh)) {
    zi_p <- tryCatch(dh$zi$p.value, error = function(e) NA_real_)
    if (is.finite(zi_p) && zi_p < 0.05) {
      warnings <- c(warnings, "Zero-inflation detected by DHARMa; refitting as ZINB.")
      model_formula <- formula(model)
      model_zi <- tryCatch({
        glmmTMB(
          formula = model_formula,
          data = model_df,
          family = nbinom2(link = "log"),
          ziformula = ~1,
          control = glmmTMBControl(
            optimizer = optim,
            optArgs = list(method = "BFGS"),
            profile = FALSE
          )
        )
      }, error = function(e) {
        warnings <<- c(warnings, paste0("ZINB refit failed: ", conditionMessage(e)))
        NULL
      })

      if (!is.null(model_zi)) {
        model <- model_zi
        used_zi <- TRUE
        selected_random <- paste0(selected_random, " + ziformula(~1)")
      }
    }
  }

  # For proportion, only enable ZI beta if DHARMa suggests zero-inflation.
  if (dv_type == "proportion" && !plan_b_used && !used_zi && !is.null(dh)) {
    zi_p <- tryCatch(dh$zi$p.value, error = function(e) NA_real_)
    if (is.finite(zi_p) && zi_p < 0.05) {
      warnings <- c(warnings, "Zero-inflation detected by DHARMa for proportion DV; refitting as ZI-beta (ziformula = ~1).")
      model_formula <- formula(model)
      model_zi <- tryCatch({
        glmmTMB(
          formula = model_formula,
          data = model_df,
          family = beta_family(link = "logit"),
          ziformula = ~1,
          control = glmmTMBControl(
            optimizer = optim,
            optArgs = list(method = "BFGS"),
            profile = FALSE
          )
        )
      }, error = function(e) {
        warnings <<- c(warnings, paste0("ZI-beta refit failed: ", conditionMessage(e)))
        NULL
      })

      if (!is.null(model_zi)) {
        model <- model_zi
        used_zi <- TRUE
        selected_random <- paste0(selected_random, " + ziformula(~1)")
      }
    } else {
      warnings <- c(warnings, "Proportion model kept without ZI component (DHARMa did not indicate zero-inflation).")
    }
  }

  if (dv_type == "rt") {
    rt_df <- model_df
    rt_df$.__log_rt <- log(rt_df[[response_col]])

    rt_model <- NULL
    rt_selected_random <- NULL
    rt_fit_errors <- c()
    rt_warnings <- c()

    for (rand in random_candidates) {
      rt_formula <- as.formula(paste0(".__log_rt ~ ", fixed_formula_text, " + ", rand))
      rt_fit <- tryCatch({
        glmmTMB(
          formula = rt_formula,
          data = rt_df,
          family = gaussian(link = "identity"),
          ziformula = ~0,
          control = glmmTMBControl(
            optimizer = optim,
            optArgs = list(method = "BFGS"),
            profile = FALSE
          )
        )
      }, error = function(e) {
        rt_fit_errors <<- c(rt_fit_errors, paste0("[", rand, "] ", conditionMessage(e)))
        NULL
      })

      if (!is.null(rt_fit)) {
        rt_model <- rt_fit
        rt_selected_random <- rand
        break
      }
    }

    if (is.null(rt_model)) {
      rt_log_lmm <- list(
        status = "failed",
        message = "log(RT) LMM companion analysis failed to converge.",
        selected_random = "N/A",
        fit_errors = rt_fit_errors,
        warnings = rt_warnings
      )
    } else {
      rt_dh <- tryCatch({
        sim <- DHARMa::simulateResiduals(fittedModel = rt_model, n = 500)
        list(
          sim = sim,
          uniformity = DHARMa::testUniformity(sim),
          dispersion = DHARMa::testDispersion(sim),
          zi = DHARMa::testZeroInflation(sim),
          outliers = DHARMa::testOutliers(sim)
        )
      }, error = function(e) {
        rt_warnings <<- c(rt_warnings, paste0("DHARMa diagnostics unavailable for log(RT) LMM: ", conditionMessage(e)))
        NULL
      })

      rt_fixed <- tryCatch({
        broom.mixed::tidy(rt_model, effects = "fixed", conf.int = TRUE)
      }, error = function(e) NULL)

      rt_r2 <- tryCatch({
        r2 <- performance::r2_nakagawa(rt_model)
        list(marginal = unname(r2$R2_marginal), conditional = unname(r2$R2_conditional))
      }, error = function(e) {
        list(marginal = NA_real_, conditional = NA_real_)
      })

      rt_log_lmm <- list(
        status = "ok",
        message = "Companion log(RT) LMM fitted on log-transformed RT values.",
        selected_random = rt_selected_random,
        fit_errors = rt_fit_errors,
        warnings = rt_warnings,
        fixed_effects = rt_fixed,
        r2 = rt_r2,
        dharma = if (is.null(rt_dh)) NULL else list(
          uniformity_p = tryCatch(unname(rt_dh$uniformity$p.value), error = function(e) NA_real_),
          dispersion_p = tryCatch(unname(rt_dh$dispersion$p.value), error = function(e) NA_real_),
          zi_p = tryCatch(unname(rt_dh$zi$p.value), error = function(e) NA_real_),
          outlier_p = tryCatch(unname(rt_dh$outliers$p.value), error = function(e) NA_real_)
        )
      )
    }
  }

  # Fixed effects table
  fixed_tbl <- NULL
  if (!plan_b_used) {
    fixed_tbl <- tryCatch({
      broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
    }, error = function(e) {
      warnings <<- c(warnings, paste0("Could not extract fixed-effects table: ", conditionMessage(e)))
      NULL
    })
  } else if (!is.null(robust_table)) {
    fixed_tbl <- robust_table
    names(fixed_tbl)[names(fixed_tbl) == "std.error"] <- "std.error"
  } else {
    fixed_tbl <- tryCatch({
      broom::tidy(model, conf.int = TRUE)
    }, error = function(e) NULL)
  }

  # R2
  r2_values <- tryCatch({
    r2 <- performance::r2_nakagawa(model)
    list(marginal = unname(r2$R2_marginal), conditional = unname(r2$R2_conditional))
  }, error = function(e) {
    list(marginal = NA_real_, conditional = NA_real_)
  })

  # EMMeans: between-subject IV plot + within-subject IV plot (if present)
  emmeans_tbl <- NULL
  contrasts_tbl <- NULL
  emmeans_plot_path <- NULL
  emmeans_within_plot_path <- NULL

  between_fixed_ivs <- unique(iv_names[iv_roles %in% c("fixed_effect", "random_slope") & tolower(iv_scopes) == "between_subject"])
  within_fixed_ivs  <- unique(iv_names[iv_roles %in% c("fixed_effect", "random_slope") & tolower(iv_scopes) == "within_subject"])
  between_iv <- if (length(between_fixed_ivs) > 0) between_fixed_ivs[1] else if (length(fixed_vars) > 0) fixed_vars[1] else NA_character_
  within_iv  <- if (length(within_fixed_ivs)  > 0) within_fixed_ivs[1]  else NA_character_
  first_iv   <- between_iv

  if (!is.na(first_iv) && nzchar(first_iv)) {
    emm <- tryCatch({
      # Use character spec to avoid formula/backtick edge cases with names like EE/SE.
      emmeans::emmeans(model, specs = first_iv, type = "response")
    }, error = function(e) {
      warnings <<- c(warnings, paste0("EMMeans failed: ", conditionMessage(e)))
      NULL
    })

    if (!is.null(emm)) {
      emmeans_tbl <- tryCatch(as.data.frame(emm), error = function(e) NULL)
      contrasts_tbl <- tryCatch(as.data.frame(contrast(emm, method = "pairwise", adjust = "none")), error = function(e) NULL)

      if (!is.null(contrasts_tbl) && "estimate" %in% names(contrasts_tbl)) {
        contrasts_tbl$effect_size <- exp(contrasts_tbl$estimate)
      }

      if (!is.null(emmeans_tbl) && first_iv %in% names(emmeans_tbl)) {
        emmeans_plot_path <- file.path(output_dir, "figures", paste0("emmeans_", dv_name, ".png"))
        dir.create(dirname(emmeans_plot_path), recursive = TRUE, showWarnings = FALSE)

        estimate_col <- pick_first_existing(c("response", "prob", "emmean"), names(emmeans_tbl))
        lower_col <- pick_first_existing(c("lower.CL", "asymp.LCL", "LCL"), names(emmeans_tbl))
        upper_col <- pick_first_existing(c("upper.CL", "asymp.UCL", "UCL"), names(emmeans_tbl))

        if (!is.null(estimate_col)) {
          p <- ggplot(emmeans_tbl, aes(x = .data[[first_iv]], y = .data[[estimate_col]])) +
            geom_point(size = 2) +
            labs(
              title = paste0("Estimated Marginal Means: ", dv_name, " \u2014 between-subject"),
              x = first_iv,
              y = dv_name
            ) +
            theme_minimal(base_size = 12) +
            coord_cartesian(ylim = c(0, NA))

          if (!is.null(lower_col) && !is.null(upper_col)) {
            p <- p + geom_errorbar(aes(ymin = .data[[lower_col]], ymax = .data[[upper_col]]), width = 0.2)
          } else {
            warnings <- c(warnings, "EMMeans CI columns were not found; plot was generated without error bars.")
          }

          tryCatch({
            ggsave(emmeans_plot_path, plot = p, width = 7, height = 4.5, dpi = 150)
          }, error = function(e) {
            warnings <<- c(warnings, paste0("EMMeans plot save failed: ", conditionMessage(e)))
            emmeans_plot_path <<- NULL
          })
        } else {
          warnings <- c(warnings, "EMMeans estimate column was not found; EMMeans plot skipped.")
          emmeans_plot_path <- NULL
        }
      }
    }
  }

  if (!is.na(within_iv) && nzchar(within_iv)) {
    emm_w <- tryCatch({
      emmeans::emmeans(model, specs = within_iv, type = "response")
    }, error = function(e) {
      warnings <<- c(warnings, paste0("EMMeans within-subject failed: ", conditionMessage(e)))
      NULL
    })

    if (!is.null(emm_w)) {
      emmeans_w_tbl <- tryCatch(as.data.frame(emm_w), error = function(e) NULL)
      if (!is.null(emmeans_w_tbl) && within_iv %in% names(emmeans_w_tbl)) {
        emmeans_within_plot_path <- file.path(output_dir, "figures", paste0("emmeans_", dv_name, "_within.png"))
        dir.create(dirname(emmeans_within_plot_path), recursive = TRUE, showWarnings = FALSE)

        estimate_col_w <- pick_first_existing(c("response", "prob", "emmean"), names(emmeans_w_tbl))
        lower_col_w    <- pick_first_existing(c("lower.CL", "asymp.LCL", "LCL"), names(emmeans_w_tbl))
        upper_col_w    <- pick_first_existing(c("upper.CL", "asymp.UCL", "UCL"), names(emmeans_w_tbl))

        if (!is.null(estimate_col_w)) {
          pw <- ggplot(emmeans_w_tbl, aes(x = .data[[within_iv]], y = .data[[estimate_col_w]])) +
            geom_point(size = 2) +
            labs(
              title = paste0("Estimated Marginal Means: ", dv_name, " \u2014 within-subject"),
              x = within_iv,
              y = dv_name
            ) +
            theme_minimal(base_size = 12) +
            coord_cartesian(ylim = c(0, NA))

          if (!is.null(lower_col_w) && !is.null(upper_col_w)) {
            pw <- pw + geom_errorbar(aes(ymin = .data[[lower_col_w]], ymax = .data[[upper_col_w]]), width = 0.2)
          } else {
            warnings <- c(warnings, "EMMeans within-subject CI columns not found; plot generated without error bars.")
          }

          tryCatch({
            ggsave(emmeans_within_plot_path, plot = pw, width = 7, height = 4.5, dpi = 150)
          }, error = function(e) {
            warnings <<- c(warnings, paste0("EMMeans within-subject plot save failed: ", conditionMessage(e)))
            emmeans_within_plot_path <<- NULL
          })
        } else {
          warnings <- c(warnings, "EMMeans within-subject estimate column not found; plot skipped.")
          emmeans_within_plot_path <- NULL
        }
      }
    }
  }

  # DHARMa figure
  qq_plot_path <- file.path(output_dir, "figures", paste0("dharma_qq_", dv_name, ".png"))
  dir.create(dirname(qq_plot_path), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(dh)) {
    png(qq_plot_path, width = 1000, height = 700)
    plot(dh$sim)
    dev.off()
  } else {
    qq_plot_path <- NULL
  }

  method_text <- paste0(
    "Data for ", dv_name, " were analyzed with a ",
    if (plan_b_used) "generalized linear model" else "generalized linear mixed model",
    " using a ",
    if (dv_type == "counts") "Negative Binomial" else if (dv_type == "rt") "Gamma" else "Beta",
    " distribution and a ", link_name,
    " link. Random-effects strategy followed a funnel approach from maximal to random intercept-only when needed. ",
    "DHARMa diagnostics were used for model validation (uniformity, outliers, residual structure), not as the primary model-selection engine, except when an explicit refit rule is stated. ",
    if (dv_type == "rt") {
      "For RT, a companion log(RT) LMM is also fitted for sensitivity analysis. Gamma models cannot include zero latencies; zero values are removed and should be audited for potential censoring/non-response coding (e.g., timeout values). "
    } else {
      ""
    },
    if (dv_type == "proportion") {
      paste0(boundary_method_note, " ZI-beta is only enabled if DHARMa indicates zero-inflation. ")
    } else {
      ""
    },
    "Post-hoc contrasts are FDR-corrected across all DVs (cross-DV Benjamini-Hochberg)."
  )

  list(
    dv = dv_name,
    type = dv_type,
    status = status,
    n_total = n_total,
    n_used = nrow(model_df),
    n_na = n_na,
    selected_random = selected_random,
    zero_one_detected = zero_one_detected,
    boundary_method_note = boundary_method_note,
    used_zi = used_zi,
    rt_log_lmm = rt_log_lmm,
    plan_b_used = plan_b_used,
    plan_b_note = plan_b_note,
    fit_errors = fit_errors,
    warnings = warnings,
    method_text = method_text,
    fixed_effects = fixed_tbl,
    r2 = r2_values,
    emmeans = emmeans_tbl,
    contrasts = contrasts_tbl,
    effect_metric = effect_metric,
    qq_plot_path = qq_plot_path,
    emmeans_plot_path = emmeans_plot_path,
    emmeans_within_plot_path = emmeans_within_plot_path,
    dharma = if (is.null(dh)) NULL else list(
      uniformity_p = tryCatch(unname(dh$uniformity$p.value), error = function(e) NA_real_),
      dispersion_p = tryCatch(unname(dh$dispersion$p.value), error = function(e) NA_real_),
      zi_p = tryCatch(unname(dh$zi$p.value), error = function(e) NA_real_)
    )
  )
}

analyze_multinomial_glmm <- function(
  df,
  dv_list,
  animal_id,
  iv_names,
  output_dir,
  observation_cycle = 30
) {
  warnings <- character(0)

  selected_dvs <- dv_list[vapply(dv_list, function(x) tolower(list_get(x, "type", "")) == "proportion_without_non_visible", logical(1))]
  if (length(selected_dvs) < 2) {
    return(list(
      status = "skipped",
      message = "Multinomial GLMM requires at least two DVs of type proportion_without_non_visible.",
      warnings = warnings
    ))
  }

  behavior_cols <- vapply(selected_dvs, function(x) list_get(x, "name", ""), character(1))
  behavior_cols <- behavior_cols[behavior_cols != ""]
  if (length(behavior_cols) < 2) {
    return(list(status = "skipped", message = "No valid behavior columns found for multinomial GLMM.", warnings = warnings))
  }

  nv_candidates <- unique(vapply(selected_dvs, function(x) list_get(x, "non_visible_column", ""), character(1)))
  nv_candidates <- nv_candidates[nv_candidates != ""]
  if (length(nv_candidates) == 0) {
    return(list(status = "skipped", message = "Multinomial GLMM needs non_visible_column in selected DVs.", warnings = warnings))
  }
  non_visible_col <- nv_candidates[[1]]
  if (length(nv_candidates) > 1) {
    warnings <- c(warnings, "Multiple non_visible_column values detected across DVs; the first one was used.")
  }

  needed_cols <- unique(c(animal_id, iv_names, behavior_cols, non_visible_col))
  model_df <- df %>% dplyr::select(all_of(needed_cols))
  n_total <- nrow(model_df)

  for (b in behavior_cols) {
    coerced <- coerce_numeric_response(model_df[[b]])
    model_df[[b]] <- coerced$values
    if (coerced$converted_bad > 0) {
      warnings <- c(warnings, paste0("Behavior column '", b, "' had ", coerced$converted_bad, " non-numeric value(s) converted to NA."))
    }
  }

  nv_coerced <- coerce_numeric_response(model_df[[non_visible_col]])
  model_df[[non_visible_col]] <- nv_coerced$values
  if (nv_coerced$converted_bad > 0) {
    warnings <- c(warnings, paste0("Non-visible column '", non_visible_col, "' had ", nv_coerced$converted_bad, " non-numeric value(s) converted to NA."))
  }

  obs_cycle_num <- suppressWarnings(as.numeric(observation_cycle))
  if (!is.finite(obs_cycle_num) || obs_cycle_num <= 0) {
    obs_cycle_num <- 30
    warnings <- c(warnings, "Invalid observation_cycle; defaulting to 30.")
  }

  model_df$trials <- obs_cycle_num - model_df[[non_visible_col]]
  invalid <- is.na(model_df$trials) | model_df$trials <= 0
  for (b in behavior_cols) {
    invalid <- invalid | is.na(model_df[[b]]) | model_df[[b]] < 0 | model_df[[b]] > model_df$trials
  }
  if (any(invalid)) {
    warnings <- c(warnings, paste0("Dropped ", sum(invalid), " row(s) with invalid trials/counts for multinomial GLMM."))
    model_df <- model_df[!invalid, , drop = FALSE]
  }

  if (nrow(model_df) < 10) {
    return(list(
      status = "failed",
      message = "Not enough usable rows (<10) for multinomial GLMM.",
      warnings = warnings,
      n_total = n_total,
      n_used = nrow(model_df)
    ))
  }

  model_df <- as_factor_if_needed(model_df, unique(c(animal_id, iv_names)))

  totals <- vapply(behavior_cols, function(b) sum(model_df[[b]], na.rm = TRUE), numeric(1))
  ref_behavior <- names(sort(totals, decreasing = TRUE))[1]
  behavior_levels <- c(ref_behavior, setdiff(behavior_cols, ref_behavior))

  model_df$.__obs_row <- seq_len(nrow(model_df))
  long_parts <- lapply(behavior_cols, function(b) {
    out <- model_df[, unique(c(animal_id, iv_names, "trials", ".__obs_row")), drop = FALSE]
    out$behavior <- b
    out$count <- model_df[[b]]
    out
  })
  long_df <- do.call(rbind, long_parts)
  long_df$behavior <- factor(long_df$behavior, levels = behavior_levels)
  long_df$obs_id <- factor(long_df$.__obs_row)

  fixed_main <- if (length(iv_names) > 0) paste(quote_terms(iv_names), collapse = " + ") else "1"
  inter_terms <- if (length(iv_names) > 0) paste0("behavior:", quote_terms(iv_names), collapse = " + ") else NULL
  rhs <- c(
    "behavior",
    if (fixed_main != "1") fixed_main else NULL,
    inter_terms,
    "offset(log(trials))",
    paste0("(1|`", animal_id, "`)") ,
    "(1|obs_id)"
  )
  model_formula <- as.formula(paste0("count ~ ", paste(rhs, collapse = " + ")))

  fit_warnings <- character(0)
  fit_errors <- character(0)
  model <- tryCatch({
    withCallingHandlers({
      glmmTMB(
        formula = model_formula,
        data = long_df,
        family = poisson(link = "log"),
        ziformula = ~0,
        control = glmmTMBControl(
          optimizer = optim,
          optArgs = list(method = "BFGS"),
          profile = FALSE
        )
      )
    }, warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    fit_errors <<- c(fit_errors, conditionMessage(e))
    NULL
  })

  if (is.null(model)) {
    return(list(
      status = "failed",
      message = "Multinomial GLMM model fitting failed.",
      warnings = unique(c(warnings, fit_warnings)),
      fit_errors = fit_errors,
      reference_behavior = ref_behavior,
      n_total = n_total,
      n_used = nrow(model_df)
    ))
  }

  if (!is.null(model$sdr) && isFALSE(model$sdr$pdHess)) {
    warnings <- c(warnings, "Model convergence warning: Hessian is not positive definite.")
  }
  if (length(fit_warnings) > 0) {
    warnings <- c(warnings, paste0("Model fit warning: ", unique(fit_warnings)))
  }

  dh <- tryCatch({
    sim <- DHARMa::simulateResiduals(fittedModel = model, n = 500)
    list(
      sim = sim,
      uniformity = DHARMa::testUniformity(sim),
      dispersion = DHARMa::testDispersion(sim),
      zi = DHARMa::testZeroInflation(sim)
    )
  }, error = function(e) {
    warnings <<- c(warnings, paste0("DHARMa diagnostics unavailable for multinomial GLMM: ", conditionMessage(e)))
    NULL
  })

  r2_values <- tryCatch({
    r2 <- performance::r2_nakagawa(model)
    list(marginal = unname(r2$R2_marginal), conditional = unname(r2$R2_conditional))
  }, error = function(e) {
    list(marginal = NA_real_, conditional = NA_real_)
  })

  contrasts_all <- list()
  for (iv in iv_names) {
    if (!(iv %in% names(long_df))) next
    if (!is.factor(long_df[[iv]])) next
    if (nlevels(long_df[[iv]]) < 2) next

    emm <- tryCatch({
      emmeans::emmeans(model, specs = c(iv, "behavior"), type = "link")
    }, error = function(e) NULL)
    if (is.null(emm)) next

    ct <- tryCatch({
      as.data.frame(contrast(emm, method = "pairwise", by = "behavior", adjust = "none"))
    }, error = function(e) NULL)
    if (is.null(ct) || nrow(ct) == 0) next

    ct$predictor <- iv
    if ("estimate" %in% names(ct)) {
      ct$effect_size <- exp(ct$estimate)
    }
    if ("lower.CL" %in% names(ct)) {
      ct$effect_lcl <- exp(ct$lower.CL)
    }
    if ("upper.CL" %in% names(ct)) {
      ct$effect_ucl <- exp(ct$upper.CL)
    }
    contrasts_all[[length(contrasts_all) + 1]] <- ct
  }

  contrasts_tbl <- if (length(contrasts_all) > 0) do.call(rbind, contrasts_all) else NULL
  if (!is.null(contrasts_tbl) && "p.value" %in% names(contrasts_tbl)) {
    contrasts_tbl$p.value_fdr <- p.adjust(contrasts_tbl$p.value, method = "fdr")
  }

  qq_plot_path <- file.path(output_dir, "figures", "dharma_qq_multinomial.png")
  if (!is.null(dh)) {
    png(qq_plot_path, width = 1000, height = 700)
    plot(dh$sim)
    dev.off()
  } else {
    qq_plot_path <- NULL
  }

  list(
    status = "ok",
    model_type = "Multinomial GLMM (Poisson trick)",
    message = "Multinomial GLMM fitted using Poisson-multinomial equivalence.",
    reference_behavior = ref_behavior,
    non_visible_column = non_visible_col,
    observation_cycle = obs_cycle_num,
    n_total = n_total,
    n_used = nrow(model_df),
    n_long = nrow(long_df),
    formula = deparse(model_formula),
    warnings = unique(warnings),
    fit_errors = fit_errors,
    r2 = r2_values,
    contrasts = contrasts_tbl,
    qq_plot_path = qq_plot_path,
    dharma = if (is.null(dh)) NULL else list(
      uniformity_p = tryCatch(unname(dh$uniformity$p.value), error = function(e) NA_real_),
      dispersion_p = tryCatch(unname(dh$dispersion$p.value), error = function(e) NA_real_),
      zi_p = tryCatch(unname(dh$zi$p.value), error = function(e) NA_real_)
    )
  )
}

build_markdown_report <- function(input_cfg, analysis_results, output_dir, timestamp_str, multinomial_result = NULL) {
  report_md <- file.path(output_dir, paste0("report_", timestamp_str, ".md"))

  iv_lines <- vapply(input_cfg$iv, function(x) {
    scope <- list_get(x, "subject_scope", list_get(x, "type", "within"))
    role <- list_get(x, "model_role", "fixed_effect")
    paste0("- ", x$name, " (scope: ", scope, "; role: ", role, ")")
  }, character(1))

  dv_lines <- vapply(input_cfg$dv, function(x) {
    paste0("- ", x$name, " (", x$type, ")")
  }, character(1))

  lines <- c(
    "# LabStats GLMM Report",
    "",
    "## Section 0 - Analysis Parameters",
    "",
    paste0("- Input file: ", input_cfg$input),
    paste0("- Sheet: ", input_cfg$sheet),
    paste0("- Animal ID: ", input_cfg$animal_id),
    paste0("- Run timestamp: ", timestamp_str),
    paste0("- R version: ", R.version.string),
    "",
    "### Independent Variables",
    iv_lines,
    "",
    "### Dependent Variables",
    dv_lines,
    "",
    "## Section 1 - Methods"
  )

  for (res in analysis_results) {
    lines <- c(lines, "", paste0("### ", res$dv), "", res$method_text)
  }

  if (!is.null(multinomial_result)) {
    lines <- c(lines, "", "### Multinomial GLMM (Poisson trick)")
    lines <- c(lines, list_get(multinomial_result, "message", "Multinomial GLMM section."))
    if (!is.null(multinomial_result$reference_behavior)) {
      lines <- c(lines, paste0("Reference behavior category: ", multinomial_result$reference_behavior, "."))
    }
  }

  lines <- c(lines, "", "## Section 2 - Results by DV")

  for (res in analysis_results) {
    res_dv <- list_get(res, "dv", "Unknown DV")
    res_type <- list_get(res, "type", "unknown")
    res_status <- list_get(res, "status", "unknown")
    res_n_used <- list_get(res, "n_used", NA_integer_)
    res_n_total <- list_get(res, "n_total", NA_integer_)
    res_n_na <- list_get(res, "n_na", NA_integer_)
    res_selected_random <- list_get(res, "selected_random", "N/A")

    lines <- c(lines, "", paste0("### ", res_dv, " (", res_type, ")"))
    lines <- c(lines, paste0("- Status: ", res_status))
    lines <- c(lines, paste0("- Used rows: ", res_n_used, " / ", res_n_total, " (NA dropped: ", res_n_na, ")"))
    lines <- c(lines, paste0("- Random structure/model: ", res_selected_random))

    if (!is.null(res$analyses)) {
      if (!is.null(res$non_visible_column)) {
        lines <- c(lines, paste0("- Non-visible column: ", res$non_visible_column))
      }
      if (!is.null(res$observation_cycle)) {
        lines <- c(lines, paste0("- Observation cycle: ", res$observation_cycle))
      }

      for (nm in names(res$analyses)) {
        sub <- res$analyses[[nm]]
        lines <- c(lines, "", paste0("#### Sub-analysis: ", nm))
        lines <- c(lines, paste0("- Status: ", list_get(sub, "status", "unknown")))
        lines <- c(lines, paste0("- Family used: ", list_get(sub, "family_used", "N/A")))
        lines <- c(lines, paste0("- Used rows: ", list_get(sub, "n_used", NA_integer_)))
        lines <- c(lines, paste0("- Random structure/model: ", list_get(sub, "selected_random", "N/A")))

        if (!is.null(sub$r2)) {
          lines <- c(lines, paste0("- Nakagawa R2 marginal: ", signif(sub$r2$marginal, 4)))
          lines <- c(lines, paste0("- Nakagawa R2 conditional: ", signif(sub$r2$conditional, 4)))
        }

        if (!is.null(sub$fixed_effects) && nrow(sub$fixed_effects) > 0) {
          lines <- c(lines, "", "Fixed effects:", "")
          fe <- sub$fixed_effects
          fe_lines <- apply(fe, 1, function(row) {
            term <- ifelse("term" %in% names(fe), as.character(row[["term"]]), "")
            est <- ifelse("estimate" %in% names(fe), as.numeric(row[["estimate"]]), NA_real_)
            pval <- ifelse("p.value" %in% names(fe), as.numeric(row[["p.value"]]), NA_real_)
            paste0("- ", term, ": estimate=", signif(est, 4), ", p=", format_p(pval), " ", sig_code(pval))
          })
          lines <- c(lines, fe_lines)
        }

        if (!is.null(sub$contrasts) && nrow(sub$contrasts) > 0) {
          lines <- c(lines, "", paste0("Post-hoc contrasts (", list_get(sub, "effect_metric", "effect"), ") -- raw p / cross-DV FDR p:"), "")
          ct <- sub$contrasts
          ct_lines <- apply(ct, 1, function(row) {
            contrast_name <- ifelse("contrast" %in% names(ct), as.character(row[["contrast"]]), "contrast")
            pval_raw <- ifelse("p.value"     %in% names(ct), as.numeric(row[["p.value"]]),     NA_real_)
            pval_fdr <- ifelse("p.value_fdr" %in% names(ct), as.numeric(row[["p.value_fdr"]]), pval_raw)
            eff <- ifelse("effect_size" %in% names(ct), as.numeric(row[["effect_size"]]), NA_real_)
            paste0("- ", contrast_name, ": effect=", signif(eff, 4),
                   ", p_raw=", format_p(pval_raw),
                   ", p_fdr=", format_p(pval_fdr), " ", sig_code(pval_fdr))
          })
          lines <- c(lines, ct_lines)
        }

        if (!is.null(sub$qq_plot_path)) {
          lines <- c(lines, paste0("- DHARMa plot: ", sub$qq_plot_path))
          lines <- c(lines, report_image_markdown(sub$qq_plot_path))
        }
        if (!is.null(sub$emmeans_plot_path)) {
          lines <- c(lines, paste0("- EMMeans plot (between-subject): ", sub$emmeans_plot_path))
          lines <- c(lines, report_image_markdown(sub$emmeans_plot_path))
        }
        if (!is.null(sub$emmeans_within_plot_path)) {
          lines <- c(lines, paste0("- EMMeans plot (within-subject): ", sub$emmeans_within_plot_path))
          lines <- c(lines, report_image_markdown(sub$emmeans_within_plot_path))
        }
      }
      next
    }

    if (!is.null(res$r2)) {
      lines <- c(lines, paste0("- Nakagawa R2 marginal: ", signif(res$r2$marginal, 4)))
      lines <- c(lines, paste0("- Nakagawa R2 conditional: ", signif(res$r2$conditional, 4)))
    }

    if (!is.null(res$rt_log_lmm)) {
      rt <- res$rt_log_lmm
      lines <- c(lines, "", "Companion analysis: log(RT) LMM")
      lines <- c(lines, paste0("- Status: ", list_get(rt, "status", "unknown")))
      lines <- c(lines, paste0("- Message: ", list_get(rt, "message", "")))
      lines <- c(lines, paste0("- Random structure/model: ", list_get(rt, "selected_random", "N/A")))

      if (!is.null(rt$r2)) {
        lines <- c(lines, paste0("- log(RT) LMM Nakagawa R2 marginal: ", signif(rt$r2$marginal, 4)))
        lines <- c(lines, paste0("- log(RT) LMM Nakagawa R2 conditional: ", signif(rt$r2$conditional, 4)))
      }

      if (!is.null(rt$fixed_effects) && nrow(rt$fixed_effects) > 0) {
        lines <- c(lines, "", "log(RT) LMM fixed effects:", "")
        fe_rt <- rt$fixed_effects
        fe_rt_lines <- apply(fe_rt, 1, function(row) {
          term <- ifelse("term" %in% names(fe_rt), as.character(row[["term"]]), "")
          est <- ifelse("estimate" %in% names(fe_rt), as.numeric(row[["estimate"]]), NA_real_)
          pval <- ifelse("p.value" %in% names(fe_rt), as.numeric(row[["p.value"]]), NA_real_)
          paste0("- ", term, ": estimate=", signif(est, 4), ", p=", format_p(pval), " ", sig_code(pval))
        })
        lines <- c(lines, fe_rt_lines)
      }
    }

    if (!is.null(res$fixed_effects) && nrow(res$fixed_effects) > 0) {
      lines <- c(lines, "", "Fixed effects:", "")
      fe <- res$fixed_effects
      fe_lines <- apply(fe, 1, function(row) {
        term <- ifelse("term" %in% names(fe), as.character(row[["term"]]), "")
        est <- ifelse("estimate" %in% names(fe), as.numeric(row[["estimate"]]), NA_real_)
        pval <- ifelse("p.value" %in% names(fe), as.numeric(row[["p.value"]]), NA_real_)
        paste0("- ", term, ": estimate=", signif(est, 4), ", p=", format_p(pval), " ", sig_code(pval))
      })
      lines <- c(lines, fe_lines)
    }

    if (!is.null(res$contrasts) && nrow(res$contrasts) > 0) {
      lines <- c(lines, "", paste0("Post-hoc contrasts (", res$effect_metric, ") -- raw p / cross-DV FDR p:"), "")
      ct <- res$contrasts
      ct_lines <- apply(ct, 1, function(row) {
        contrast_name <- ifelse("contrast" %in% names(ct), as.character(row[["contrast"]]), "contrast")
        pval_raw <- ifelse("p.value"     %in% names(ct), as.numeric(row[["p.value"]]),     NA_real_)
        pval_fdr <- ifelse("p.value_fdr" %in% names(ct), as.numeric(row[["p.value_fdr"]]), pval_raw)
        eff <- ifelse("effect_size" %in% names(ct), as.numeric(row[["effect_size"]]), NA_real_)
        paste0("- ", contrast_name, ": effect=", signif(eff, 4),
               ", p_raw=", format_p(pval_raw),
               ", p_fdr=", format_p(pval_fdr), " ", sig_code(pval_fdr))
      })
      lines <- c(lines, ct_lines)
    }

    if (!is.null(res$qq_plot_path)) {
      lines <- c(lines, paste0("- DHARMa plot: ", res$qq_plot_path))
      lines <- c(lines, report_image_markdown(res$qq_plot_path))
    }
    if (!is.null(res$emmeans_plot_path)) {
      lines <- c(lines, paste0("- EMMeans plot (between-subject): ", res$emmeans_plot_path))
      lines <- c(lines, report_image_markdown(res$emmeans_plot_path))
    }
    if (!is.null(res$emmeans_within_plot_path)) {
      lines <- c(lines, paste0("- EMMeans plot (within-subject): ", res$emmeans_within_plot_path))
      lines <- c(lines, report_image_markdown(res$emmeans_within_plot_path))
    }
  }

  if (!is.null(multinomial_result)) {
    lines <- c(lines, "", "## Section 2B - Multinomial GLMM (Poisson trick)", "")
    lines <- c(lines, paste0("- Status: ", list_get(multinomial_result, "status", "unknown")))
    lines <- c(lines, paste0("- Message: ", list_get(multinomial_result, "message", "")))

    if (!is.null(multinomial_result$reference_behavior)) {
      lines <- c(lines, paste0("- Reference behavior category: ", multinomial_result$reference_behavior))
    }
    if (!is.null(multinomial_result$n_used) && !is.null(multinomial_result$n_total)) {
      lines <- c(lines, paste0("- Used rows: ", multinomial_result$n_used, " / ", multinomial_result$n_total))
    }
    if (!is.null(multinomial_result$n_long)) {
      lines <- c(lines, paste0("- Long-format rows: ", multinomial_result$n_long))
    }
    if (!is.null(multinomial_result$formula)) {
      lines <- c(lines, paste0("- Model formula: ", multinomial_result$formula))
    }

    if (!is.null(multinomial_result$r2)) {
      lines <- c(lines, paste0("- Nakagawa R2 marginal: ", signif(multinomial_result$r2$marginal, 4)))
      lines <- c(lines, paste0("- Nakagawa R2 conditional: ", signif(multinomial_result$r2$conditional, 4)))
    }

    if (!is.null(multinomial_result$contrasts) && nrow(multinomial_result$contrasts) > 0) {
      lines <- c(lines, "", "Multinomial post-hoc contrasts (raw p / multinomial-section FDR p):", "")
      ct <- multinomial_result$contrasts
      ct_lines <- apply(ct, 1, function(row) {
        contrast_name <- ifelse("contrast" %in% names(ct), as.character(row[["contrast"]]), "contrast")
        behavior_name <- ifelse("behavior" %in% names(ct), as.character(row[["behavior"]]), "behavior")
        predictor_name <- ifelse("predictor" %in% names(ct), as.character(row[["predictor"]]), "predictor")
        pval_raw <- ifelse("p.value" %in% names(ct), as.numeric(row[["p.value"]]), NA_real_)
        pval_fdr <- ifelse("p.value_fdr" %in% names(ct), as.numeric(row[["p.value_fdr"]]), pval_raw)
        eff <- ifelse("effect_size" %in% names(ct), as.numeric(row[["effect_size"]]), NA_real_)
        lcl <- ifelse("effect_lcl" %in% names(ct), as.numeric(row[["effect_lcl"]]), NA_real_)
        ucl <- ifelse("effect_ucl" %in% names(ct), as.numeric(row[["effect_ucl"]]), NA_real_)
        paste0(
          "- ", predictor_name, " | ", behavior_name, " | ", contrast_name,
          ": OR=", signif(eff, 4), " [", signif(lcl, 4), ", ", signif(ucl, 4), "]",
          ", p_raw=", format_p(pval_raw), ", p_fdr=", format_p(pval_fdr), " ", sig_code(pval_fdr)
        )
      })
      lines <- c(lines, ct_lines)
    }

    if (!is.null(multinomial_result$qq_plot_path)) {
      lines <- c(lines, paste0("- DHARMa plot: ", multinomial_result$qq_plot_path))
      lines <- c(lines, report_image_markdown(multinomial_result$qq_plot_path))
    }
  }

  lines <- c(lines, "", "## Section 3 - Warnings and Diagnostics", "")
  for (res in analysis_results) {
    res_warnings <- list_get(res, "warnings", character(0))
    if (!is.null(res$analyses)) {
      for (nm in names(res$analyses)) {
        sub_warn <- list_get(res$analyses[[nm]], "warnings", character(0))
        if (length(sub_warn) > 0) {
          lines <- c(lines, paste0("### ", list_get(res, "dv", "Unknown DV"), " / ", nm))
          lines <- c(lines, paste0("- WARNING: ", sub_warn))
        }
      }
    }
    if (length(res_warnings) > 0) {
      lines <- c(lines, paste0("### ", list_get(res, "dv", "Unknown DV")))
      lines <- c(lines, paste0("- WARNING: ", res_warnings))
    }

    if (!is.null(res$rt_log_lmm)) {
      rt_warn <- list_get(res$rt_log_lmm, "warnings", character(0))
      rt_err <- list_get(res$rt_log_lmm, "fit_errors", character(0))
      if (length(rt_warn) > 0) {
        lines <- c(lines, paste0("### ", list_get(res, "dv", "Unknown DV"), " / log(RT) LMM"))
        lines <- c(lines, paste0("- WARNING: ", rt_warn))
      }
      if (length(rt_err) > 0) {
        lines <- c(lines, paste0("### ", list_get(res, "dv", "Unknown DV"), " / log(RT) LMM fit errors"))
        lines <- c(lines, paste0("- ERROR: ", rt_err))
      }
    }
  }

  if (!is.null(multinomial_result)) {
    mw <- list_get(multinomial_result, "warnings", character(0))
    if (length(mw) > 0) {
      lines <- c(lines, "### Multinomial GLMM")
      lines <- c(lines, paste0("- WARNING: ", mw))
    }
  }

  lines <- c(lines, "", "## Section 4 - Reporting Example", "")
  for (res in analysis_results) {
    if (!is.null(res$analyses)) {
      lines <- c(
        lines,
        paste0(
          "- For ", list_get(res, "dv", "this dependent variable"),
          ": the binomial model uses a denominator adjusted to exclude non-visible observations",
          " (denominator = cycle - non-visible);",
          " if overdispersion is detected, a beta-binomial model is preferred over ZOIB for count-trial data.",
          " Post-hoc contrasts are FDR-corrected across all DVs (cross-DV Benjamini-Hochberg)."
        )
      )
      next
    }

    res_plan_b_used <- isTRUE(list_get(res, "plan_b_used", FALSE))
    res_dv <- list_get(res, "dv", "this dependent variable")
    res_type <- toupper(list_get(res, "type", "unknown"))
    res_effect_metric <- list_get(res, "effect_metric", "effect size")

    lines <- c(
      lines,
      paste0(
        "- ",
        "For ", res_dv, ": The treatment effect was evaluated with a ",
        if (res_plan_b_used) "GLM" else "GLMM",
        " (", res_type, ") and reported as ", res_effect_metric,
        " with 95% confidence intervals and cross-DV FDR-corrected p-values (Benjamini-Hochberg across all DVs)."
      )
    )
  }

  if (!is.null(multinomial_result)) {
    lines <- c(
      lines,
      paste0(
        "- Multinomial GLMM (Poisson trick): counts were modeled jointly across behaviors with reference category '",
        list_get(multinomial_result, "reference_behavior", "most frequent behavior"),
        "', using denominator = cycle - non-visible and random effects (1|",
        input_cfg$animal_id,
        ") + (1|obs_id)."
      )
    )
  }

  writeLines(lines, report_md, useBytes = TRUE)
  report_md
}

attempt_render <- function(report_md_path, output_dir, timestamp_str) {
  outputs <- list(pdf = NULL, html = NULL)

  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    html_out <- file.path(output_dir, paste0("report_", timestamp_str, ".html"))
    try({
      rmarkdown::pandoc_convert(report_md_path, to = "html", output = html_out)
      outputs$html <- html_out
    }, silent = TRUE)

    pdf_out <- file.path(output_dir, paste0("report_", timestamp_str, ".pdf"))
    try({
      rmarkdown::pandoc_convert(report_md_path, to = "pdf", output = pdf_out)
      if (file.exists(pdf_out)) {
        outputs$pdf <- pdf_out
      }
    }, silent = TRUE)
  }

  outputs
}

write_markdown_pdf_copy <- function(report_md_path, output_dir, timestamp_str) {
  out_pdf <- file.path(output_dir, paste0("report_", timestamp_str, "_mdcopy.pdf"))

  md_lines <- tryCatch(readLines(report_md_path, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  if (length(md_lines) == 0) {
    md_lines <- c("(Markdown report was empty or unreadable.)")
  }

  wrapped <- unlist(strwrap(md_lines, width = 95, simplify = FALSE), use.names = FALSE)
  if (length(wrapped) == 0) {
    wrapped <- c("(No printable markdown content.)")
  }

  lines_per_page <- 55
  chunks <- split(wrapped, ceiling(seq_along(wrapped) / lines_per_page))

  ok <- tryCatch({
    grDevices::pdf(out_pdf, width = 8.5, height = 11, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)

    for (i in seq_along(chunks)) {
      grid::grid.newpage()
      header <- paste0("Markdown Copy - Page ", i, "/", length(chunks))
      body <- paste(chunks[[i]], collapse = "\n")

      grid::grid.text(
        label = header,
        x = grid::unit(0.03, "npc"),
        y = grid::unit(0.98, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 10, fontface = "bold", fontfamily = "mono")
      )

      grid::grid.text(
        label = body,
        x = grid::unit(0.03, "npc"),
        y = grid::unit(0.94, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 8.5, fontfamily = "mono")
      )
    }

    TRUE
  }, error = function(e) FALSE)

  if (ok && file.exists(out_pdf)) out_pdf else NULL
}

cfg <- fromJSON(config_path, simplifyVector = FALSE)
input_file <- cfg$input
sheet_name <- cfg$sheet
animal_id <- cfg$animal_id
iv_list <- cfg$iv
dv_list <- cfg$dv

if (is.null(input_file) || is.null(sheet_name) || is.null(animal_id)) {
  stop("Config must define input, sheet, and animal_id.")
}
if (!file.exists(input_file)) {
  stop(paste0("Input file not found: ", input_file))
}

if (is.null(iv_list) || length(iv_list) == 0) {
  stop("At least one independent variable is required.")
}
if (is.null(dv_list) || length(dv_list) == 0) {
  stop("At least one dependent variable is required.")
}

output_dir <- cfg$output_dir
if (is.null(output_dir) || !nzchar(output_dir)) {
  output_dir <- file.path(getwd(), "outputs", format(Sys.time(), "%Y%m%d_%H%M%S"))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

message_line("Reading Excel file...")
raw_df <- read_excel(input_file, sheet = sheet_name)
raw_df <- as.data.frame(raw_df)

iv_names <- vapply(iv_list, function(x) x$name, character(1))
iv_types <- vapply(iv_list, function(x) tolower(x$type), character(1))
iv_roles <- vapply(iv_list, function(x) tolower(list_get(x, "model_role", "fixed_effect")), character(1))
iv_scopes <- vapply(iv_list, function(x) tolower(list_get(x, "subject_scope", ifelse(tolower(list_get(x, "type", "within")) == "between", "between_subject", "within_subject"))), character(1))
observation_cycle <- list_get(cfg, "observation_cycle", 30)
run_multinomial <- isTRUE(run_multinomial_flag) || isTRUE(parse_bool_arg(list_get(cfg, "run_multinomial", FALSE), default = FALSE))

needed_cols <- unique(c(animal_id, iv_names, vapply(dv_list, function(x) x$name, character(1))))
dv_nonvisible_cols <- vapply(dv_list, function(x) list_get(x, "non_visible_column", ""), character(1))
needed_cols <- unique(c(needed_cols, dv_nonvisible_cols[dv_nonvisible_cols != ""]))
missing_cols <- setdiff(needed_cols, names(raw_df))
if (length(missing_cols) > 0) {
  stop(paste0("Columns not found in selected sheet: ", paste(missing_cols, collapse = ", ")))
}

timestamp_str <- tools::file_path_sans_ext(basename(input_file))
results <- list()

for (dv in dv_list) {
  message_line("Analyzing DV: ", dv$name, " [", dv$type, "]")
  res <- analyze_single_dv(
    df = raw_df,
    dv_name = dv$name,
    dv_type = tolower(dv$type),
    animal_id = animal_id,
    iv_names = iv_names,
    iv_types = iv_types,
    iv_roles = iv_roles,
    iv_scopes = iv_scopes,
    observation_cycle = observation_cycle,
    non_visible_column = list_get(dv, "non_visible_column", NULL),
    output_dir = output_dir
  )
  results[[length(results) + 1]] <- res
}

multinomial_result <- NULL
if (run_multinomial) {
  message_line("Running additive section: Multinomial GLMM (Poisson trick)")
  multinomial_result <- analyze_multinomial_glmm(
    df = raw_df,
    dv_list = dv_list,
    animal_id = animal_id,
    iv_names = iv_names,
    observation_cycle = observation_cycle,
    output_dir = output_dir
  )
}

# ---- Cross-DV FDR correction (Benjamini-Hochberg across all DVs) ------------
# emmeans contrasts were collected with adjust="none" (raw p-values).
# We now pool every pairwise-contrast p-value across all DVs and apply a single
# BH-FDR pass, storing the result in p.value_fdr within each contrasts table.
{
  all_pvals  <- numeric(0)
  pval_index <- list()

  for (i in seq_along(results)) {
    res <- results[[i]]
    if (!is.null(res$analyses)) {
      for (sub_nm in names(res$analyses)) {
        ct <- res$analyses[[sub_nm]]$contrasts
        if (!is.null(ct) && is.data.frame(ct) && "p.value" %in% names(ct) && nrow(ct) > 0) {
          start      <- length(all_pvals) + 1
          all_pvals  <- c(all_pvals, ct$p.value)
          pval_index[[length(pval_index) + 1]] <- list(
            res_idx = i, sub = sub_nm, rows = start:length(all_pvals)
          )
        }
      }
    } else {
      ct <- res$contrasts
      if (!is.null(ct) && is.data.frame(ct) && "p.value" %in% names(ct) && nrow(ct) > 0) {
        start      <- length(all_pvals) + 1
        all_pvals  <- c(all_pvals, ct$p.value)
        pval_index[[length(pval_index) + 1]] <- list(
          res_idx = i, sub = NULL, rows = start:length(all_pvals)
        )
      }
    }
  }

  if (length(all_pvals) > 0) {
    adj_pvals <- p.adjust(all_pvals, method = "fdr")
    for (pm in pval_index) {
      adj <- adj_pvals[pm$rows]
      if (!is.null(pm$sub)) {
        results[[pm$res_idx]]$analyses[[pm$sub]]$contrasts$p.value_fdr <- adj
      } else {
        results[[pm$res_idx]]$contrasts$p.value_fdr <- adj
      }
    }
    message_line(
      "Cross-DV FDR correction applied across ", length(all_pvals),
      " contrast(s) from ", length(pval_index), " DV(s)."
    )
  }
}
# ---- End cross-DV FDR -------------------------------------------------------

report_md <- build_markdown_report(cfg, results, output_dir, timestamp_str, multinomial_result = multinomial_result)
rendered <- attempt_render(report_md, output_dir, timestamp_str)
report_pdf_copy <- write_markdown_pdf_copy(report_md, output_dir, timestamp_str)

summary_obj <- list(
  status = "ok",
  generated_at = as.character(Sys.time()),
  input = list(file = input_file, sheet = sheet_name),
  output_dir = output_dir,
  report_md = report_md,
  report_pdf = rendered$pdf,
  report_html = rendered$html,
  report_pdf_copy = report_pdf_copy,
  results = results,
  multinomial = multinomial_result,
  session = list(
    r_version = R.version.string,
    packages = sapply(required_packages, function(p) as.character(utils::packageVersion(p)), simplify = FALSE)
  )
)

summary_path <- file.path(output_dir, "results_summary.json")
writeLines(safe_json(summary_obj), con = summary_path, useBytes = TRUE)

message_line("Done.")
message_line("Summary JSON: ", summary_path)
if (!is.null(rendered$pdf)) message_line("Report PDF: ", rendered$pdf)
if (!is.null(rendered$html)) message_line("Report HTML: ", rendered$html)
if (!is.null(report_pdf_copy)) message_line("Markdown PDF Copy: ", report_pdf_copy)
message_line("Report Markdown: ", report_md)
