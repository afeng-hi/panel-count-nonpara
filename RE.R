
library(splines2)



len1 <- INN1 + nknots1 - 1
len2 <- INN2 + nknots2 - 1
total_len <- len1 + len2 + length(alpha)

id.alpha <- (len1 + len2 + 1):total_len



id <- which(
  !is.na(theta.est[, id.alpha[1]]) &
    abs(theta.est[, id.alpha[1]] - alpha[1]) < 3 &
    !is.na(theta.est[, id.alpha[2]]) &
    abs(theta.est[, id.alpha[2]] - alpha[2]) < 3
)

cat(
  "Valid simulations after filtering:",
  length(id),
  "out of",
  nrow(theta.est),
  "\n"
)

if (length(id) == 0) {
  stop("No valid simulation replications remain after filtering.")
}


num_grid <- 1000


z_grid <- seq(0, 1, length.out = num_grid)
dz <- z_grid[2] - z_grid[1]


beta_true <- 2 * z_grid^2


BB2_grid <- bSpline(
  z_grid,
  df = INN2,
  knots = knots2,
  degree = INN2 - 1,
  intercept = FALSE,
  Boundary.knots = c(0, 1)
)



s_grid <- seq(0, tau, length.out = num_grid)
ds <- s_grid[2] - s_grid[1]

# 真实 baseline reversed mean function: Lambda_0(s) = s
Lambda_true <- s_grid

# Lambda(s) 的 B-spline 基函数矩阵
BB1_grid <- bSpline(
  s_grid,
  df = INN1,
  knots = knots1,
  degree = INN1 - 1,
  intercept = FALSE,
  Boundary.knots = c(0, tau)
)


# =====================================================================
#计算 beta(x) 的估计曲线
# =====================================================================

# beta(x) 的样条系数索引
id_xi2_start <- len1 + 1
id_xi2_end <- len1 + len2

# 仅使用通过过滤的模拟结果
xi2_est_matrix <- theta.est[
  id,
  id_xi2_start:id_xi2_end,
  drop = FALSE
]

# 每一行对应一次模拟得到的 beta(x) 估计曲线
beta_hat_matrix <- xi2_est_matrix %*% t(BB2_grid)

# 将真实 beta(x) 复制为矩阵
beta_true_matrix <- matrix(
  beta_true,
  nrow = nrow(beta_hat_matrix),
  ncol = num_grid,
  byrow = TRUE
)


# =====================================================================
# 计算 beta(x) 的 ISE、MISE 和相对 MISE
# =====================================================================

# 每次模拟在各网格点上的平方误差
squared_diff_beta <- (
  beta_hat_matrix - beta_true_matrix
)^2

# 每次模拟的 integrated squared error
ISE_beta_vec <- rowSums(squared_diff_beta) * dz

# 普通 MISE
MISE_beta <- mean(ISE_beta_vec, na.rm = TRUE)

# ISE 在模拟重复间的标准差
SD_ISE_beta <- sd(ISE_beta_vec, na.rm = TRUE)

# 真实 beta 函数平方的积分
denom_beta <- sum(beta_true^2) * dz

if (!is.finite(denom_beta) || denom_beta <= 0) {
  stop("The integrated squared magnitude of beta_true must be positive.")
}

# 每次模拟的 relative integrated squared error
RE_beta_vec <- ISE_beta_vec / denom_beta

# 平均 RE，即 beta(x) 的相对 MISE
ReMISE_beta <- mean(RE_beta_vec, na.rm = TRUE)

# RE 在模拟重复间的标准差
SD_RE_beta <- sd(RE_beta_vec, na.rm = TRUE)


# =====================================================================
#计算 Lambda(s) 的估计曲线
# =====================================================================

# Lambda(s) 的样条系数索引
id_xi1_start <- 1
id_xi1_end <- len1

# 仅使用通过过滤的模拟结果
xi1_est_matrix <- theta.est[
  id,
  id_xi1_start:id_xi1_end,
  drop = FALSE
]


Lambda_hat_matrix <- xi1_est_matrix %*% t(BB1_grid)


Lambda_true_matrix <- matrix(
  Lambda_true,
  nrow = nrow(Lambda_hat_matrix),
  ncol = num_grid,
  byrow = TRUE
)


# =====================================================================
# 计算 Lambda(s) 的 ISE、MISE 和相对 MISE
# =====================================================================

# 每次模拟在各网格点上的平方误差
squared_diff_Lambda <- (
  Lambda_hat_matrix - Lambda_true_matrix
)^2

# 每次模拟的 integrated squared error
ISE_Lambda_vec <- rowSums(squared_diff_Lambda) * ds

# 普通 MISE
MISE_Lambda <- mean(ISE_Lambda_vec, na.rm = TRUE)

# ISE 在模拟重复间的标准差
SD_ISE_Lambda <- sd(ISE_Lambda_vec, na.rm = TRUE)

# 真实 Lambda 函数平方的积分
denom_Lambda <- sum(Lambda_true^2) * ds

if (!is.finite(denom_Lambda) || denom_Lambda <= 0) {
  stop("The integrated squared magnitude of Lambda_true must be positive.")
}

# 每次模拟的 relative integrated squared error
RE_Lambda_vec <- ISE_Lambda_vec / denom_Lambda

# 平均 RE，即 Lambda(s) 的相对 MISE
ReMISE_Lambda <- mean(RE_Lambda_vec, na.rm = TRUE)

# RE 在模拟重复间的标准差
SD_RE_Lambda <- sd(RE_Lambda_vec, na.rm = TRUE)


# =====================================================================
# 9. 输出最终结果
# =====================================================================
cat("\n=======================================================\n")
cat("Sample Size (n)                       :", n, "\n")
cat("Total simulations                    :", nrow(theta.est), "\n")
cat("Valid simulations after filtering     :", length(id), "\n")
cat("-------------------------------------------------------\n")

cat(">>> Nonparametric Function beta(x) <<<\n")
cat("Valid RE values                       :",
    sum(is.finite(RE_beta_vec)), "\n")
cat("MISE                                  :",
    round(MISE_beta, 6), "\n")
cat("SD of ISE                             :",
    round(SD_ISE_beta, 6), "\n")
cat("Average RE (ReMISE)                   :",
    round(ReMISE_beta, 6), "\n")
cat("SD of RE                              :",
    round(SD_RE_beta, 6), "\n")

cat("-------------------------------------------------------\n")

cat(">>> Baseline Reversed Mean Lambda(s) <<<\n")
cat("Valid RE values                       :",
    sum(is.finite(RE_Lambda_vec)), "\n")
cat("MISE                                  :",
    round(MISE_Lambda, 6), "\n")
cat("SD of ISE                             :",
    round(SD_ISE_Lambda, 6), "\n")
cat("Average RE (ReMISE)                   :",
    round(ReMISE_Lambda, 6), "\n")
cat("SD of RE                              :",
    round(SD_RE_Lambda, 6), "\n")

cat("=======================================================\n")

