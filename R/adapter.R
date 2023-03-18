USE_CPP=TRUE
print(paste0("Loading package with USE_CPP=", USE_CPP))

EM_Estep <- function(Zs, is_masked, X_f, w_f, beta_f, sigma2) {
    if (USE_CPP) {
        return(cpp_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2))
    } else {
        return(R_EM_Estep(Zs, is_masked, X_f, w_f, beta_f, sigma2))
    }
}


masked_moments <- function(zs, pi, mu, sigma) {
    if (USE_CPP) {
        return(cpp_masked_moments(zs, pi, mu, sigma))
    } else {
        return(R_masked_moments(zs, pi, mu, sigma))
    }
}


unmasked_moments <- function(z, pi, mu, sigma) {
    if (USE_CPP) {
        return(cpp_unmasked_moments(z, pi, mu, sigma))
    } else {
        return(R_unmasked_moments(z, pi, mu, sigma))
    }
}


pi_matrix <- function(X_f, w_f) {
    if (USE_CPP) {
        return(cpp_pi_matrix(X_f, w_f))
    } else {
        return(R_pi_matrix(X_f, w_f))
    }
}


gating_update <- function(X_f, tau, w_f, gamma, use_proximal_newton=FALSE) {
    if (USE_CPP) {
        return(cpp_gating_update(X_f, tau, w_f, gamma, use_proximal_newton))
    } else {
        return(R_gating_update(X_f, tau, w_f, gamma, use_proximal_newton))
    }
}


weight_marginal_CD <- function(Y, X_f, tau, wk_f, gammak) {
    if (USE_CPP) {
        return(cpp_weight_marginal_CD(Y, X_f, tau, wk_f, gammak))
    } else {
        return(R_weight_marginal_CD(Y, X_f, tau, wk_f, gammak))
    }
}


obj_gating <- function(Y, X_f, tau, wk_f, gammak) {
    if (USE_CPP) {
        return(cpp_obj_gating(Y, X_f, tau, wk_f, gammak))
    } else {
        return(R_obj_gating(Y, X_f, tau, wk_f, gammak))
    }
}


Fs <- function(X_f, tau, w_f, gamma) {
    if (USE_CPP) {
        return(cpp_Fs(X_f, tau, w_f, gamma))
    } else {
        return(R_Fs(X_f, tau, w_f, gamma))
    }
}

beta_update <- function(X_f, D0, D1, D2, beta_f, sigma2, lambda) {
    if (USE_CPP) {
        return(cpp_beta_update(X_f, D0, D1, D2, beta_f, sigma2, lambda))
    } else {
        return(R_beta_update(X_f, D0, D1, D2, beta_f, sigma2, lambda))
    }
}


beta_marginal_CD <- function(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k) {
    if (USE_CPP) {
        return(cpp_beta_marginal_CD(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k))
    } else {
        return(R_beta_marginal_CD(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k))
    }
}


obj_expert <- function(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k) {
    if (USE_CPP) {
        return(cpp_obj_expert(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k))
    } else {
        return(R_obj_expert(X_f, D0k, D1k, D2k, betak_f, sigma2_k, lambda_k))
    }
}


sigma2_update <- function(X_f, D0, D1, D2, beta_f) {
    if (USE_CPP) {
        return(cpp_sigma2_update(X_f, D0, D1, D2, beta_f))
    } else {
        return(R_sigma2_update(X_f, D0, D1, D2, beta_f))
    }
}

SoTh <- function(x, lambda) {
    if (USE_CPP) {
        return(cpp_SoTh(x, lambda))
    } else {
        return(R_SoTh(x, lambda))
    }
}
