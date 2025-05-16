def fun_critical_state_line(mean_stress, e_csref, lambda_cs, xi_cs):
    return e_csref - lambda_cs * (mean_stress / 100) ** xi_cs