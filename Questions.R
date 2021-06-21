##============
# Questions
#============

## Form matrices (Unscaled pressure)
S11_true <- ifelse(model_num > 4,P_scale^2,1) * S11
S12_true <- P_scale * S12
S21_true <- P_scale * S21
S22_true <-  ifelse(model_num < 5, P_scale^2,1) *S22








