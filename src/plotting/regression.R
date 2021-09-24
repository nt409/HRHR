# library(dplyr)
# setwd("C:/Users/user/Documents/Python/PhD/HRHR/src")


df = read.csv("./param_scan/outputs/regression/r_data.csv")

# summary(df)

df$lowHighEqual = df$diff==0

delta_filt = 0
omega_filt = 0.8
theta_filt = 5
sr_filt = 0

df = df %>%
    filter(sr_prop>sr_filt) %>%
    filter(delta_1>delta_filt, delta_2>delta_filt) %>%
    filter(omega_1>omega_filt, omega_2>omega_filt) %>%
    filter(theta_1>theta_filt, theta_2>theta_filt)
    # filter(diff==0)

# head(df)
# summary(df)

# model = lm(diff ~ RRlower + omega_1 + omega_2 + theta_1 + theta_2 + delta_1 + delta_2 + sr_prop + RR + RS + SR, data=df)
# summary(model)

# model = glm(lowHighEqual ~ RRlower + omega_1 + omega_2 + theta_1 + theta_2 + delta_1 + delta_2 + sr_prop + RR + RS + SR, data=df, family="binomial")
# summary(model)

plot(lowHighEqual~sr_prop, data=df)
plot(diff~sr_prop, data=df)
