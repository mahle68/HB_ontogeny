#density plot for the learning curve in BLS8 talk
#Feb 20, 2024


#example data for intense learning peirod
# dd <- data.frame( time = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), learning = c(0,0,1.5,2,2.5,3.5,4,4,4,4,4,4,4,4,4))
# #dd <- data.frame( time = c(1:20), learning = c(0,1,2,2,3.5,4,4,4,4,4,4,4,4,4))
# dd <- data.frame(time = c(1:7), learning = c(rep(0,2),rep(1,5)))
# 
# ggplot(dd, aes(x=learning)) + 
#   geom_density()
# 
# ggplot(dd, aes(x = time, y = learning)) +
#   stat_smooth(method="glm", color="black", se=FALSE,
#               method.args = list(family=binomial)) +
#   theme_void()

#based on an example on stack overflow
data <-  data.frame("dv1" = rep(c(rnorm(15, mean = 20, sd = 7),
                                  rnorm(15, mean = 40, sd = 7)), times = 3),
                    outcome = rep(c(rbinom(15, 0, prob = .95),
                                    rbinom(15, 1, prob = .95)), times = 3))

data$id <- as.factor(data$id)
data$intervention <- as.factor(data$intervention)
data$area <- as.factor(data$area)
data$outcome <- as.factor(data$outcome)

ggplot(data, aes(x = dv1, y = as.numeric(outcome) - 1, color = factor(area))) +
  stat_smooth(method="glm", color="black", se=FALSE,
              method.args = list(family=binomial), lwd = 3) +
  theme_void()

