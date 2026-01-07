#library(bigsnpr) 
library(ppmSuite)
library(tidyverse)
library(ggfortify)
library(gt)
library(patchwork)
library(gtsummary)





meanModel = 1
M=1
similarity_function = 1
draws = 2000
burn = 500
thin = 10


ppmxsummary <- function(mod, df, ...) {
  df$predicted <- apply(mod$fitted.values, 2, mean)
  labelings <- t(mod$Si)
  df$label <- as.data.frame(sapply(1:max(mod$Si), function(x) rowMeans(labelings ==x)),
                            col.names = c("l1", "l2", "l3"))  %>%
             mutate(label = max.col(., ties.method = "first")) %>%
             pull(label) %>%
             as.factor()
  g <- df %>%
    ggplot(aes(!!!ensyms(...))) +
    geom_point() +
    geom_point(aes(y = predicted), color = "red") 
  print(g)
 
  test<- table(df[c("label", "riskGroups", "subj_ancestries")])
  test<-c(rowSums(test[,,1]), rowSums(test[,,1]))
  test <- sum((test - 125)**2)
  return(test )
}

# load data
df <- read_table("temp/simulation.fam", col_names = c("FID", "IID", "...1", "...2", "...3", "Y"))  %>%
    left_join(read_table("temp/simulation.best"), by = c("FID", "IID")) %>%
    left_join(read_table("temp/simulation.covar"), by = c("FID", "IID")) %>%
    mutate(
        Y = (Y - mean(Y, na.rm = T)) / sd(Y, na.rm = T),
        Y_anc = lm(Y ~ subj_ancestries, data= .)$resid,
        PRS_anc = lm(PRS ~ subj_ancestries, data= .)$resid)

conf <- c()
for (i in df$riskGroups) {
   if (i == 0) {
    confi = sample(c(0,1), size = 1, prob = c(0.25, 0.75), replace =T)
   } else {
    confi = sample(c(0,1), size = 1, prob = c(0.75, 0.25), replace = T)
   }
   conf <- append(conf, confi) 
}
df$conf <- conf
df <- df %>%
    mutate(
        Yconf = Y + conf,
        Yconf_conf = lm(Yconf ~ conf, data= .)$resid,
        PRS_conf = lm(Yconf_conf ~ conf, data= .)$resid
        )


m1 = gaussian_ppmx(y = df$Yconf, X = df[c("PRS", "subj_ancestries", "conf")], meanModel = meanModel, similarity_function = similarity_function, draws = draws, burn = burn, thin = thin, M = M)
test1 <- ppmxsummary(m1, df, x ="PRS", y= "Yconf")
ggsave("temp/m1.png", width = 10, height = 10)

m2 = gaussian_ppmx(y = df$Y_anc, X = df$PRS_anc, meanModel = meanModel, similarity_function = similarity_function, draws = draws, burn = burn, thin = thin, M = M)
test2 <- ppmxsummary(m2, df, x ="PRS_anc", y= "Y_anc")
ggsave("temp/m2.png", width = 10, height = 10)


# append to test1 and test2 to columns "a" and "b" of "temp/simu.out"
write.table(data.frame("full" = test1, "proj" = test2), file = "temp/simu.out", append = T, col.names = F, row.names = F)
