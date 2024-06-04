#Change input/output directory as needed. 
polmuts <- read.csv('~/all_3D_consensus_mutations.csv')

#equation that relates the coefficients of the GLM and y-value to corresponding x-value
getInt <- function(p,m,b){
  (log(p/(1-p)) - b/m)
}

fixdf <- data.frame(matrix(nrow=0, ncol=10))
colnames(fixdf) <- c('cell', 'temperature', 'strain', 'aaSub', 'passLow', 'passLowCIlow', 'passLowCIhigh', 'passHigh', 'passHighCIlow', 'passHighCIhigh')
#loop through all three variables of all conditions and subset the mutants for those variables
for (type in unique(polmuts$cell)){
  for (temp in unique(polmuts$temperature)){
    for (virus in unique(polmuts$strain)){
      subdf <- subset(polmuts, cell==type & temperature==temp & strain==virus)
      bindf <- data.frame(matrix(nrow=0, ncol=3))
      colnames(bindf) <- c('aaSub', 'pass', 'present')
      #make a df that has how many replicates have a given mutation in each passage
      for (mut in unique(subdf[,12])){
        for (pass in unique(subdf[,11])){
          subdf2 <- subset(subdf, mutation==mut & passage==pass)
          bins <- c(rep(1, sum(subdf2$mutation==mut)),rep(0, (8-sum(subdf2$mutation==mut)))) #counts presence of each mutation across 8 reps
          newdf <- data.frame(mut=rep(mut, 8), pass=rep(pass, 8), present=bins)
          bindf <- rbind(bindf, newdf)
        }
        newdf <- data.frame(mut=rep(mut, 8), pass=rep(0, 8), present=rep(0, 8))
        bindf <- rbind(bindf, newdf)
      }
      #for each mutation generate a glm across passage and then calculate passage thresholds for passing probability > 0.5 and > 0.75 with CIs for each
      for (subst in unique(bindf[,1])){
        subdf <- subset(bindf, mut==subst)
        model <- glm(subdf$present~subdf$pass, family='binomial')
        b <- coefficients(model)[1]
        m <- coefficients(model)[2]
        cilowB <- b-summary(model)$coefficients[,2][1]
        cilowM <- m-summary(model)$coefficients[,2][2]
        cihighB <- b+summary(model)$coefficients[,2][1]
        cihighM <- m+summary(model)$coefficients[,2][2]
        passA <- getInt(0.5, m, b)
        passALow <- getInt(0.5, cilowM, cilowB)
        passAHigh <- getInt(0.5, cihighM, cihighB)
        passB <- getInt(0.75, m, b)
        passBLow <- getInt(0.75, cilowM, cilowB)
        passBHigh <- getInt(0.75, cihighM, cihighB)
        newdf <- data.frame(cell=type, temperature=temp, strain=virus, aaSub=subst, passLow=passA, passLowCIlow=passALow, passLowCIhigh=passAHigh, passHigh=passB, passHighCIlow=passBLow, passHighCIhigh=passBHigh)
        fixdf <- rbind(fixdf, newdf)
        #plot the presence of each mutation over passages
        plot <- ggplot(subdf, aes(x=pass, y=present))+geom_jitter(width=0.2, height=0)+labs(title=paste(type, temp, virus, ':', subst)) + stat_smooth(method="glm", color="green", se=TRUE, method.args = list(family=binomial))
        ggsave(paste('~/PassGLM_plots/', type, temp, virus, subst, '.pdf', sep=''), plot)
      }
    }
  }
}

write.csv(fixdf, '~/passageMutGLM_v3_StdEr_Zeros_noFilt.csv', row.names = FALSE, quote=FALSE)