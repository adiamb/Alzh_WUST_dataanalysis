require(readr)
require(dplyr)
require(plyr)
al_wust = read_csv('~/Documents/Alzheimers_WUST_analysis_Jan6/Copy of t1503_mignot_1216.csv')
##read in the HCRT values
hcrt_vals = read_csv('~/Documents/Alzheimers_WUST_analysis_Jan6/Copy of Alzheimer CSF HCRT result-9-17-15.csv') %>% select(MAP_ID, `HCRT(pg/mL)`)
##merge the HCRT values into the main sheet 
intersect(hcrt_vals$MAP_ID, al_wust$id) %>% length() ##check if all ids are matching then merge the data frame
AL_HCRT = merge.data.frame(hcrt_vals, al_wust, by.x = "MAP_ID", by.y = "id")
colnames(AL_HCRT)[2] = c("HCRT_vals")
##remove unneseccary columns
AL_HCRT=AL_HCRT[-c(5, 12:18)]
require(dplyr)
apply(AL_HCRT, 2, FUN = function(x)(unique(x)%>% length()))
AL_HCRT$APOE = as.factor(AL_HCRT$APOE)
AL_HCRT$ab42_stat = as.factor(AL_HCRT$ab42_stat)
AL_HCRT$ratiotau_ab42 = AL_HCRT$inno_Tau/AL_HCRT$inno_Ab42
glm(HCRT_vals~., data =AL_HCRT[-c(1, 9, 11)]) %>% summary()
require(nnet)
neufit=multinom(CDR ~., data = AL_HCRT[-c(1, 9, 11)])
require(randomForest)
AL_HCRT2 = AL_HCRT
AL_HCRT2$CDR=as.factor(AL_HCRT2$CDR)
AL_HCRT2$Gender = as.factor(AL_HCRT2$Gender)
rf_model1=randomForest(HCRT_vals ~inno_pTau+inno_Tau+inno_Ab42, data = AL_HCRT2[-c(1, 4, 11)])
varImpPlot(rf_model1)
glm(HCRT_vals~inno_pTau+inno_Tau+inno_Ab42, data =AL_HCRT2[-c(1, 4, 11)]) %>% summary()

require(rpart)
require(rattle)
r_model1=rpart(HCRT_vals ~., data = AL_HCRT2[-c(1, 4, 11)], control = rpart.control(minsplit=2, cp=0))
fancyRpartPlot(r_model1)

##plots
require(reshape2)
al_melt=melt(AL_HCRT2, id.vars = c("MAP_ID", "CDR"), measure.vars = c("inno_Ab42","inno_Tau", "inno_pTau", "ratiotau_ab42", "bmi", "age_at_LP", "HCRT_vals"))
require(ggplot2)
ggplot(al_melt, aes(factor(CDR), factor(variable), fill = log(value)))+geom_raster()+scale_fill_gradient(low = "yellow", high = "red")
ggplot(al_melt, aes(factor(CDR), value))+geom_boxplot()+facet_wrap(~variable, scales = "free")+theme(strip.text.x = element_text(size = 15, colour = "red"), axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 15, face = "bold"))

corr_eqn <- function(x,y) {
  corr_coef <- cor.test(x, y)
  paste("p.value =", round(corr_coef$p.value, digits = 2), "estimate = ", round(corr_coef$estimate, digits = 2))
}

ggplot(AL_HCRT2, aes(HCRT_vals, inno_pTau))+geom_point(size = 5)+geom_smooth(method = lm, se =F)+annotate(x = 100, y = 50, colour = "red", size = 5, "text", label = corr_eqn(AL_HCRT2$HCRT_vals, AL_HCRT2$inno_pTau))
ggplot(AL_HCRT2, aes(HCRT_vals, inno_Tau))+geom_point(size = 5)+geom_smooth(method = lm, se =F)
ggplot(AL_HCRT2, aes(HCRT_vals, inno_Ab42))+geom_point(size = 5)+geom_smooth(method = lm, se =F)


ggplot(AL_HCRT2, aes(inno_pTau,HCRT_vals, color = CDR))+geom_point(size = 5)+geom_smooth(method = lm, se =F)+facet_grid(~CDR, scales = "free")
ggplot(AL_HCRT2, aes(HCRT_vals, inno_Tau, color = CDR))+geom_point(size = 5)+geom_smooth(method = lm, se =F)+facet_grid(~CDR, scales = "free")
ggplot(AL_HCRT2, aes(HCRT_vals, inno_Ab42, color = CDR))+geom_point(size = 5)+geom_smooth(method = lm, se=F)+facet_grid(~CDR, scales = "free")


ggplot(filter(AL_HCRT2, CDR == 0.5), aes(HCRT_vals, inno_pTau, color = CDR))+geom_point(size = 5)


filter(AL_HCRT2, CDR == 0.5) %>% cor.test(HCRT_vals, inno_pTau)

for (i in unique(AL_HCRT2$CDR)){
  temp=filter(AL_HCRT2, CDR == i)
  temp2 = cor.test(temp$HCRT_vals, temp$inno_Tau)
  temp3 = glm(temp$HCRT_vals ~ temp$inno_Tau)
  print(paste("this is the correlation", temp2$p.value))
  print(summary(temp3))
}


for (i in colnames(AL_HCRT2)[8:10]){
  temp4 = cor.test(AL_HCRT2$HCRT_vals, AL_HCRT2[[i]])
  print(paste(i, temp4$estimate))
}






importance    <- importance(rf_model1)
varImportance <- data.frame(Variables = row.names(importance), 
                            Importance = round(importance[ ,'MeanDecreaseGini'],2))

# Create a rank variable based on importance
rankImportance <- varImportance %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))

# Use ggplot2 to visualize the relative importance of variables
ggplot(rankImportance, aes(x = reorder(Variables, Importance), 
                           y = Importance, fill = Importance)) +
  geom_bar(stat='identity') + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),
            hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables') +
  coord_flip()+theme(axis.text.y = element_text(size = 12, face = "bold"))
