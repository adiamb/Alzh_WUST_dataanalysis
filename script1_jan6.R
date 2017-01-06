require(readr)
al_wust = read_csv('~/Desktop/Copy of t1503_mignot_1216.csv')
##remove unneseccary columns
al_wust=al_wust[-c(3, 4, 11:17)]
require(dplyr)
apply(al_wust, 2, FUN = function(x)(unique(x)%>% length()))
al_wust$APOE = as.factor(al_wust$APOE)
al_wust$ab42_stat = as.factor(al_wust$ab42_stat)
al_wust$ratiotau_ab42 = al_wust$inno_Tau/al_wust$inno_Ab42
glm(CDR ~., data = al_wust[-c(1, 9)]) %>% summary()
require(nnet)
neufit=multinom(CDR ~., data = al_wust[-c(1, 9)])
require(randomForest)
al_wust2 = al_wust
al_wust2$CDR=as.factor(al_wust2$CDR)
al_wust2$Gender = as.factor(al_wust2$Gender)
rf_model1=randomForest(CDR ~., data = al_wust2[-c(1, 9)])
varImpPlot(rf_model1)
require(rpart)
require(rattle)
r_model1=rpart(CDR ~., data = al_wust2[-c(1, 9)], method = "class", control = rpart.control(minsplit=2, cp=0))
fancyRpartPlot(r_model1)

##plots
require(reshape2)
al_melt=melt(al_wust, id.vars = c("id", "CDR"), measure.vars = c("inno_Ab42","inno_Tau", "inno_pTau", "ratiotau_ab42", "bmi", "age_at_LP"))
ggplot(al_melt, aes(factor(CDR), factor(variable), fill = log(value)))+geom_raster()+scale_fill_gradient(low = "yellow", high = "red")
ggplot(al_melt, aes(factor(CDR), value))+geom_boxplot()+facet_wrap(~variable, scales = "free")+theme(strip.text.x = element_text(size = 15, colour = "red"), axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 15, face = "bold"))

ggplot(al_wust2, aes(factor(APOE)))+geom_bar()+facet_wrap(~CDR)+theme(strip.text.x = element_text(size = 15, colour = "red"), axis.text.x = element_text(size = 15),  axis.text.y = element_text(size = 15, face = "bold"))


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
