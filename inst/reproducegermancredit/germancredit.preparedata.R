# German Credit Data
data <- read.csv('german_credit.csv')
Y <- data[,"Creditability"]
x.categorical <- c('Account.Balance', 'Payment.Status.of.Previous.Credit', 'Purpose', 'Value.Savings.Stocks',
                   'Length.of.current.employment', 'Sex...Marital.Status', 'Guarantors', 'Most.valuable.available.asset',
                   'Concurrent.Credits', 'Type.of.apartment', 'Occupation', 'Telephone', 'Foreign.Worker')
x.quant <- c('Duration.of.Credit..month.', 'Credit.Amount', 'Instalment.per.cent', 'Duration.in.Current.address',
             'Age..years.', 'No.of.Credits.at.this.Bank', 'No.of.dependents')
for(x in x.categorical){
  data[,x] = as.factor(data[,x])
}


fmla <- paste('~',paste(c(x.quant,x.categorical),collapse ='+'))
X <- model.matrix(formula(fmla), data=data)

save(Y, X, file = "germancredit.RData")
